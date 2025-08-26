suppressPackageStartupMessages({
  library(dplyr)
  library(mice)
})


ckd_epi_formula <- function(version = c("2009", "2021")) {
  version <- match.arg(version)
  if (version == "2009") {
    paste0(
      "~I(pmin(141 * (0.993^age) * ",
      "(pmin(creat / ifelse(sex_num==1,0.7,0.9), 1)^ifelse(sex_num==1,-0.329,-0.411)) * ",
      "(pmax(creat / ifelse(sex_num==1,0.7,0.9), 1)^-1.209) * ",
      "ifelse(sex_num==1,1.018,1), 120))"
    )
  } else {
    paste0(
      "~I(pmin(142 * (0.9938^age) * ",
      "(pmin(creat / ifelse(sex_num==1,0.7,0.9), 1)^ifelse(sex_num==1,-0.241,-0.302)) * ",
      "(pmax(creat / ifelse(sex_num==1,0.7,0.9), 1)^-1.2) * ",
      "ifelse(sex_num==1,1.012,1), 120))"
    )
  }
}

dataset_ram_estimate <- function(df,
                                 m = 20,
                                 sample_n = 1e4,
                                 fudge = getOption("mi_ram_fudge", 4)) {
  n_total <- nrow(df)
  if (n_total < 1e4 && fudge == 4) fudge <- 2
  n_samp <- min(sample_n, n_total)
  is_posixt <- function(x) inherits(x, "POSIXt")
  df_samp <- df[sample(seq_len(n_total), n_samp), , drop = FALSE] |>
    mutate(
      across(where(is.factor), as.character),
      across(where(is_posixt), ~ as.numeric(.x))
    )
  bytes_per_row <- as.numeric(object.size(df_samp)) / n_samp
  bytes_total <- bytes_per_row * n_total * m * fudge
  bytes_total
}

with_progress_if <- function(expr, ...) {
  if (requireNamespace("progressr", quietly = TRUE)) {
    progressr::with_progress(expr, ...)
  } else {
    force(expr)
  }
}



mice_with_ridge <- function(data, ridge_seq = c(1e-3, 1e-2, 1e-1, 1),
                            ..., .engine = c("mice", "parlmice")) {
  .engine <- match.arg(.engine)
  FUN     <- if (.engine == "mice") mice::mice else mice::parlmice
  for (r in ridge_seq) {
    out <- tryCatch(
      FUN(data, ridge = r, 
          remove.collinear  = TRUE, remove.constant = TRUE, ...),
      error = function(e) e
    )
    if (!inherits(out, "error")) {
      message("✓ mice() succeeded with ridge = ", r)
      return(out)
    }
    message("… mice() failed with ridge = ", r, "  — trying next value")
  }
  stop("mice() failed for all ridge values")
}



impute_data <- function(df,
                        m        = getOption("mi_m", 20),
                        maxit    = getOption("mi_maxit", 10),
                        seed     = 1000,
                        printFlag = TRUE) {
  
  stopifnot(m >= 1)
  message(sprintf(" Multiple imputation started  (requested m = %d, maxit = %d)", m, maxit))
  
  df <- df %>%
    mutate(
      age     = leeftijd_jaar,
      across(c(study_grp, sex), as.factor),
      sex_num = ifelse(sex == "Vrouw", 1, 0)
    )
  
  nzv <- vapply(df, function(x) length(unique(na.omit(x))) > 1L, logical(1))
  if (any(!nzv)) {
    const_vars <- names(df)[!nzv]
    message("✓ Removing constant variable(s) from imputation: ",
            paste(const_vars, collapse = ", "))
    df <- df[, nzv, drop = FALSE]
  }
  
  pred <- mice::make.predictorMatrix(df)
  meth <- mice::make.method(df)
  post <- setNames(rep("", ncol(df)), names(df))
  derived_vars <- intersect(names(df),
                            c("ckd", "egfr_capped", "cr_resid", "renal_not_ordered"))
  if (length(derived_vars)) {
    meth[derived_vars]        <- ""
    pred[derived_vars, ]      <- 0
    pred[,   derived_vars  ]  <- 0
    message("✓ Excluding derived variable(s) from imputation: ",
            paste(derived_vars, collapse = ", "))
  }
  if ("creat" %in% names(df)) {
    message("✓ Imputing 'creat' by simple random sampling (fallback).")
    meth["creat"]  <- "sample"
    pred["creat",] <- 0
  }
  
  limit_predictors <- function(var, keep) {
    pred[var, ] <- 0
    pred[var, intersect(colnames(pred), keep)] <- 1
  }
  if ("egfr" %in% rownames(pred)) {
    limit_predictors("egfr", c("age", "sex", "bmi", "creat"))
  }
  
  drop_vars <- intersect(names(df),
                         c("creat_raw", "creat_conv", "egfr_raw",
                           "egfr_capped",
                           "length_of_stay_totaal_dagen"))
  
  if (length(drop_vars)) {
    message("✓ Dropping unused raw variable(s) from mice: ",
            paste(drop_vars, collapse = ", "))
    meth[drop_vars]       <- ""
    pred[drop_vars, ]     <- 0
    pred[,   drop_vars ]  <- 0
  }
  
  na_cnt   <- colSums(is.na(df))
  tiny_var <- names(na_cnt[na_cnt > 0 & na_cnt <= 2])
  for (v in tiny_var) {
    message(sprintf("✓ Imputing '%s' by random sampling (%d NA).", v, na_cnt[v]))
    meth[v]   <- "sample"
    pred[v, ] <- 0
  }
  
  if ("sex" %in% names(df) && anyNA(df$sex) && !("sex" %in% tiny_var)) {
    message("✓ Imputing 'sex' by random sampling (", sum(is.na(df$sex)), " NA).")
    meth["sex"]  <- "sample"
    pred["sex",] <- 0
  }
  
  if ("study_grp" %in% colnames(pred)) {
    pred[, "study_grp"] <- 0
    message("✓ Excluding 'study_grp' from predictor matrix to avoid collinearity.")
  }
  
  outcome_vars <- intersect(names(df), c("icu_death", "status", "mort1m", "mort3m"))
  if (length(outcome_vars)) {
    meth[outcome_vars]    <- ""
    pred[outcome_vars, ]  <- 0
    pred[, outcome_vars ] <- 0
    message("✓ Outcome variables excluded from imputation model: ",
            paste(outcome_vars, collapse = ", "))
  }
  
  meth["sex_num"]   <- ""
  pred["sex_num", ] <- 0
  pred[, "sex_num"] <- 0
  
  meth["bmi"] <- "pmm"   
  
  post["creat"] <- "imp[[j]][,i] <- pmin(pmax(imp[[j]][,i], 0.23), 17)"
  post["egfr"]  <- "imp[[j]][,i] <- pmin(pmax(imp[[j]][,i], 1),   120)"
  
  est_bytes <- dataset_ram_estimate(df, m = m)
  big_data <- est_bytes > 8e9
  if (big_data) {
    warning(sprintf(" Estimated RAM %.1f GiB > 8 GiB, reducing m to 5.",
                    est_bytes/2^30), call. = FALSE, immediate. = TRUE)
    m <- 5
  }
  
  can_future   <- requireNamespace("future", quietly = TRUE)
  has_multisess <- can_future && inherits(future::plan(), "multisession") &&
    future::nbrOfWorkers() > 1
  has_parlmice <- "parlmice" %in% utils::lsf.str("package:mice")
  allow_parallel <- getOption("mi_allow_parallel", TRUE)
  
  use_parlmice <- allow_parallel && !big_data && m <= 10 &&
    has_multisess  && has_parlmice
  
  engine <- if (use_parlmice) "parlmice" else "mice"
  
  imp <- with_progress_if({
    mice_with_ridge(
      df,
      m               = m,
      maxit           = maxit,
      predictorMatrix = pred,
      method          = meth,
      post            = post,
      seed            = seed,
      printFlag       = printFlag,
      .engine         = engine
    )
  })
  
  if (!inherits(imp, "mids")) {
    stop("Multiple imputation failed – `mice` did not return a mids object",
         call. = FALSE)
  }
  
  message("✓ Multiple imputation completed successfully.")
  imp
}

