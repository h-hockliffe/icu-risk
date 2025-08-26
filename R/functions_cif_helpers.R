suppressPackageStartupMessages({
  library(survival)
  library(riskRegression)
})

.dedup_patients <- function(df) {
  if (!"patient" %in% names(df)) stop("newdata must contain 'patient'")
  df[!duplicated(df$patient), , drop = FALSE]
}
.get_age_mu <- function(fit, nd) {
  am <- attr(fit, "age_mean", exact = TRUE)
  if (is.finite(am)) am else mean(nd$age, na.rm = TRUE)
}
.prep_covariates <- function(nd, age_mu) {
  if (!"egfr" %in% names(nd) || !"bmi" %in% names(nd) || !"age" %in% names(nd))
    stop("newdata must include 'egfr', 'bmi', and 'age'")
  nd$log_egfr <- log(nd$egfr)
  nd$log_bmi  <- log(nd$bmi)
  nd$age_c    <- nd$age - age_mu
  if (!"mega_pool" %in% names(nd)) nd$mega_pool <- 0L
  if (!is.factor(nd$sex)) nd$sex <- factor(nd$sex)
  nd
}
.norm_pr <- function(pr, n_expected) {
  if (is.null(pr)) return(rep(NA_real_, n_expected))
  if (is.matrix(pr)) pr <- pr[, 1, drop = TRUE]
  if (length(pr) != n_expected)
    stop("Prediction length mismatch: got ", length(pr),
         " expected ", n_expected)
  as.numeric(pr)
}


add_age_splines <- function(df, knots_int, knots_bd) {
  if (is.null(knots_int) || is.null(knots_bd)) {
    knots_int <- stats::quantile(df$age_c, c(.25, .5, .75), na.rm = TRUE)
    knots_bd  <- range(df$age_c, na.rm = TRUE)
    warning("add_age_splines(): knots unavailable – using empirical quantiles.")
  }
  spl <- splines::ns(df$age_c, knots = knots_int, Boundary.knots = knots_bd) |>
    as.data.frame()
  names(spl) <- paste0("age_spline", 1:3)
  dplyr::bind_cols(df, spl)
}

lp_inla <- function(df, fe, re_vec) {
  if (is.data.frame(fe)) fe <- setNames(as.numeric(fe$mean), rownames(fe))
  if (is.matrix(fe) && !is.null(rownames(fe))) fe <- setNames(as.numeric(fe[, "mean"]), rownames(fe))
  stopifnot(!is.null(names(fe)))
  
  X <- if ("(Intercept)" %in% names(fe)) {
    m <- matrix(1, nrow(df), 1); colnames(m) <- "(Intercept)"; m
  } else NULL
  
  add_col <- function(nm, v) {
    if (is.null(X)) { X <<- cbind(v); colnames(X) <<- nm }
    else            { X <<- cbind(X, v); colnames(X)[ncol(X)] <<- nm }
  }
  
  for (nm in setdiff(names(fe), "(Intercept)")) {
    if (nm %in% names(df)) { add_col(nm, df[[nm]]); next }
    
    if (grepl("^sex", nm) && "sex" %in% names(df)) {
      lvl <- sub("^sex", "", nm)
      add_col(nm, as.integer(as.character(df$sex) == lvl))
      next
    }
    
    add_col(nm, 0)
  }
  
  fe_aligned <- fe[match(colnames(X), names(fe))]
  lp_fe <- as.numeric(X %*% fe_aligned)
  
  re_term <- re_vec[as.character(df$study_grp)]
  re_term[is.na(re_term)] <- 0
  lp_fe + re_term
}

cif_piecewise <- function(cuts, lambda1, lambda2,
                          scale1, scale2, horizon) {
  S <- 1
  cif <- 0
  for (j in seq_along(lambda1)) {
    t_start <- cuts[j]
    t_end <- cuts[j + 1]
    if (horizon <= t_start) break
    L <- min(horizon, t_end) - t_start
    if (L <= 0) next
    h1 <- lambda1[j] * scale1
    h2 <- lambda2[j] * scale2
    ht <- h1 + h2
    if (ht > 0) {
      cif <- cif + S * (h1 / ht) * (1 - exp(-ht * L))
      S <- S * exp(-ht * L)
    }
  }
  cif
}


predict_inla_cif <- function(inla_fit, newdata, horizon = 30) {
  d_fit <- inla_fit$death
  s_fit <- inla_fit$discharge
  if (is.null(d_fit) || is.null(s_fit)) {
    return(list(death = rep(NA_real_, nrow(newdata)),
                discharge = rep(NA_real_, nrow(newdata))))
  }
  if (is.null(d_fit$age_mean)) {
    stop("predict_inla_cif(): training age_mean not found – re-fit model.")
  }
  age_ctr <- d_fit$age_mean
  
  newdata <- newdata %>%
    dplyr::mutate(
      study_grp = as.character(study_grp),
      log_egfr  = log(egfr),
      log_bmi   = log(bmi),
      age_c     = age - age_ctr,
      mega_pool = as.integer(study_grp == "Merged-small-sites")
    ) %>%
    add_age_splines(d_fit$age_knots_int, d_fit$age_knots_bd)
  
  sf1 <- d_fit$model$summary.fixed
  sf2 <- s_fit$model$summary.fixed
  fe1 <- setNames(sf1$mean, rownames(sf1))
  fe2 <- setNames(sf2$mean, rownames(sf2))
  
  needs_cr <- any(grepl("^cr_resid$", names(fe1))) || any(grepl("^cr_resid$", names(fe2)))
  if (needs_cr && !("cr_resid" %in% names(newdata))) {
    warning("predict_inla_cif(): model includes 'cr_resid' but it is missing in newdata; treating as 0.",
            call. = FALSE, immediate. = TRUE)
  }
  
  re1 <- setNames(d_fit$model$summary.random$study_grp$mean,
                  d_fit$model$summary.random$study_grp$ID)
  re2 <- setNames(s_fit$model$summary.random$study_grp$mean,
                  s_fit$model$summary.random$study_grp$ID)
  
  if (any(!newdata$study_grp %in% names(re1))) {
    unknown <- unique(newdata$study_grp[!newdata$study_grp %in% names(re1)])
    warning("predict_inla_cif(): unseen study_grp: ", paste(unknown, collapse = ", "),
            " – assigning mean random effect.")
  }
  
  lp1 <- lp_inla(newdata, fe1, re1)
  lp2 <- lp_inla(newdata, fe2, re2)
  
  cuts <- d_fit$cuts
  lam1 <- exp(d_fit$model$summary.random$interval$mean)
  lam2 <- exp(s_fit$model$summary.random$interval$mean)
  
  cif_death <- vapply(
    seq_len(nrow(newdata)), function(i) {
      cif_piecewise(cuts, lam1, lam2, scale1 = exp(lp1[i]), scale2 = exp(lp2[i]), horizon = horizon)
    }, numeric(1)
  )
  
  cif_discharge <- vapply(
    seq_len(nrow(newdata)), function(i) {
      cif_piecewise(cuts, lam2, lam1, scale1 = exp(lp2[i]), scale2 = exp(lp1[i]), horizon = horizon)
    }, numeric(1)
  )
  
  list(death = as.numeric(cif_death), discharge = as.numeric(cif_discharge))
}


cox_cif <- function(cox_fit, newdata, horizon = 30, cause = 1) {
  if (is.null(cox_fit$death) || is.null(cox_fit$discharge))
    return(rep(NA_real_, nrow(newdata)))
  
  nd <- .dedup_patients(newdata)
  age_mu <- .get_age_mu(cox_fit$death, nd)
  nd <- .prep_covariates(nd, age_mu)
  
  csc <- tryCatch(
    riskRegression::CSC(list(cox_fit$death, cox_fit$discharge)),
    error = function(e) NULL
  )
  if (is.null(csc)) return(rep(NA_real_, nrow(nd)))
  nd_pred <- nd
  if ("status" %in% names(nd_pred)) nd_pred$status <- NULL  
  out <- suppressWarnings(
    tryCatch(
      riskRegression::predictRisk(csc, newdata = nd_pred, times = horizon,
                                  cause = 1, product.limit = FALSE)[, 1],
      error = function(e) rep(NA_real_, nrow(newdata))
    )
  )
  pr <- tryCatch(
    riskRegression::predictRisk(csc, newdata = nd, times = horizon, cause = cause),
    error = function(e) NULL
  )
  .norm_pr(pr, nrow(nd))
}

predict_cox_cif <- function(cox_list, newdata, horizon = 30,
                            include_center_fixed = FALSE) {
  if (is.null(cox_list$death) || is.null(cox_list$discharge)) {
    return(rep(NA_real_, nrow(newdata)))
  }
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  base_df <- attr(cox_list, "csc_base",        exact = TRUE)
  f1      <- attr(cox_list, "csc_form1_nore") %||% attr(cox_list, "csc_form1")
  f2      <- attr(cox_list, "csc_form2_nore") %||% attr(cox_list, "csc_form2")
  
  age_mu <- attr(cox_list$death, "age_mean", exact = TRUE)
  if (!is.finite(age_mu)) age_mu <- mean(newdata$age, na.rm = TRUE)
  
  nd <- newdata |>
    dplyr::mutate(
      log_egfr  = log(egfr),
      log_bmi   = log(bmi),
      age_c     = age - age_mu,
      sex       = if (is.factor(sex)) sex else factor(sex),
      study_grp = if (is.factor(study_grp)) study_grp else factor(study_grp)
    )
  
  if (is.data.frame(base_df)) {
    if ("sex" %in% names(base_df) && is.factor(base_df$sex)) {
      nd$sex <- factor(nd$sex, levels = levels(base_df$sex))
    }
    if ("study_grp" %in% names(base_df) && is.factor(base_df$study_grp)) {
      nd$study_grp <- factor(nd$study_grp, levels = levels(base_df$study_grp))
    } else {
      nd$study_grp <- factor(nd$study_grp)
    }
  } else {
    nd$sex       <- if (is.factor(nd$sex)) nd$sex else factor(nd$sex)
    nd$study_grp <- if (is.factor(nd$study_grp)) nd$study_grp else factor(nd$study_grp)
  }
  
  if (is.null(base_df) || is.null(f1) || is.null(f2)) {
    time_col   <- intersect(c("time", "stop", "tstop", "ftime"), names(newdata))[1]
    status_col <- intersect(c("status_cr", "status", "event"),   names(newdata))[1]
    if (is.na(time_col) || is.na(status_col)) {
      return(rep(NA_real_, nrow(newdata)))
    }
    
    base_df <- nd |>
      dplyr::mutate(
        time      = as.numeric(.data[[time_col]]),
        status_cr = factor(as.integer(.data[[status_col]]), levels = c(0L, 1L, 2L))
      ) |>
      dplyr::distinct(patient, .keep_all = TRUE)
    
    if (include_center_fixed) {
      f1 <- prodlim::Hist(time, status_cr) ~ log_egfr + sex + splines::ns(age_c, 3) + log_bmi + study_grp
      f2 <- prodlim::Hist(time, status_cr) ~ log_egfr + sex + splines::ns(age_c, 3) + log_bmi + study_grp
    } else {
      base_df <- base_df |>
        dplyr::mutate(mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites"))
      nd <- nd |>
        dplyr::mutate(mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites"))
      
      f1 <- prodlim::Hist(time, status_cr) ~ log_egfr + sex + splines::ns(age_c, 3) + log_bmi + mega_pool
      f2 <- prodlim::Hist(time, status_cr) ~ log_egfr + sex + splines::ns(age_c, 3) + log_bmi + mega_pool
    }
  } else {
    if (include_center_fixed) {
      f1 <- tryCatch(update(f1, ~ . - mega_pool + study_grp), error = function(e) f1)
      f2 <- tryCatch(update(f2, ~ . - mega_pool + study_grp), error = function(e) f2)
    } else {
      f_terms <- unique(c(all.vars(f1), all.vars(f2)))
      if ("mega_pool" %in% f_terms && !"mega_pool" %in% names(nd)) {
        nd <- nd |>
          dplyr::mutate(mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites"))
      }
    }
  }
  
  csc <- tryCatch(
    riskRegression::CSC(list("1" = f1, "2" = f2), data = base_df, ties = "breslow"),
    error = function(e) NULL
  )
  if (is.null(csc)) return(rep(NA_real_, nrow(newdata)))
  
  out <- suppressWarnings(
    tryCatch(
      riskRegression::predictRisk(csc, newdata = nd, times = horizon, cause = 1, product.limit = FALSE)[, 1],
      error = function(e) rep(NA_real_, nrow(newdata))
    )
  )
  out[!is.finite(out)] <- NA_real_
  out <- pmin(pmax(out, 0), 1)  
  as.numeric(out)
}

fg_cif_fgr <- function(fg_fit, newdata, horizon = 30, cause = 1) {
  nd <- .dedup_patients(newdata)
  age_mu <- .get_age_mu(fg_fit, nd)
  nd <- .prep_covariates(nd, age_mu)
  
  pr <- tryCatch(
    riskRegression::predictRisk(fg_fit, newdata = nd, times = horizon, cause = cause),
    error = function(e) NULL
  )
  .norm_pr(pr, nrow(nd))
}

fg_cif_crr <- function(fg_fit, newdata, horizon = 30) {
  if (!requireNamespace("cmprsk", quietly = TRUE))
    return(rep(NA_real_, nrow(newdata)))
  
  nd <- .dedup_patients(newdata)
  age_mu <- .get_age_mu(fg_fit, nd)
  nd <- .prep_covariates(nd, age_mu)
  
  df_spl <- attr(fg_fit, "df_spline", exact = TRUE); if (is.null(df_spl)) df_spl <- 3
  egfr_ki <- attr(fg_fit, "egfr_knots_int", exact = TRUE)
  egfr_kb <- attr(fg_fit, "egfr_knots_bd",  exact = TRUE)
  age_ki  <- attr(fg_fit, "age_knots_int",  exact = TRUE)
  age_kb  <- attr(fg_fit, "age_knots_bd",   exact = TRUE)
  
  spl_egfr <- splines::ns(nd$log_egfr, df = df_spl, knots = egfr_ki, Boundary.knots = egfr_kb)
  colnames(spl_egfr) <- paste0("log_egfr_s", seq_len(ncol(spl_egfr)))
  spl_age  <- splines::ns(nd$age_c,   df = df_spl, knots = age_ki,  Boundary.knots = age_kb)
  colnames(spl_age)  <- paste0("age_c_s",    seq_len(ncol(spl_age)))
  dd_spl <- cbind(nd, as.data.frame(spl_egfr), as.data.frame(spl_age))
  
  rhs <- c(colnames(spl_egfr), "log_bmi", "sex", colnames(spl_age))
  if ("mega_pool" %in% names(fg_fit$coef) && !"mega_pool" %in% rhs) rhs <- c(rhs, "mega_pool")
  
  mm <- stats::model.matrix(
    stats::as.formula(paste("~", paste(rhs, collapse = " + "))),
    data = dd_spl,
    contrasts.arg = list(sex = stats::contr.treatment)
  )
  X <- mm[, -1, drop = FALSE]
  X <- X[, intersect(colnames(X), names(fg_fit$coef)), drop = FALSE]
  
  failcode_val <- if (!is.null(fg_fit$failcode)) fg_fit$failcode else 1
  
  pr <- tryCatch(
    cmprsk::predict.crr(object = fg_fit, cov1 = X, failcode = failcode_val,
                        times = as.numeric(horizon), se.fit = FALSE),
    error = function(e) NULL
  )
  if (is.matrix(pr)) pr <- pr[1, , drop = TRUE]
  .norm_pr(pr, nrow(nd))
}

predict_fg_cif <- function(fg_fit, newdata, horizon = 30) {
  if (is.null(fg_fit)) return(rep(NA_real_, nrow(newdata)))
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  nd <- newdata %>% dplyr::distinct(patient, .keep_all = TRUE)
  age_mu <- attr(fg_fit, "age_mean", exact = TRUE)
  if (!is.finite(age_mu)) age_mu <- mean(nd$age, na.rm = TRUE)
  nd <- nd %>%
    dplyr::mutate(
      log_egfr = log(egfr),
      log_bmi  = log(bmi),
      age_c    = age - age_mu
    )
  
  if (inherits(fg_fit, "FGR")) {
    df_spl  <- attr(fg_fit, "df_spline",      exact = TRUE) %||% 3
    egfr_ki <- attr(fg_fit, "egfr_knots_int", exact = TRUE)
    egfr_kb <- attr(fg_fit, "egfr_knots_bd",  exact = TRUE)
    age_ki  <- attr(fg_fit, "age_knots_int",  exact = TRUE)
    age_kb  <- attr(fg_fit, "age_knots_bd",   exact = TRUE)
    
    spl_egfr <- splines::ns(nd$log_egfr, df = df_spl, knots = egfr_ki, Boundary.knots = egfr_kb)
    spl_age  <- splines::ns(nd$age_c,    df = df_spl, knots = age_ki,  Boundary.knots = age_kb)
    colnames(spl_egfr) <- paste0("log_egfr_s", seq_len(ncol(spl_egfr)))
    colnames(spl_age)  <- paste0("age_c_s",    seq_len(ncol(spl_age)))
    nd2 <- cbind(nd, as.data.frame(spl_egfr), as.data.frame(spl_age))
    
    pr <- try(riskRegression::predictRisk(fg_fit, newdata = nd2, times = horizon), silent = TRUE)
    if (!inherits(pr, "try-error")) {
      v <- as.numeric(pr[, 1])
      if (length(v) == nrow(nd)) return(pmin(pmax(v, 0), 1))
    }
  }
  
  if (!requireNamespace("cmprsk", quietly = TRUE)) {
    warning("predict_fg_cif(): cmprsk not available – returning NA.")
    return(rep(NA_real_, nrow(nd)))
  }
  
  df_spl  <- attr(fg_fit, "df_spline",      exact = TRUE) %||% 3
  egfr_ki <- attr(fg_fit, "egfr_knots_int", exact = TRUE)
  egfr_kb <- attr(fg_fit, "egfr_knots_bd",  exact = TRUE)
  age_ki  <- attr(fg_fit, "age_knots_int",  exact = TRUE)
  age_kb  <- attr(fg_fit, "age_knots_bd",   exact = TRUE)
  
  spl_egfr <- splines::ns(nd$log_egfr, df = df_spl, knots = egfr_ki, Boundary.knots = egfr_kb)
  spl_age  <- splines::ns(nd$age_c,    df = df_spl, knots = age_ki,  Boundary.knots = age_kb)
  colnames(spl_egfr) <- paste0("log_egfr_s", seq_len(ncol(spl_egfr)))
  colnames(spl_age)  <- paste0("age_c_s",    seq_len(ncol(spl_age)))
  dd <- cbind(nd, as.data.frame(spl_egfr), as.data.frame(spl_age))
  
  rhs <- c(colnames(spl_egfr), "log_bmi", "sex", colnames(spl_age))
  if ("mega_pool" %in% names(fg_fit$coef) && !"mega_pool" %in% rhs) rhs <- c(rhs, "mega_pool")
  if (!"mega_pool" %in% names(dd)) dd$mega_pool <- 0L
  if (!is.factor(dd$sex)) dd$sex <- factor(dd$sex)
  
  mm <- stats::model.matrix(
    stats::as.formula(paste("~", paste(rhs, collapse = " + "))),
    data = dd,
    contrasts.arg = list(sex = stats::contr.treatment)
  )
  X  <- mm[, -1, drop = FALSE]
  X  <- X[, names(fg_fit$coef), drop = FALSE]
  nP <- nrow(X)
  
  pr <- cmprsk::predict.crr(object = fg_fit, cov1 = X, times = as.numeric(horizon), se.fit = FALSE)
  
  safe_get <- function(x, nm) if (is.list(x) && !is.null(x[[nm]])) x[[nm]] else NULL
  
  if (is.matrix(pr)) {
    M <- pr
    has_overall <- ncol(M) == (nP + 1L)
    
  } else if (is.list(pr)) {
    M <- safe_get(pr, "pred") %||% safe_get(pr, "cif") %||% safe_get(pr, "P1") %||% safe_get(pr, "est")
    if (is.null(M)) M <- pr[[1]]
    if (!is.matrix(M)) M <- as.matrix(M)
    has_overall <- ncol(M) == (nP + 1L)
    
  } else {
    v <- as.numeric(pr)
    L <- length(v)
    
    if (L %% (nP + 1L) == 0L) {
      Tn <- L / (nP + 1L)
      M  <- matrix(v, nrow = Tn, ncol = nP + 1L, byrow = FALSE)
      has_overall <- TRUE
    } else if (L %% nP == 0L) {
      Tn <- L / nP
      M  <- matrix(v, nrow = Tn, ncol = nP, byrow = FALSE)
      has_overall <- FALSE
    } else {
      stop("predict_fg_cif(): cannot reshape predictions (length = ", L, ").")
    }
  }
  
  tgrid <- attr(pr, "times") %||% attr(pr, "time") %||% rownames(M) %||%
    fg_fit$uftime %||% fg_fit$ftime %||% NULL
  if (!is.null(tgrid)) tgrid <- as.numeric(tgrid)
  
  pick_idx <- function(ts, h, nrowM) {
    if (!is.null(ts) && length(ts) == nrowM) {
      idx <- max(which(ts <= h))
      if (!is.finite(idx)) idx <- which.min(abs(ts - h))
      if (!is.finite(idx) || length(idx) == 0) idx <- 1L
      return(as.integer(idx))
    }
    as.integer(nrowM)
  }
  
  ridx <- pick_idx(tgrid, horizon, nrow(M))
  vv   <- M[ridx, , drop = TRUE]
  if (has_overall && length(vv) == nP + 1L) vv <- vv[-1]
  
  if (length(vv) != nP) {
    stop("predict_fg_cif(): final length mismatch: got ", length(vv), " expected ", nP)
  }
  
  v <- as.numeric(vv)
  v[!is.finite(v)] <- NA_real_
  pmin(pmax(v, 0), 1)
}
