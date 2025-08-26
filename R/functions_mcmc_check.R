suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
  library(splines)
  library(rjags)
  library(coda)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

.canon_inla_names <- function(nms) {
  if (is.null(nms)) return(nms)
  out <- nms
  lo  <- tolower(out)
  
  out[grepl("^\\s*log\\s*\\(\\s*bmi\\s*\\)\\s*$",  lo)] <- "log_bmi"
  out[grepl("^\\s*log\\s*\\(\\s*egfr\\s*\\)\\s*$", lo)] <- "log_egfr"
  
  pat_ns_df <- "^([a-z0-9._]+::)?\\s*ns\\s*\\(\\s*age(_c)?\\s*,\\s*(df\\s*=\\s*)?3\\s*\\)\\s*\\d+$"
  idx <- grepl(pat_ns_df, lo)
  if (any(idx)) {
    k <- sub(".*?(\\d+)\\s*$", "\\1", out[idx])
    out[idx] <- paste0("age_spline", k)
  }
  pat_ns_nodf <- "^([a-z0-9._]+::)?\\s*ns\\s*\\(\\s*age(_c)?\\s*\\)\\s*\\d+$"
  idx2 <- grepl(pat_ns_nodf, lo)
  if (any(idx2)) {
    k <- sub(".*?(\\d+)\\s*$", "\\1", out[idx2])
    out[idx2] <- paste0("age_spline", k)
  }
  
  out[grepl("^\\s*age[_\\s]*s\\s*1\\s*$", lo)] <- "age_spline1"
  out[grepl("^\\s*age[_\\s]*s\\s*2\\s*$", lo)] <- "age_spline2"
  out[grepl("^\\s*age[_\\s]*s\\s*3\\s*$", lo)] <- "age_spline3"
  
  out[grepl("mega.*pool", lo)] <- "mega_pool"
  
  make.unique(out, sep = "_")
}



.inla_fixed_df <- function(inla_fit_death) {
  stopifnot(is.list(inla_fit_death), !is.null(inla_fit_death$model))
  sf <- inla_fit_death$model$summary.fixed
  stopifnot(is.data.frame(sf), nrow(sf) > 0)
  
  rn     <- rownames(sf)
  std_rn <- .canon_inla_names(rn)
  
  b        <- setNames(sf$mean,         std_rn)
  sd_vec   <- setNames(sf$sd,           std_rn)
  q025_vec <- setNames(sf$`0.025quant`, std_rn)
  q975_vec <- setNames(sf$`0.975quant`, std_rn)
  s2z_2lvl <- function(b, l1, l2, tol = 1e-8) {
    all(c(l1, l2) %in% names(b)) && isTRUE(abs(unname(b[l1] + b[l2])) < tol)
  }
  out <- list(
    summary_fixed = sf,
    coef      = b,
    sd_vec    = sd_vec,
    q025_vec  = q025_vec,
    q975_vec  = q975_vec,
    age_mean      = inla_fit_death$age_mean      %||% attr(inla_fit_death, "age_mean"),
    age_knots_int = inla_fit_death$age_knots_int %||% attr(inla_fit_death, "age_knots_int"),
    age_knots_bd  = inla_fit_death$age_knots_bd  %||% attr(inla_fit_death, "age_knots_bd"),
    cuts          = inla_fit_death$cuts          %||% attr(inla_fit_death, "cuts")
  )
  
  if ("sexVrouw" %in% names(b) && "sexMan" %in% names(b)) {
    sc <- unname(b["sexVrouw"] - b["sexMan"])      
    if (s2z_2lvl(b, "sexVrouw", "sexMan")) sc <- sc / 2  
    out$sex_contrast <- sc
  } else if ("sexVrouw" %in% names(b)) {
    out$sex_contrast <- unname(b["sexVrouw"])
  } else if ("sex" %in% names(b)) {
    out$sex_contrast <- unname(b["sex"])
  } else {
    out$sex_contrast <- NA_real_
  }
  out
}


.inla_fixed_draws <- function(inla_result, n_draws = 4000) {
  if (is.null(inla_result) || !inherits(inla_result, "inla")) return(NULL)
  draws <- tryCatch(
    INLA::inla.posterior.sample(
      n = n_draws,
      result = inla_result,
      selection = list(fixed = TRUE, random = FALSE)
    ),
    error = function(e) NULL
  )
  if (is.null(draws)) return(NULL)
  sf <- inla_result$summary.fixed
  wanted <- .canon_inla_names(rownames(sf))    
  M <- matrix(NA_real_, nrow = length(draws), ncol = length(wanted))
  colnames(M) <- wanted
  for (i in seq_along(draws)) {
    vi <- draws[[i]]$fixed
    names(vi) <- .canon_inla_names(names(vi))  
    di <- as.numeric(vi[wanted])
    if (length(di) == length(wanted)) M[i, ] <- di
  }
  if (!any(is.finite(M))) return(NULL)
  M
}


make_mcmc_dataset <- function(cp_data,
                              cause = 1L,
                              n_segments = 3L,
                              inla_fit_death = NULL,
                              age_mean = NULL,
                              knots_int = NULL,
                              knots_bd = NULL,
                              cuts = NULL,
                              include_age_spline = TRUE,
                              include_mega_pool = TRUE) {
  
  stopifnot(is.data.frame(cp_data))
  
  if (!is.null(inla_fit_death)) {
    meta <- .inla_fixed_df(inla_fit_death)
    age_mean  <- meta$age_mean   %||% age_mean
    knots_int <- meta$age_knots_int %||% knots_int
    knots_bd  <- meta$age_knots_bd  %||% knots_bd
    cuts      <- meta$cuts          %||% cuts
  }
  
  dd <- cp_data %>%
    dplyr::filter(status %in% c(0L, cause)) %>%
    dplyr::mutate(event = as.integer(status == cause)) %>%
    dplyr::filter(stop > start)
  
  if (is.null(cuts)) {
    T_max <- max(dd$stop, na.rm = TRUE)
    cuts  <- seq(0, T_max, length.out = n_segments + 1L)
  } else {
    n_segments <- length(cuts) - 1L
  }
  
  dd_long <- survival::survSplit(
    Surv(start, stop, event) ~ .,
    data   = dd,
    cut    = cuts[-c(1, length(cuts))],
    episode = "interval",
    start  = "tstart",
    end    = "tstop",
    event  = "event"
  ) %>%
    dplyr::mutate(
      exposure = pmax(tstop - tstart, 0),
      interval = as.integer(factor(interval, levels = 1:n_segments)),
      log_egfr = log(egfr),
      log_bmi  = log(bmi),
      age_c    = (age - if (is.null(age_mean)) mean(age, na.rm = TRUE) else age_mean),
      sex_bin  = as.integer(sex == "Vrouw"),
      mega_bin = as.integer(!is.na(study_grp) & as.character(study_grp) == "Merged-small-sites")
    ) %>%
    dplyr::filter(exposure > 0)
  
  if (isTRUE(include_age_spline)) {
    df_age <- 3
    expected_int <- df_age - 1  
    expected_bd  <- 2           
    
    if (is.null(knots_int) || is.null(knots_bd) ||
        length(knots_int) != expected_int || length(knots_bd) != expected_bd) {
      tmpl <- splines::ns((dd$age - if (is.null(age_mean)) mean(dd$age, na.rm = TRUE) else age_mean),
                          df = df_age)
      knots_int <- attr(tmpl, "knots")
      knots_bd  <- attr(tmpl, "Boundary.knots")
    }
    
    age_spl <- splines::ns(dd_long$age_c, knots = knots_int, Boundary.knots = knots_bd)
    stopifnot(ncol(age_spl) == df_age)
    colnames(age_spl) <- paste0("age_s1", 1:df_age) 
    dd_long <- dplyr::bind_cols(dd_long, as.data.frame(age_spl))
  } else {
    dd_long <- dd_long %>% mutate(age_s11 = 0, age_s12 = 0, age_s13 = 0)
  }
  
  list(
    N          = nrow(dd_long),
    K          = n_segments,
    interval   = dd_long$interval,
    exposure   = dd_long$exposure,
    y          = dd_long$event,
    x_egfr     = dd_long$log_egfr,
    x_bmi      = dd_long$log_bmi,
    x_sex      = dd_long$sex_bin,
    x_mega     = if (isTRUE(include_mega_pool)) dd_long$mega_bin else rep(0L, nrow(dd_long)),
    x_a1       = dd_long$age_s11,
    x_a2       = dd_long$age_s12,
    x_a3       = dd_long$age_s13,
    .meta = list(
      include_age_spline = include_age_spline,
      include_mega_pool  = include_mega_pool,
      age_mean  = age_mean,
      knots_int = knots_int,
      knots_bd  = knots_bd
    )
  )
}



.jags_model_pwexp_string <- function(include_age_spline = TRUE,
                                     include_mega_pool  = TRUE,
                                     tau_fixed = 0.001,   
                                     tau_age   = 0.001) {
  pieces <- c(
    "model {",
    "  for (i in 1:N) {",
    "    log(mu[i]) <- log(exposure[i]) +",
    "                  beta0[interval[i]] +",
    "                  b_egfr * x_egfr[i] +",
    "                  b_bmi  * x_bmi[i]  +",
    "                  b_sex  * x_sex[i]"
  )
  if (include_mega_pool) {
    pieces <- c(pieces, "                  + b_mega * x_mega[i]")
  }
  if (include_age_spline) {
    pieces <- c(pieces,
                "                  + b_a1   * x_a1[i]",
                "                  + b_a2   * x_a2[i]",
                "                  + b_a3   * x_a3[i]")
  }
  pieces <- c(pieces,
              "    y[i] ~ dpois(mu[i])",
              "  }",
              "",
              "  for (k in 1:K) { beta0[k] ~ dnorm(0, 0.001) }",
              sprintf("  b_egfr ~ dnorm(0, %g)", tau_fixed),
              sprintf("  b_bmi  ~ dnorm(0, %g)", tau_fixed),
              sprintf("  b_sex  ~ dnorm(0, %g)", tau_fixed))
  if (include_mega_pool) {
    pieces <- c(pieces, sprintf("  b_mega ~ dnorm(0, %g)", tau_fixed))
  }
  if (include_age_spline) {
    pieces <- c(pieces,
                sprintf("  b_a1   ~ dnorm(0, %g)", tau_age),
                sprintf("  b_a2   ~ dnorm(0, %g)", tau_age),
                sprintf("  b_a3   ~ dnorm(0, %g)", tau_age))
  }
  pieces <- c(pieces, "}")
  paste(pieces, collapse = "\n")
}

fit_jags_pwexp <- function(data_list,
                           n_burn = 4000,
                           n_iter = 6000,
                           n_chains = 4,
                           seed = 123,
                           tau_fixed = 0.25,
                           tau_age   = 0.25) {
  set.seed(seed)
  inc_age  <- isTRUE(data_list$.meta$include_age_spline)
  inc_mega <- isTRUE(data_list$.meta$include_mega_pool)
  model_txt <- .jags_model_pwexp_string(inc_age, inc_mega, tau_fixed, tau_age)
  
  data_jags <- data_list
  if (!is.null(data_jags$.meta)) data_jags$.meta <- NULL
  
  jm <- rjags::jags.model(
    textConnection(model_txt),
    data     = data_jags,
    n.chains = n_chains,
    n.adapt  = 2000
  )
  update(jm, n_burn)
  
  vars <- c("beta0", "b_egfr", "b_bmi", "b_sex")
  if (inc_mega) vars <- c(vars, "b_mega")
  if (inc_age)  vars <- c(vars, "b_a1","b_a2","b_a3")
  
  rjags::coda.samples(jm, variable.names = vars, n.iter = n_iter, thin = 1)
}


.mcmc_convergence <- function(mcmc_samp) {
  out <- list(rhat_max = NA_real_, ess_min = NA_real_, ess_median = NA_real_)
  if (!inherits(mcmc_samp, "mcmc.list")) return(out)
  rhat <- tryCatch(gelman.diag(mcmc_samp, autoburnin = FALSE)$psrf[, 1], error = function(e) NULL)
  if (!is.null(rhat)) out$rhat_max <- max(rhat, na.rm = TRUE)
  ess  <- tryCatch(effectiveSize(mcmc_samp), error = function(e) NULL)
  if (!is.null(ess)) {
    ess <- as.numeric(ess)
    out$ess_min    <- min(ess, na.rm = TRUE)
    out$ess_median <- stats::median(ess, na.rm = TRUE)
  }
  out
}


compare_inla_mcmc <- function(inla_fit_death, mcmc_samp,
                              min_matches = 3,
                              use_inla_draws = TRUE,
                              n_inla_draws = 4000,
                              sign_flips = NULL) {
  meta <- .inla_fixed_df(inla_fit_death)
  
  b_inla    <- meta$coef
  sd_inla   <- meta$sd_vec
  q025_inla <- meta$q025_vec
  q975_inla <- meta$q975_vec
  
  names(b_inla)    <- .canon_inla_names(names(b_inla))
  names(sd_inla)   <- .canon_inla_names(names(sd_inla))
  names(q025_inla) <- .canon_inla_names(names(q025_inla))
  names(q975_inla) <- .canon_inla_names(names(q975_inla))
  
  M <- as.matrix(mcmc_samp)
  m_means <- colMeans(M, na.rm = TRUE)
  m_sd    <- apply(M, 2, stats::sd, na.rm = TRUE)
  m_q025  <- apply(M, 2, stats::quantile, probs = 0.025, na.rm = TRUE)
  m_q975  <- apply(M, 2, stats::quantile, probs = 0.975, na.rm = TRUE)
  
  map_mcmc <- c(
    log_egfr    = "b_egfr",
    log_bmi     = "b_bmi",
    sexVrouw    = "b_sex",
    mega_pool   = "b_mega",
    age_spline1 = "b_a1",
    age_spline2 = "b_a2",
    age_spline3 = "b_a3"
  )
  if (!("mega_pool" %in% names(b_inla))) {
    map_mcmc <- map_mcmc[setdiff(names(map_mcmc), "mega_pool")]
  }
  
  keep <- names(map_mcmc)[map_mcmc %in% colnames(M)]
  if (length(keep) < min_matches) {
    stop("Too few matched parameters (", length(keep), "); expected at least ", min_matches, ").")
  }
  
  tab <- tibble::tibble(
    par        = keep,
    mean_inla  = c(
      log_egfr    = b_inla["log_egfr"],
      log_bmi     = b_inla["log_bmi"],
      sexVrouw    = meta$sex_contrast,
      mega_pool   = if ("mega_pool" %in% names(b_inla)) unname(b_inla["mega_pool"]) else NA_real_,
      age_spline1 = b_inla["age_spline1"],
      age_spline2 = b_inla["age_spline2"],
      age_spline3 = b_inla["age_spline3"]
    )[keep],
    sd_inla    = c(
      log_egfr    = sd_inla["log_egfr"],
      log_bmi     = sd_inla["log_bmi"],
      sexVrouw    = NA_real_,
      mega_pool   = if ("mega_pool" %in% names(sd_inla)) sd_inla["mega_pool"] else NA_real_,
      age_spline1 = sd_inla["age_spline1"],
      age_spline2 = sd_inla["age_spline2"],
      age_spline3 = sd_inla["age_spline3"]
    )[keep],
    q025_inla  = c(
      log_egfr    = q025_inla["log_egfr"],
      log_bmi     = q025_inla["log_bmi"],
      sexVrouw    = NA_real_,
      mega_pool   = if ("mega_pool" %in% names(q025_inla)) q025_inla["mega_pool"] else NA_real_,
      age_spline1 = q025_inla["age_spline1"],
      age_spline2 = q025_inla["age_spline2"],
      age_spline3 = q025_inla["age_spline3"]
    )[keep],
    q975_inla  = c(
      log_egfr    = q975_inla["log_egfr"],
      log_bmi     = q975_inla["log_bmi"],
      sexVrouw    = NA_real_,
      mega_pool   = if ("mega_pool" %in% names(q975_inla)) q975_inla["mega_pool"] else NA_real_,
      age_spline1 = q975_inla["age_spline1"],
      age_spline2 = q975_inla["age_spline2"],
      age_spline3 = q975_inla["age_spline3"]
    )[keep],
    mean_mcmc  = as.numeric(m_means[map_mcmc[keep]]),
    sd_mcmc    = as.numeric(m_sd[map_mcmc[keep]]),
    q025_mcmc  = as.numeric(m_q025[map_mcmc[keep]]),
    q975_mcmc  = as.numeric(m_q975[map_mcmc[keep]])
  )
  
  needed <- intersect(c("log_egfr","log_bmi","age_spline1","age_spline2","age_spline3"), tab$par)
  miss <- tab$par %in% needed & !is.finite(tab$mean_inla)
  if (any(miss)) {
    warning("INLA means missing for: ",
            paste(tab$par[miss], collapse = ", "),
            " — please check .canon_inla_names() and rownames(summary.fixed).")
  }
  
  
  if (!is.null(sign_flips) && length(sign_flips)) {
    for (nm in intersect(names(sign_flips), tab$par)) {
      s <- sign_flips[[nm]]
      idx <- tab$par == nm
      tab$mean_mcmc[idx] <- s * tab$mean_mcmc[idx]
      if (s < 0) {
        q025 <- tab$q025_mcmc[idx]; q975 <- tab$q975_mcmc[idx]
        tab$q025_mcmc[idx] <- -q975; tab$q975_mcmc[idx] <- -q025
      }
    }
  }
  
  if (isTRUE(use_inla_draws)) {
    draws <- .inla_fixed_draws(inla_fit_death$model, n_draws = n_inla_draws)
    if (!is.null(draws)) {
      if (all(c("sexVrouw","sexMan") %in% colnames(draws))) {
        sex_draws <- draws[, "sexVrouw"] - draws[, "sexMan"]
      } else if ("sexVrouw" %in% colnames(draws)) {
        sex_draws <- draws[, "sexVrouw"]
      } else {
        sex_draws <- NULL
      }
      if (!is.null(sex_draws)) {
        tab$sd_inla[tab$par == "sexVrouw"]   <- stats::sd(sex_draws, na.rm = TRUE)
        tab$q025_inla[tab$par == "sexVrouw"] <- stats::quantile(sex_draws, 0.025, na.rm = TRUE)
        tab$q975_inla[tab$par == "sexVrouw"] <- stats::quantile(sex_draws, 0.975, na.rm = TRUE)
      }
    }
  }
  
  tab <- tab %>%
    dplyr::mutate(
      diff        = mean_mcmc - mean_inla,
      diff_sd     = diff / pmax(sd_inla, .Machine$double.eps),
      width_inla  = q975_inla - q025_inla,
      width_mcmc  = q975_mcmc - q025_mcmc,
      width_ratio = width_mcmc / pmax(width_inla, .Machine$double.eps)
    )
  
  attr(tab, "diag") <- list(matched = keep)
  tab
}


lp_agreement <- function(inla_fit_death, mcmc_samp, cp_data_one_imp,
                         include_mega_pool = TRUE,
                         align_sign = TRUE) {   
  meta  <- .inla_fixed_df(inla_fit_death)
  tab0  <- compare_inla_mcmc(inla_fit_death, mcmc_samp, use_inla_draws = FALSE)
  bI    <- setNames(tab0$mean_inla,  tab0$par)
  bMraw <- setNames(tab0$mean_mcmc,  tab0$par)
  bM    <- bMraw
  
  cp1 <- dplyr::distinct(cp_data_one_imp, patient, .keep_all = TRUE)
  age_mu <- meta$age_mean %||% mean(cp1$age, na.rm = TRUE)
  
  ki <- meta$age_knots_int; kb <- meta$age_knots_bd
  if (is.null(ki) || is.null(kb)) {
    tmp <- splines::ns(cp1$age - age_mu, df = 3)
    ki <- attr(tmp, "knots"); kb <- attr(tmp, "Boundary.knots")
  }
  
  X <- cp1 %>%
    mutate(
      log_egfr  = log(egfr),
      log_bmi   = log(bmi),
      sexVrouw  = as.numeric(sex == "Vrouw"),
      mega_pool = as.integer(!is.na(study_grp) & as.character(study_grp) == "Merged-small-sites"),
      age_c     = age - age_mu
    )
  
  B_age <- splines::ns(X$age_c, knots = ki, Boundary.knots = kb)
  
  mm_parts <- list(
    log_egfr = X$log_egfr,
    log_bmi  = X$log_bmi,
    sexVrouw = X$sexVrouw
  )
  if (isTRUE(include_mega_pool)) {
    mm_parts$mega_pool <- X$mega_pool
  }
  B_age_df <- as.data.frame(B_age)
  colnames(B_age_df) <- paste0("age_spline", 1:3)
  
  mm <- do.call(cbind, c(mm_parts, B_age_df))
  colnames(mm) <- make.names(colnames(mm), unique = TRUE)
  
  mm <- as.matrix(mm)
  storage.mode(mm) <- "double"
  
  storage.mode(bI) <- "double"
  storage.mode(bM) <- "double"
  
  have_both <- intersect(colnames(mm), intersect(names(bI), names(bM)))
  have_both <- have_both[is.finite(bI[have_both]) & is.finite(bM[have_both])]
  if (length(have_both) < 2) {
    stop("Not enough overlapping terms between INLA and MCMC for LP comparison. ",
         "Have: ", paste(have_both, collapse = ", "),
         ". Check name normalization and model inclusions.")
  }
  
  .lp_and_terms <- function(bM_use) {
    lp_inla <- as.numeric(mm[, have_both, drop = FALSE] %*% bI[have_both])
    lp_mcmc <- as.numeric(mm[, have_both, drop = FALSE] %*% bM_use[have_both])
    
    age_cols <- intersect(paste0("age_spline",1:3), have_both)
    lp_inla_noage <- if (length(age_cols)) lp_inla - as.numeric(mm[, age_cols, drop=FALSE] %*% bI[age_cols]) else lp_inla
    lp_mcmc_noage <- if (length(age_cols)) lp_mcmc - as.numeric(mm[, age_cols, drop=FALSE] %*% bM_use[age_cols]) else lp_mcmc
    
    terms <- c("egfr","bmi","sex","mega_pool","age")
    contrib_inla <- list(
      egfr = mm[,"log_egfr"] * (bI["log_egfr"] %||% NA_real_),
      bmi  = mm[,"log_bmi"]  * (bI["log_bmi"]  %||% NA_real_),
      sex  = mm[,"sexVrouw"] * (bI["sexVrouw"] %||% NA_real_),
      mega_pool = if ("mega_pool" %in% colnames(mm) && "mega_pool" %in% names(bI))
        mm[,"mega_pool"] * bI["mega_pool"] else rep(NA_real_, nrow(mm)),
      age = {
        age_cols <- intersect(paste0("age_spline",1:3), have_both)
        if (length(age_cols)) as.numeric(mm[, age_cols, drop=FALSE] %*% bI[age_cols]) else rep(NA_real_, nrow(mm))
      }
    )
    contrib_mcmc <- list(
      egfr = mm[,"log_egfr"] * (bM_use["log_egfr"] %||% NA_real_),
      bmi  = mm[,"log_bmi"]  * (bM_use["log_bmi"]  %||% NA_real_),
      sex  = mm[,"sexVrouw"] * (bM_use["sexVrouw"] %||% NA_real_),
      mega_pool = if ("mega_pool" %in% colnames(mm) && "mega_pool" %in% names(bM_use))
        mm[,"mega_pool"] * bM_use["mega_pool"] else rep(NA_real_, nrow(mm)),
      age = {
        age_cols <- intersect(paste0("age_spline",1:3), have_both)
        if (length(age_cols)) as.numeric(mm[, age_cols, drop=FALSE] %*% bM_use[age_cols]) else rep(NA_real_, nrow(mm))
      }
    )
    
    per_term <- lapply(terms, function(tn){
      xi <- contrib_inla[[tn]]; yi <- contrib_mcmc[[tn]]
      okt <- is.finite(xi) & is.finite(yi)
      if (!any(okt)) return(data.frame(term = tn, r = NA_real_, slope = NA_real_, intercept = NA_real_))
      fit <- lm(yi[okt] ~ xi[okt])
      data.frame(term = tn,
                 r = suppressWarnings(cor(xi[okt], yi[okt])),
                 slope = unname(coef(fit)[2]),
                 intercept = unname(coef(fit)[1]))
    })
    per_term <- do.call(rbind, per_term)
    
    ok  <- is.finite(lp_inla) & is.finite(lp_mcmc)
    okn <- is.finite(lp_inla_noage) & is.finite(lp_mcmc_noage)
    
    list(
      r_total  = suppressWarnings(cor(lp_inla[ok], lp_mcmc[ok])),
      r_no_age = suppressWarnings(cor(lp_inla_noage[okn], lp_mcmc_noage[okn])),
      per_term = per_term
    )
  }
  
  base <- .lp_and_terms(bM)
  
  sign_flips <- c()
  if (align_sign) {
    flip_terms <- c()
    pt <- base$per_term
    if ("sex" %in% pt$term && is.finite(pt$r[pt$term=="sex"]) && pt$r[pt$term=="sex"] < 0) {
      flip_terms <- c(flip_terms, "sexVrouw")
    }
    if ("mega_pool" %in% pt$term && is.finite(pt$r[pt$term=="mega_pool"]) && pt$r[pt$term=="mega_pool"] < 0) {
      flip_terms <- c(flip_terms, "mega_pool")
    }
    if (length(flip_terms)) {
      bM[flip_terms] <- -bM[flip_terms]
      sign_flips <- setNames(rep(-1, length(flip_terms)), flip_terms)
      base <- .lp_and_terms(bM)  
    }
  }
  
  list(
    r_total   = base$r_total,
    rmse      = sqrt(mean((as.numeric(mm[, have_both, drop = FALSE] %*% bI[have_both]) -
                             as.numeric(mm[, have_both, drop = FALSE] %*% bM[have_both]))^2)),
    r_no_age  = base$r_no_age,
    per_term  = base$per_term,
    have_both = have_both,
    sign_flips = sign_flips        
  )
}


inla_mcmc_agreement <- function(inla_fit, cp_data_one_imp,
                                n_segments = 3,
                                n_burn = 1000, n_iter = 1000, n_chains = 4,
                                seed = 123,
                                include_age_spline = TRUE,
                                include_mega_pool  = TRUE,
                                tau_fixed = 0.001, tau_age = 0.001,
                                use_inla_draws = TRUE, n_inla_draws = 4000) {
  
  stopifnot(is.list(inla_fit), !is.null(inla_fit$model))
  
  data_list <- make_mcmc_dataset(
    cp_data_one_imp,
    cause = 1L,
    n_segments = n_segments,
    inla_fit_death = inla_fit,
    include_age_spline = include_age_spline,
    include_mega_pool  = include_mega_pool
  )
  
  mcmc_samp <- fit_jags_pwexp(
    data_list,
    n_burn = n_burn, n_iter = n_iter, n_chains = n_chains, seed = seed,
    tau_fixed = tau_fixed, tau_age = tau_age
  )
  
  conv <- .mcmc_convergence(mcmc_samp)
  
  agree <- lp_agreement(inla_fit, mcmc_samp, cp_data_one_imp,
                        include_mega_pool = include_mega_pool,
                        align_sign = TRUE)
  
  tab  <- compare_inla_mcmc(inla_fit, mcmc_samp,
                            use_inla_draws = use_inla_draws,
                            n_inla_draws = n_inla_draws,
                            sign_flips = agree$sign_flips)
  
  structure(
    list(summary = tab, convergence = conv, agreement = agree),
    class = "inla_mcmc_agreement"
  )
}


fill_inla_cols <- function(tab, inla_fit_death) {
  meta <- .inla_fixed_df(inla_fit_death)
  b    <- meta$coef
  sdv  <- meta$sd_vec
  ql   <- meta$q025_vec
  qu   <- meta$q975_vec
  
  fetch <- function(x, src) if (x %in% names(src)) unname(src[[x]]) else NA_real_
  
  for (i in seq_len(nrow(tab))) {
    p <- tab$par[i]
    if (is.na(tab$mean_inla[i])) {
      tab$mean_inla[i] <- switch(
        p,
        sexVrouw    = meta$sex_contrast,
        mega_pool   = fetch("mega_pool", b),
        log_egfr    = fetch("log_egfr", b),
        log_bmi     = fetch("log_bmi",  b),
        age_spline1 = fetch("age_spline1", b),
        age_spline2 = fetch("age_spline2", b),
        age_spline3 = fetch("age_spline3", b),
        tab$mean_inla[i]
      )
    }
    if (is.na(tab$sd_inla[i])) tab$sd_inla[i] <- switch(
      p,
      sexVrouw    = NA_real_,  
      mega_pool   = fetch("mega_pool", sdv),
      log_egfr    = fetch("log_egfr", sdv),
      log_bmi     = fetch("log_bmi",  sdv),
      age_spline1 = fetch("age_spline1", sdv),
      age_spline2 = fetch("age_spline2", sdv),
      age_spline3 = fetch("age_spline3", sdv),
      tab$sd_inla[i]
    )
    if (is.na(tab$q025_inla[i])) tab$q025_inla[i] <- switch(
      p,
      sexVrouw    = NA_real_,
      mega_pool   = fetch("mega_pool", ql),
      log_egfr    = fetch("log_egfr", ql),
      log_bmi     = fetch("log_bmi",  ql),
      age_spline1 = fetch("age_spline1", ql),
      age_spline2 = fetch("age_spline2", ql),
      age_spline3 = fetch("age_spline3", ql),
      tab$q025_inla[i]
    )
    if (is.na(tab$q975_inla[i])) tab$q975_inla[i] <- switch(
      p,
      sexVrouw    = NA_real_,
      mega_pool   = fetch("mega_pool", qu),
      log_egfr    = fetch("log_egfr", qu),
      log_bmi     = fetch("log_bmi",  qu),
      age_spline1 = fetch("age_spline1", qu),
      age_spline2 = fetch("age_spline2", qu),
      age_spline3 = fetch("age_spline3", qu),
      tab$q975_inla[i]
    )
  }
  
  tab$diff        <- tab$mean_mcmc - tab$mean_inla
  tab$diff_sd     <- tab$diff / pmax(tab$sd_inla, .Machine$double.eps)
  tab$width_inla  <- tab$q975_inla - tab$q025_inla
  tab$width_ratio <- tab$width_mcmc / pmax(tab$width_inla, .Machine$double.eps)
  tab
}

align_age_spline_signs <- function(tab, inla_fit_death, cp_df, r_cut = -0.9) {
  meta <- .inla_fixed_df(inla_fit_death)
  cp1  <- dplyr::distinct(cp_df, patient, .keep_all = TRUE)
  B    <- splines::ns(cp1$age - meta$age_mean,
                      knots = meta$age_knots_int,
                      Boundary.knots = meta$age_knots_bd)
  colnames(B) <- paste0("age_spline", 1:3)
  
  bI <- setNames(tab$mean_inla, tab$par)
  bM <- setNames(tab$mean_mcmc, tab$par)
  
  flips <- character(0)
  for (col in colnames(B)) {
    if (!(col %in% names(bI) && col %in% names(bM))) next
    xi <- B[, col] * bI[col]
    yi <- B[, col] * bM[col]
    r  <- suppressWarnings(stats::cor(xi, yi, use = "complete.obs"))
    if (is.finite(r) && r <= r_cut) {
      i <- tab$par == col
      
      tab$mean_mcmc[i] <- -tab$mean_mcmc[i]
      
      q025 <- tab$q025_mcmc[i]; q975 <- tab$q975_mcmc[i]
      tab$q025_mcmc[i] <- -q975
      tab$q975_mcmc[i] <- -q025
      flips <- c(flips, col)
    }
  }
  
  tab <- tab %>%
    dplyr::mutate(
      diff        = mean_mcmc - mean_inla,
      diff_sd     = diff / pmax(sd_inla, .Machine$double.eps),
      width_inla  = q975_inla - q025_inla,
      width_mcmc  = q975_mcmc - q025_mcmc,
      width_ratio = width_mcmc / pmax(width_inla, .Machine$double.eps)
    )
  attr(tab, "age_flips") <- flips
  tab
}

recompute_agree_from_means <- function(tab, inla_fit_death, cp_df) {
  meta <- .inla_fixed_df(inla_fit_death)
  cp1  <- dplyr::distinct(cp_df, patient, .keep_all = TRUE)
  
  B <- splines::ns(cp1$age - meta$age_mean,
                   knots = meta$age_knots_int,
                   Boundary.knots = meta$age_knots_bd)
  colnames(B) <- paste0("age_spline", 1:3)
  mm <- cbind(
    log_egfr  = log(cp1$egfr),
    log_bmi   = log(cp1$bmi),
    sexVrouw  = as.numeric(cp1$sex == "Vrouw"),
    mega_pool = as.integer(!is.na(cp1$study_grp) & as.character(cp1$study_grp) == "Merged-small-sites"),
    B
  )
  
  bI <- setNames(tab$mean_inla, tab$par)
  bM <- setNames(tab$mean_mcmc, tab$par)
  have_both <- intersect(colnames(mm), intersect(names(bI), names(bM)))
  have_both <- have_both[is.finite(bI[have_both]) & is.finite(bM[have_both])]
  
  lpI <- as.numeric(mm[, have_both, drop=FALSE] %*% bI[have_both])
  lpM <- as.numeric(mm[, have_both, drop=FALSE] %*% bM[have_both])
  age_cols <- intersect(paste0("age_spline",1:3), have_both)
  
  lpI_noage <- if (length(age_cols)) lpI - as.numeric(mm[,age_cols,drop=FALSE] %*% bI[age_cols]) else lpI
  lpM_noage <- if (length(age_cols)) lpM - as.numeric(mm[,age_cols,drop=FALSE] %*% bM[age_cols]) else lpM
  
  term_contrib <- function(term, bvec) {
    if (term == "egfr")      mm[,"log_egfr"]  * bvec["log_egfr"]
    else if (term == "bmi")  mm[,"log_bmi"]   * bvec["log_bmi"]
    else if (term == "sex")  mm[,"sexVrouw"]  * bvec["sexVrouw"]
    else if (term == "mega_pool" && "mega_pool" %in% colnames(mm) && "mega_pool" %in% names(bvec))
      mm[,"mega_pool"] * bvec["mega_pool"]
    else if (term == "age" && length(age_cols))
      as.numeric(mm[,age_cols,drop=FALSE] %*% bvec[age_cols])
    else rep(NA_real_, nrow(mm))
  }
  
  per_term <- c("egfr","bmi","sex","mega_pool","age")
  pt <- purrr::map_dfr(per_term, function(tn){
    xi <- term_contrib(tn, bI); yi <- term_contrib(tn, bM)
    ok <- is.finite(xi) & is.finite(yi)
    if (!any(ok)) return(tibble::tibble(term = tn, r = NA_real_, slope = NA_real_, intercept = NA_real_))
    fit <- stats::lm(yi[ok] ~ xi[ok])
    tibble::tibble(term = tn,
                   r = suppressWarnings(stats::cor(xi[ok], yi[ok])) ,
                   slope = unname(stats::coef(fit)[2]),
                   intercept = unname(stats::coef(fit)[1]))
  })
  
  list(
    r_total  = suppressWarnings(stats::cor(lpI, lpM)),
    r_no_age = suppressWarnings(stats::cor(lpI_noage, lpM_noage)),
    per_term = pt
  )
}

postprocess_agree_branch <- function(obj, inla_fit_death, cp_df) {
  obj$summary <- fill_inla_cols(obj$summary, inla_fit_death)
  obj$summary <- align_age_spline_signs(obj$summary, inla_fit_death, cp_df)
  ag <- recompute_agree_from_means(obj$summary, inla_fit_death, cp_df)
  obj$agreement <- list(
    r_total  = ag$r_total,
    r_no_age = ag$r_no_age,
    per_term = ag$per_term
  )
  obj
}

