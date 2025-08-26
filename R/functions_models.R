suppressPackageStartupMessages({
  library(dplyr)
  library(splines)
  library(INLA)
  library(survival)
  library(riskRegression)
  library(frailtypack) 
  requireNamespace("progressr", quietly = TRUE)
  requireNamespace("later", quietly = TRUE)
  requireNamespace("cmprsk", quietly = TRUE)
  requireNamespace("memuse", quietly = TRUE) 
})

filter_valid_intervals <- function(d, verbose = TRUE, label = "intervals") {
  stopifnot(all(c("start","stop") %in% names(d)))
  ok <- is.finite(d$start) & is.finite(d$stop) & (d$stop > d$start)
  n_bad <- sum(!ok)
  if (n_bad > 0L && isTRUE(verbose)) {
    message(sprintf("%s: dropping %d row(s) with stop <= start or missing times.",
                    label, n_bad))
  }
  d[ok, , drop = FALSE]
}

.resolve_spec <- function(spec) {
  base <- list(age_df=3L, bmi_spline=NULL, include_mega_pool=TRUE,
               include_cr_resid=FALSE, egfr_spline_df=3L,
               egfr_as_quartiles=FALSE, include_centre_re=TRUE)
  if (is.null(spec)) return(base)
  if (!is.list(spec)) stop(".resolve_spec(): spec must be a list.", call. = FALSE)
  unknown <- setdiff(names(spec), names(base))
  if (length(unknown)) stop(".resolve_spec(): unknown spec field(s): ", paste(unknown, collapse=", "), call. = FALSE)
  modifyList(base, spec)
}


inla_verb <- getOption("inla_verbose", TRUE)

log_timed <- function(stage, start_time = NULL) {
  if (is.null(start_time)) {
    message(sprintf(
      "[%s] %s …",
      format(Sys.time(), "%Y‑%m‑%d %H:%M:%S"),
      stage
    ))
    invisible(Sys.time())
  } else {
    dt <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
    message(sprintf(
      "[%s] %s finished in %.1f min",
      format(Sys.time(), "%Y‑%m‑%d %H:%M:%S"),
      stage, dt
    ))
  }
}

schedule_heartbeat <- function(msg = "… still fitting …",
                               every = 120) {
  if (!requireNamespace("later", quietly = TRUE)) {
    return(function() NULL)
  } 
  
  done <- FALSE
  ping <- function() {
    if (!done) {
      message(msg)
      later::later(ping, every)
    }
  }
  later::later(ping, every)
  
  function() done <<- TRUE 
}

add_spline_if_better <- function(base_formula, data, var_term) {
  tl <- attr(terms(base_formula), "term.labels")
  if (is.null(tl)) tl <- character(0)
  tl <- setdiff(tl, var_term)  
  f0 <- if (length(tl)) {
    reformulate(tl, response = "Surv(start, stop, event)")
  } else {
    as.formula("Surv(start, stop, event) ~ 1")
  }
  
  f_spline <- tryCatch(
    update(f0, paste0("~ . + ns(", var_term, ", df = 3)")),
    error = function(e) NULL
  )
  
  m0 <- tryCatch(
    coxph(f0, data = data, ties = "breslow", x = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(m0)) return(base_formula)
  
  m1 <- if (!is.null(f_spline)) {
    tryCatch(coxph(f_spline, data = data, ties = "breslow", x = TRUE),
             error = function(e) NULL)
  } else NULL
  
  aic0 <- tryCatch(AIC(m0), error = function(e) NA_real_)
  aic1 <- tryCatch(if (!is.null(m1)) AIC(m1) else NA_real_,
                   error = function(e) NA_real_)
  
  if (is.finite(aic0) && is.finite(aic1) && (aic0 - aic1) >= 3) {
    if (!is.null(m1)) m1$formula else f_spline
  } else {
    f0
  }
}

predict_csc_safe <- function(object, newdata, times, cause = 1L) {
  m <- riskRegression::predictRisk(object, newdata = newdata, times = times,
                                   cause = cause, product.limit = FALSE)
  m[!is.finite(m)] <- NA_real_
  pmin(pmax(m, 0), 1)
}

predict_any <- function(fit, newdata, times, cause = 1L) {
  if (inherits(fit, "CauseSpecificCox")) {
    predict_csc_safe(fit, newdata, times, cause)
  } else {
    riskRegression::predictRisk(fit, newdata = newdata, times = times)
  }
}

`%||%` <- function(a,b) if (!is.null(a)) a else b

make_base_for_csc <- function(d, fit_like) {
  age_mu  <- attr(fit_like, "age_mean") %||% mean(d$age, na.rm = TRUE)
  use_q   <- isTRUE(attr(fit_like, "egfr_as_quartiles"))
  qs_all  <- attr(fit_like, "egfr_quartiles") %||% NULL
  
  base <- d %>%
    dplyr::filter(is.finite(start), is.finite(stop), stop > start, !is.na(status)) %>%
    dplyr::distinct(patient, .keep_all = TRUE) %>%
    dplyr::mutate(
      log_egfr  = log(egfr),
      log_bmi   = log(bmi),
      age_c     = age - age_mu,
      sex       = factor(sex),
      study_grp = droplevels(factor(study_grp)),
      patient   = factor(patient),
      time      = as.numeric(stop),
      status_cr = factor(as.integer(status), levels = c(0L,1L,2L),
                         labels = c("0","1","2")),
      mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites")
    )
  
  if (use_q) {
    qs <- if (!is.null(qs_all)) qs_all else
      stats::quantile(base$egfr, probs = c(.25,.5,.75), na.rm = TRUE, names = FALSE)
    base$egfr_q <- cut(base$egfr, breaks = c(-Inf, qs, Inf),
                       labels = c("Q1","Q2","Q3","Q4"),
                       include.lowest = TRUE, right = TRUE)
    base$egfr_q <- stats::relevel(base$egfr_q, ref = "Q4")
  }
  base
}

norm_cox_fit_to_csc <- function(fit_like, data_like) {
  if (inherits(fit_like, "CauseSpecificCox")) return(fit_like)
  if (is.list(fit_like) && all(c("death","discharge") %in% names(fit_like))) {
    f1 <- attr(fit_like, "csc_form1_nore") %||% attr(fit_like, "csc_form1_center")
    f2 <- attr(fit_like, "csc_form2_nore") %||% attr(fit_like, "csc_form2_center")
    stopifnot(inherits(f1, "formula"), inherits(f2, "formula"))
    base <- make_base_for_csc(data_like, fit_like)
    return(riskRegression::CSC(list("1" = f1, "2" = f2), data = base))
  }
  stop("Unexpected cox_fit class: ", paste(class(fit_like), collapse = ", "))
}

predict_csc_safe <- function(object, newdata, times, cause = 1L) {
  m <- riskRegression::predictRisk(
    object, newdata = newdata, times = times,
    cause = cause, product.limit = FALSE
  )
  m[!is.finite(m)] <- NA_real_
  pmin(pmax(m, 0), 1)
}


fit_cox <- function(cp_data,
                    patient_cluster = FALSE,
                    spec = NULL,
                    verbose = TRUE) {
  stopifnot(is.data.frame(cp_data))
  req <- c("start","stop","status","age","bmi","egfr","sex","study_grp","patient")
  miss <- setdiff(req, names(cp_data))
  if (length(miss)) stop("fit_cox(): cp_data is missing column(s): ", paste(miss, collapse=", "),
                         call. = FALSE)
  
  if (is.null(spec)) {
    if (exists("make_spec", mode = "function")) spec <- make_spec()
    else spec <- list(age_df = 3, bmi_spline = NULL, include_mega_pool = TRUE,
                      include_cr_resid = FALSE, egfr_spline_df = 3,
                      egfr_as_quartiles = FALSE, include_centre_re = TRUE,
                      age_tvc = NULL)  
  }
  or_default <- function(x, d) if (is.null(x)) d else x
  
  age_df        <- or_default(spec$age_df, 3L)
  egfr_df       <- or_default(spec$egfr_spline_df, 3L)
  use_quart     <- isTRUE(spec$egfr_as_quartiles)
  use_bmi_spline_requested <- spec$bmi_spline  
  include_mega_pool <- isTRUE(spec$include_mega_pool)
  include_cr_resid  <- isTRUE(spec$include_cr_resid)
  use_frailty_req   <- isTRUE(spec$include_centre_re)
  age_tvc <- spec$age_tvc  
  
  tt_age <- NULL
  if (!is.null(age_tvc) && !identical(age_tvc, FALSE)) {
    if (isTRUE(age_tvc) || identical(age_tvc, "log")) {
      tt_age <- function(x, t, ...) splines::ns(x, df = age_df) * log(pmax(t, 1e-6))
    } else if (is.function(age_tvc)) {
      tt_age <- function(x, t, ...) splines::ns(x, df = age_df) * age_tvc(t)
    } else if (is.character(age_tvc) && age_tvc == "identity") {
      tt_age <- function(x, t, ...) splines::ns(x, df = age_df) * t
    } else {
      tt_age <- function(x, t, ...) splines::ns(x, df = age_df) * log(pmax(t, 1e-6))
    }
  }
  use_tt_age <- !is.null(tt_age)
  
  use_cr <- include_cr_resid && "cr_resid" %in% names(cp_data) &&
    any(is.finite(cp_data$cr_resid))
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package is required.", call. = FALSE)
  }
  if (!requireNamespace("prodlim", quietly = TRUE) && isTRUE(verbose)) {
    message("ℹ prodlim not available: CSC formulas use prodlim::Hist and may fail downstream.")
  }
  
  age_mu <- mean(cp_data$age, na.rm = TRUE)
  
  qs_all <- NULL
  if (use_quart) {
    cp_ok <- filter_valid_intervals(cp_data, verbose = FALSE)
    egfr_vec <- cp_ok$egfr[is.finite(cp_ok$egfr)]
    if (length(egfr_vec) >= 5L) {
      qs_all <- stats::quantile(egfr_vec, probs = c(.25,.5,.75), na.rm = TRUE, names = FALSE)
    }
  }
  
  if (isTRUE(verbose)) {
    n_bad_egfr <- sum(is.finite(cp_data$egfr) & cp_data$egfr <= 0, na.rm = TRUE)
    n_bad_bmi  <- sum(is.finite(cp_data$bmi)  & cp_data$bmi  <= 0, na.rm = TRUE)
    if (n_bad_egfr > 0) message(" ", n_bad_egfr, " row(s) with non-positive eGFR (log undefined).")
    if (n_bad_bmi  > 0) message(" ", n_bad_bmi,  " row(s) with non-positive BMI (log undefined).")
  }
  
  safe_coxph <- function(form, data, tt_fun = NULL) {
    has_tt <- !is.null(tt_fun)
    tryCatch(
      survival::coxph(
        formula = form, data = data, ties = "breslow",
        x = TRUE, y = TRUE, model = if (has_tt) FALSE else TRUE, tt = tt_fun
      ),
      error = function(e) { if (isTRUE(verbose)) message("   coxph failed: ", conditionMessage(e)); NULL }
    )
  }
  
  
  
  build_rhs <- function(use_bmi_spline, use_bmi, use_mega, use_frailty, use_cluster, use_quart, use_cr, use_tt_age) {
    rhs <- c(
      if (use_quart) "egfr_q" else paste0("splines::ns(log_egfr, df = ", egfr_df, ")"),
      "sex",
      paste0("splines::ns(age_c, df = ", age_df, ")")  
    )
    if (isTRUE(use_tt_age)) rhs <- c(rhs, "tt(age_c)")  
    if (use_bmi) {
      if (isTRUE(use_bmi_spline)) rhs <- c(rhs, "splines::ns(log_bmi, df = 3)")
      else                        rhs <- c(rhs, "log_bmi")
    }
    if (isTRUE(use_mega))    rhs <- c(rhs, "mega_pool")
    if (isTRUE(use_cr))      rhs <- c(rhs, "cr_resid")
    if (isTRUE(use_frailty)) rhs <- c(rhs, "frailty(study_grp, distribution = 'gaussian')")
    if (isTRUE(use_cluster)) rhs <- c(rhs, "cluster(patient)")
    stats::as.formula(paste("Surv(start, stop, event) ~", paste(rhs, collapse = " + ")))
  }
  
  fit_one_cause <- function(k) {
    dd <- cp_data %>%
      dplyr::mutate(
        log_egfr  = log(egfr),
        log_bmi   = log(bmi),
        age_c     = age - age_mu,
        mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites"),
        patient   = factor(patient),
        sex       = factor(sex),
        study_grp = droplevels(factor(study_grp))
      )
    dd <- filter_valid_intervals(dd, verbose = verbose, label = sprintf("Cox cause %s", k))
    
    dd$event <- ifelse(is.na(dd$status), NA_integer_, as.integer(dd$status == k))
    
    if (use_quart) {
      qs <- if (!is.null(qs_all)) qs_all else
        stats::quantile(dd$egfr, probs = c(.25,.5,.75), na.rm = TRUE, names = FALSE)
      dd$egfr_q <- cut(dd$egfr, breaks = c(-Inf, qs, Inf),
                       labels = c("Q1","Q2","Q3","Q4"),
                       include.lowest = TRUE, right = TRUE)
      dd$egfr_q <- stats::relevel(dd$egfr_q, ref = "Q4")
    }
    
    keep <- is.finite(dd$log_egfr) & is.finite(dd$age_c) & !is.na(dd$sex)
    if (isTRUE(use_bmi_spline_requested) || identical(use_bmi_spline_requested, FALSE)) {
      keep <- keep & is.finite(dd$log_bmi)
    }
    if (include_mega_pool) keep <- keep & is.finite(dd$mega_pool)
    if (isTRUE(use_cr))     keep <- keep & is.finite(dd$cr_resid)
    
    n_drop <- sum(!keep)
    if (n_drop > 0L && isTRUE(verbose))
      message("Cox cause ", k, ": dropping ", n_drop, " row(s) with missing predictors before fitting.")
    dd <- dd[keep, , drop = FALSE]
    
    n_evt <- sum(dd$event, na.rm = TRUE)
    if (n_evt < 5L) {
      if (isTRUE(verbose)) message("Cox: cause ", k, " has < 5 events (n = ", n_evt, ") – skipped.")
      return(NULL)
    }
    
    n_centres <- nlevels(dd$study_grp)
    use_frailty_possible <- isTRUE(use_frailty_req) && n_centres >= 2L
    use_bmi_possible     <- all(is.finite(dd$log_bmi)) && stats::sd(dd$log_bmi) > 0
    
    prefer_spline <- NULL
    if (use_bmi_possible && is.null(use_bmi_spline_requested)) {
      dd_aic <- dd[is.finite(dd$log_bmi), , drop = FALSE]
      if (nrow(dd_aic) > 0L) {
        form_lin <- build_rhs(FALSE, TRUE, include_mega_pool, FALSE, FALSE, use_quart, use_cr, use_tt_age)
        form_spl <- build_rhs(TRUE,  TRUE, include_mega_pool, FALSE, FALSE, use_quart, use_cr, use_tt_age)
        m_lin <- safe_coxph(form_lin, dd_aic, tt_fun = tt_age)
        m_spl <- safe_coxph(form_spl, dd_aic, tt_fun = tt_age)
        a_lin <- tryCatch(AIC(m_lin), error = function(e) NA_real_)
        a_spl <- tryCatch(AIC(m_spl), error = function(e) NA_real_)
        if (isTRUE(verbose))
          message("   AIC check (", if (use_quart) "egfr=quartiles" else "egfr=spline",
                  "): linear=", round(a_lin,2), " | spline=", round(a_spl,2))
        if (is.finite(a_lin) && is.finite(a_spl)) {
          prefer_spline <- ((a_lin - a_spl) >= 3)
        }
      }
    }

    
    bmi_forms <- if (!use_bmi_possible) {
      "none"
    } else if (isTRUE(use_bmi_spline_requested)) {
      c("spline","linear")
    } else if (identical(use_bmi_spline_requested, FALSE)) {
      c("linear")
    } else if (is.logical(prefer_spline)) {
      if (prefer_spline) c("spline","linear") else c("linear","spline")
    } else {
      c("spline","linear")
    }
    
    try_grid <- list()
    for (fra in c(use_frailty_possible, FALSE)) {
      for (bmi in bmi_forms) {
        for (mega in c(include_mega_pool, FALSE)) {
          try_grid[[length(try_grid) + 1L]] <- list(
            use_frailty = fra,
            bmi_form    = bmi,
            use_bmi     = (bmi != "none"),
            use_mega    = mega
          )
        }
      }
    }
    try_grid[[length(try_grid) + 1L]] <- list(
      use_frailty = FALSE, bmi_form = "none", use_bmi = FALSE, use_mega = FALSE
    )
    
    mdl <- NULL
    tried <- 0L
    for (cfg in try_grid) {
      tried <- tried + 1L
      form <- build_rhs(
        use_bmi_spline = identical(cfg$bmi_form, "spline"),
        use_bmi        = isTRUE(cfg$use_bmi),
        use_mega       = isTRUE(cfg$use_mega),
        use_frailty    = isTRUE(cfg$use_frailty),
        use_cluster    = patient_cluster,
        use_quart      = use_quart,
        use_cr         = use_cr,
        use_tt_age     = use_tt_age
      )
      if (isTRUE(verbose))
        message(sprintf("   • Try #%d: %s | BMI=%s | mega=%s | cluster=%s",
                        tried,
                        if (cfg$use_frailty) "frailty" else "no frailty",
                        cfg$bmi_form, cfg$use_mega, patient_cluster))
      mdl <- safe_coxph(form, dd, tt_fun = tt_age)
      if (!is.null(mdl)) break
    }
    
    if (is.null(mdl)) {
      if (isTRUE(verbose)) message("Cox fit failed for cause ", k, " after all fallbacks.")
      return(NULL)
    }
    
    tf   <- stats::terms(form, keep.order = TRUE)
    labs <- attr(tf, "term.labels") %||% character(0)
    
    is_special <- grepl("\\b(frailty|cluster|tt|strata)\\s*\\(", labs)
    labs <- labs[!is_special]
    
    if (!length(labs)) {
      labs <- c(
        sprintf("splines::ns(log_egfr, df = %d)", egfr_df),
        "sex",
        sprintf("splines::ns(age_c, df = %d)", age_df)
      )
    }
    
    rhs_txt <- paste(labs, collapse = " + ")
    ph_form <- stats::as.formula(paste("survival::Surv(start, stop, event) ~", rhs_txt))
    attr(mdl, "ph_form") <- ph_form
    attr(mdl, "ph_df")   <- dd
    
    if (isTRUE(verbose)) message("   PH RHS: ", rhs_txt)
    
    
    y_count <- cbind(start = dd$start, stop = dd$stop, event = dd$event)
    nx <- NROW(mdl$x); ny <- NROW(y_count); n <- min(nx, ny)
    if (is.finite(n) && n > 0L) {
      if (ny != n) y_count <- y_count[seq_len(n), , drop = FALSE]
      attr(mdl, "ph_y_counting") <- y_count
      attr(mdl, "ph_x_colnames") <- colnames(mdl$x)
    } else if (isTRUE(verbose)) {
      message("   PH store: no overlapping rows for x/y; skipping PH attributes.")
    }
    
    attr(mdl, "age_mean") <- age_mu
    if (use_quart) attr(mdl, "egfr_quartiles") <- if (!is.null(qs_all)) qs_all else qs
    attr(mdl, "egfr_df")  <- egfr_df
    attr(mdl, "age_df")   <- age_df
    attr(mdl, "age_tvc")      <- age_tvc
    attr(mdl, "age_tvc_fun")  <- if (isTRUE(age_tvc) || identical(age_tvc, "log")) "ns(age_c, df=age_df) * log(t)" else
      if (is.character(age_tvc) && age_tvc == "identity") "ns(age_c, df=age_df) * t" else
        if (is.function(age_tvc)) "ns(age_c, df=age_df) * age_tvc(t)" else NULL
    return(mdl)
  }
  
  death     <- fit_one_cause(1L)
  discharge <- fit_one_cause(2L)
  
  dd_base <- cp_data %>%
    filter_valid_intervals(verbose = verbose, label = "CSC base") %>%   
    dplyr::mutate(
      log_egfr  = log(egfr),
      log_bmi   = log(bmi),
      age_c     = age - age_mu,
      sex       = factor(sex),
      study_grp = droplevels(factor(study_grp)),
      patient   = factor(patient),
      time      = as.numeric(stop),
      status_cr = factor(as.integer(status), levels = c(0L,1L,2L), labels = c("0","1","2")),
      mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites")
    ) %>%
    dplyr::filter(is.finite(time), time > 0, !is.na(status_cr)) %>%
    dplyr::distinct(patient, .keep_all = TRUE)
  
  if (use_quart) {
    qs <- if (!is.null(qs_all)) qs_all else
      stats::quantile(cp_data$egfr, probs = c(.25,.5,.75), na.rm = TRUE, names = FALSE)
    dd_base$egfr_q <- cut(dd_base$egfr, breaks = c(-Inf, qs, Inf),
                          labels = c("Q1","Q2","Q3","Q4"),
                          include.lowest = TRUE, right = TRUE)
    dd_base$egfr_q <- stats::relevel(dd_base$egfr_q, ref = "Q4")
  }
  
  has_term <- function(m, patt) {
    if (is.null(m)) return(FALSE)
    labs <- try(attr(stats::terms(m), "term.labels"), silent = TRUE)
    if (inherits(labs, "try-error") || is.null(labs)) return(FALSE)
    any(grepl(patt, labs))
  }
  
  keep_bmi_spline_death <- has_term(death,     "ns\\(log_bmi\\s*,\\s*df\\s*=\\s*3\\)")
  keep_bmi_spline_disc  <- has_term(discharge, "ns\\(log_bmi\\s*,\\s*df\\s*=\\s*3\\)")
  keep_bmi_linear_death <- !keep_bmi_spline_death && has_term(death,     "\\blog_bmi\\b")
  keep_bmi_linear_disc  <- !keep_bmi_spline_disc  && has_term(discharge, "\\blog_bmi\\b")
  keep_mega_death       <- has_term(death,     "\\bmega_pool\\b")
  keep_mega_disc        <- has_term(discharge, "\\bmega_pool\\b")
  keep_cr_resid_death   <- isTRUE(use_cr) && has_term(death, "\\bcr_resid\\b")
  keep_cr_resid_disc    <- isTRUE(use_cr) && has_term(discharge, "\\bcr_resid\\b")
  
  base_terms <- function(keep_bmi_spline, keep_bmi_linear, keep_mega, keep_cr, use_quart) {
    rhs <- c(
      if (use_quart) "egfr_q" else paste0("splines::ns(log_egfr, df = ", egfr_df, ")"),
      "sex",
      paste0("splines::ns(age_c, df = ", age_df, ")")
    )
    if (keep_bmi_spline) rhs <- c(rhs, "splines::ns(log_bmi, df = 3)")
    if (keep_bmi_linear) rhs <- c(rhs, "log_bmi")
    if (keep_mega)       rhs <- c(rhs, "mega_pool")
    if (keep_cr)         rhs <- c(rhs, "cr_resid")
    stats::as.formula(paste("prodlim::Hist(time, status_cr) ~", paste(rhs, collapse = " + ")))
  }
  
  form_death_csc_nore <- base_terms(keep_bmi_spline_death, keep_bmi_linear_death,
                                    keep_mega_death, keep_cr_resid_death, use_quart)
  form_disc_csc_nore  <- base_terms(keep_bmi_spline_disc,  keep_bmi_linear_disc,
                                    keep_mega_disc,  keep_cr_resid_disc,  use_quart)
  
  swap_mega_with_centre <- function(form) stats::update(form, ~ . - mega_pool + study_grp)
  form_death_csc_center <- if (keep_mega_death) swap_mega_with_centre(form_death_csc_nore)
  else stats::update(form_death_csc_nore, ~ . + study_grp)
  form_disc_csc_center  <- if (keep_mega_disc)  swap_mega_with_centre(form_disc_csc_nore)
  else stats::update(form_disc_csc_nore,  ~ . + study_grp)
  
  out <- list(death = death, discharge = discharge)
  
  attr(out, "csc_base")             <- dd_base
  attr(out, "csc_form1_nore")       <- form_death_csc_nore
  attr(out, "csc_form2_nore")       <- form_disc_csc_nore
  attr(out, "csc_form1_center")     <- form_death_csc_center
  attr(out, "csc_form2_center")     <- form_disc_csc_center
  attr(out, "levels_sex")           <- levels(dd_base$sex)
  attr(out, "levels_study_grp")     <- levels(dd_base$study_grp)
  attr(out, "age_mean")             <- age_mu
  attr(out, "egfr_df")              <- egfr_df
  attr(out, "age_df")               <- age_df
  attr(out, "egfr_as_quartiles")    <- use_quart
  attr(out, "bmi_spline_requested") <- use_bmi_spline_requested
  attr(out, "egfr_quartiles")       <- qs_all
  attr(out, "age_tvc")              <- age_tvc
  return(out)
}

fit_finegray <- function(cp_data,
                         event_code = 1L,
                         df_spline  = 3,
                         cv_fold    = 5,
                         engine     = c("auto", "FGR", "crr"),
                         patient_cluster = FALSE,
                         verbose    = TRUE,
                         spec = NULL) {
  
  stopifnot(is.data.frame(cp_data))
  req_cols <- c("stop", "status", "age", "bmi", "egfr", "sex", "study_grp", "patient")
  miss <- setdiff(req_cols, names(cp_data))
  if (length(miss)) stop("fit_finegray(): cp_data is missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  
  spec <- .resolve_spec(spec)
  engine <- match.arg(engine)
  
  age_mu <- mean(cp_data$age, na.rm = TRUE)
  dd <- cp_data %>%
    distinct(patient, .keep_all = TRUE) %>%
    mutate(
      log_egfr  = log(egfr),
      log_bmi   = log(bmi),
      age_c     = age - age_mu,
      mega_pool = as.integer(as.character(study_grp) == "Merged-small-sites"),
      sex       = factor(sex)
    )
  keep_fg <- is.finite(dd$log_egfr) & is.finite(dd$log_bmi) & is.finite(dd$age_c) & !is.na(dd$sex)
  if ("mega_pool" %in% names(dd)) keep_fg <- keep_fg & is.finite(dd$mega_pool)
  n_drop_fg <- sum(!keep_fg); if (n_drop_fg > 0L && isTRUE(verbose)) message("Fine–Gray: dropping ", n_drop_fg, " row(s) with missing predictors before fitting.")
  dd <- dd[keep_fg, , drop = FALSE]
  if (all(c("start","stop") %in% names(dd))) {
    dd <- filter_valid_intervals(dd, verbose = verbose, label = "Fine–Gray")
  } else {
    ok <- is.finite(dd$stop) & dd$stop > 0 & !is.na(dd$status)
    n_bad <- sum(!ok)
    if (n_bad > 0L && isTRUE(verbose))
      message("Fine–Gray: dropping ", n_bad, " row(s) with nonpositive or missing stop/status.")
    dd <- dd[ok, , drop = FALSE]
  }
  
  n_evt <- sum(dd$status == event_code)
  if (n_evt < 5) {
    if (isTRUE(verbose)) message("Fine–Gray: < 5 events (n = ", n_evt, ") – model skipped.")
    return(NULL)
  }
  
  use_quart <- isTRUE(spec$egfr_as_quartiles)
  egfr_df   <- spec$egfr_spline_df %||% df_spline
  
  if (use_quart) {
    qs <- stats::quantile(dd$egfr, probs = c(.25, .5, .75), na.rm = TRUE, names = FALSE)
    if (any(!is.finite(qs)) || any(diff(qs) <= 0)) {
      warning("Fine–Gray: cannot form eGFR quartiles (non-monotone breaks); using spline instead.",
              call. = FALSE, immediate. = TRUE)
      use_quart <- FALSE
    } else {
      dd$egfr_q <- cut(
        dd$egfr,
        breaks = c(-Inf, qs, Inf),
        labels = c("Q1","Q2","Q3","Q4"),
        include.lowest = TRUE,
        right = TRUE,
        ordered_result = TRUE
      )
    }
  }
  if (!use_quart) {
    spl_egfr <- splines::ns(dd$log_egfr, df = egfr_df)
    colnames(spl_egfr) <- paste0("log_egfr_s", seq_len(ncol(spl_egfr)))
    dd <- cbind(dd, spl_egfr)
  }
  
  
  spl_age  <- splines::ns(dd$age_c, df = spec$age_df)
  colnames(spl_age)  <- paste0("age_c_s", seq_len(ncol(spl_age)))
  dd <- cbind(dd, spl_age)
  
  if (isTRUE(spec$bmi_spline)) {
    spl_bmi <- splines::ns(dd$log_bmi, df = 3)
    colnames(spl_bmi) <- paste0("log_bmi_s", seq_len(ncol(spl_bmi)))
    dd <- cbind(dd, spl_bmi)
    bmi_terms <- colnames(spl_bmi)
  } else {
    bmi_terms <- "log_bmi"
  }
  
  rhs <- c(
    if (use_quart) "egfr_q" else colnames(spl_egfr),
    bmi_terms,
    "sex",
    colnames(spl_age)
  )
  if (isTRUE(spec$include_mega_pool) && any(dd$mega_pool == 1, na.rm = TRUE)) rhs <- c(rhs, "mega_pool")
  if (isTRUE(spec$include_cr_resid) && "cr_resid" %in% names(dd)) rhs <- c(rhs, "cr_resid")
  
  has_censor <- any(dd$status == 0)
  use_FGR <- switch(engine,
                    "FGR"  = TRUE,
                    "crr"  = FALSE,
                    "auto" = has_censor)
  
  if (use_FGR && !requireNamespace("riskRegression", quietly = TRUE)) {
    if (isTRUE(verbose)) message("riskRegression not available – falling back to crr().")
    use_FGR <- FALSE
  }
  if (use_FGR && !requireNamespace("prodlim", quietly = TRUE)) {
    stop("prodlim package is required for the FGR engine.", call. = FALSE)
  }
  if (!use_FGR && !requireNamespace("cmprsk", quietly = TRUE)) {
    stop("cmprsk package is required for crr() engine but is not installed.", call. = FALSE)
  }
  
  if (use_FGR) {
    form <- stats::as.formula(paste(
      "prodlim::Hist(stop, status) ~", paste(rhs, collapse = " + ")
    ))
    obj <- tryCatch(
      riskRegression::FGR(
        formula = form,
        data    = dd,
        cause   = event_code,
        CV      = cv_fold,
        cluster = if (patient_cluster) dd$patient else NULL
      ),
      error = function(e) { if (isTRUE(verbose)) message("Fine–Gray (FGR) fit failed: ", e$message); NULL }
    )
    engine_used <- "FGR"
    
    if (!is.null(obj)) {
      attr(obj, "age_mean")      <- age_mu
      attr(obj, "fg_engine")     <- engine_used
      attr(obj, "df_spline")     <- egfr_df
      if (!use_quart) {
        attr(obj, "egfr_knots_int")<- attr(splines::ns(dd$log_egfr, df = egfr_df), "knots")
        attr(obj, "egfr_knots_bd") <- attr(splines::ns(dd$log_egfr, df = egfr_df), "Boundary.knots")
      } else {
        attr(obj, "egfr_quartiles") <- qs
      }
      attr(obj, "age_knots_int") <- attr(spl_age, "knots")
      attr(obj, "age_knots_bd")  <- attr(spl_age, "Boundary.knots")
    }
    return(obj)
  }
  
  if (patient_cluster) {
    warning("fit_finegray(): patient_cluster ignored for crr(); use engine='FGR' for clustering.", call. = FALSE, immediate. = TRUE)
  }
  
  contr <- list(sex = "contr.treatment")
  if (use_quart && "egfr_q" %in% names(dd)) contr$egfr_q <- "contr.treatment"
  
  mm <- stats::model.matrix(
    stats::as.formula(paste("~", paste(rhs, collapse = " + "))),
    data = dd,
    contrasts.arg = contr
  )
  X <- mm[, -1, drop = FALSE]
  
  
  if (ncol(X) > 0) {
    nzv <- which(apply(X, 2, var, na.rm = TRUE) < 1e-12); if (length(nzv)) X <- X[, -nzv, drop = FALSE]
    if (ncol(X) > 1) {
      dup <- which(duplicated(t(X))); if (length(dup)) X <- X[, -dup, drop = FALSE]
    }
  }
  if (ncol(X) == 0) {
    warning("Design matrix became empty after cleaning – Fine–Gray skipped.", call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
  
  obj <- tryCatch(
    cmprsk::crr(
      ftime    = dd$stop,
      fstatus  = dd$status,
      cov1     = X,
      failcode = event_code,
      cencode  = 0
    ),
    error = function(e) { if (isTRUE(verbose)) message("Fine–Gray (crr) fit failed: ", e$message); NULL }
  )
  engine_used <- "crr"
  
  if (!is.null(obj)) {
    attr(obj, "age_mean")      <- age_mu
    attr(obj, "fg_engine")     <- engine_used
    attr(obj, "df_spline")     <- egfr_df
    if (!use_quart) {
      tmp_ns <- splines::ns(dd$log_egfr, df = egfr_df)
      attr(obj, "egfr_knots_int")<- attr(tmp_ns, "knots")
      attr(obj, "egfr_knots_bd") <- attr(tmp_ns, "Boundary.knots")
    } else {
      attr(obj, "egfr_quartiles") <- qs
    }
    attr(obj, "age_knots_int") <- attr(spl_age, "knots")
    attr(obj, "age_knots_bd")  <- attr(spl_age, "Boundary.knots")
  }
  return(obj)
}

fit_bayesian <- function(cp_data,
                         cause_codes = c(1L, 2L),
                         n_segments = 5,
                         min_events = 10,
                         inla_cores = getOption("inla_cores", 8),
                         pc_scale = getOption("inla_pc_scale", 1),
                         patient_cluster = FALSE,
                         spec = NULL) {
  
  spec <- .resolve_spec(spec)
  
  avail_mem_gb <- tryCatch({
    smi <- memuse::Sys.meminfo()
    tot <- if (is.numeric(smi)) {
      smi
    } else if ("totalram" %in% slotNames(smi)) {
      slot(smi, "totalram")
    } else {
      NA_real_
    }
    as.numeric(tot) / 1e9
  }, error = function(e) NA_real_)
  
  cores_by_ram <- if (is.finite(avail_mem_gb)) max(1, floor(avail_mem_gb / 8)) else inla_cores
  n_cores <- max(1, min(as.integer(inla_cores[1]), cores_by_ram))
  thread_arg <- sprintf("%d:%d", n_cores, n_cores)
  message("ℹ INLA threading: using ", n_cores, " core(s).")
  
  cuts <- choose_intervals(cp_data, n_max = n_segments, min_events = min_events)
  
  fit_one <- function(cause_code) {
    stage <- paste0("INLA (cause ", cause_code, ")")
    t0 <- log_timed(stage)
    cancel_hb <- schedule_heartbeat(sprintf("%s – still fitting …", stage), every = 120)
    on.exit(cancel_hb(), add = TRUE)
    
    age_mu <- mean(cp_data$age, na.rm = TRUE)
    dd <- cp_data %>%
      mutate(
        event     = (status == cause_code),
        log_egfr  = log(egfr),
        log_bmi   = log(bmi),
        age_c     = age - age_mu,
        mega_pool = as.integer(study_grp == "Merged-small-sites"),
        patient   = factor(patient)
      ) %>%
      filter(stop > start, is.finite(start), is.finite(stop))
    
    spl_raw <- ns(dd$age_c, df = spec$age_df)
    knots_int <- attr(spl_raw, "knots")
    knots_bd  <- attr(spl_raw, "Boundary.knots")
    spl <- as.data.frame(spl_raw)
    names(spl) <- paste0("age_spline", seq_len(ncol(spl)))
    dd <- bind_cols(dd, spl)
    
    fixed_terms <- c("log_egfr", "log_bmi", "sex",
                     paste0("age_spline", seq_len(ncol(spl))))
    
    if (isTRUE(spec$include_mega_pool) && any(dd$mega_pool == 1)) {
      fixed_terms <- c(fixed_terms, "mega_pool")
    }
    if (isTRUE(spec$include_cr_resid) && "cr_resid" %in% names(dd)) {
      fixed_terms <- c(fixed_terms, "cr_resid")
    }
    
    fixed_part <- paste(fixed_terms, collapse = " + ")
    
    re_patient <- if (patient_cluster) {
      sprintf(
        " + f(patient, model = 'iid',
          hyper = list(prec = list(prior = 'pc.prec', param = c(%f,0.01))))",
        pc_scale
      )
    } else ""
    
    re_centre <- if (isTRUE(spec$include_centre_re)) {
      sprintf(
        " + f(study_grp, model = 'iid',
          hyper = list(prec = list(prior = 'pc.prec', param = c(%f, 0.01))))",
        pc_scale
      )
    } else ""
    
    formula <- as.formula(sprintf(
      "event ~ -1 +
     f(interval, model = 'rw2',
       hyper = list(prec = list(prior = 'pc.prec', param = c(%f, 0.01)))) +
     %s%s%s",
      pc_scale, fixed_part, re_centre, re_patient
    ))
    
    dd_long <- survival::survSplit(
      Surv(start, stop, event) ~ .,
      data = dd, cut = cuts[-1], episode = "interval",
      start = "tstart", end = "tstop", event = "event"
    ) %>%
      mutate(
        exposure = tstop - tstart,
        interval = as.integer(interval),
        sex      = factor(sex)
      )
    bad_exp <- !is.finite(dd_long$exposure) | dd_long$exposure <= 0
    if (any(bad_exp)) dd_long <- dd_long[!bad_exp, ]
    
    res <- tryCatch(
      INLA::inla(
        formula,
        family = "poisson",
        data = dd_long, E = exposure,
        num.threads = thread_arg,
        verbose = inla_verb,
        control.fixed = list(prec = 1 / 25, expand.factor.strategy = "model.matrix"),
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE)
      ),
      error = function(e) {
        message("INLA fit failed (cause ", cause_code, "): ", e$message)
        NULL
      }
    )
    
    log_timed(stage, t0); cancel_hb()
    if (is.null(res)) return(NULL)
    
    list(
      model = res,
      cuts = cuts,
      age_knots_int = knots_int,
      age_knots_bd = knots_bd,
      age_mean = age_mu
    )
  }
  
  cause_fits <- lapply(cause_codes, fit_one)
  names(cause_fits) <- dplyr::recode(as.character(cause_codes), `1` = "death", `2` = "discharge")
  cause_fits
}
