suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(splines)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

make_spec <- function(age_df = 3,
                      bmi_spline = TRUE,
                      include_mega_pool = TRUE,
                      include_cr_resid = FALSE,
                      egfr_spline_df = NULL,
                      egfr_as_quartiles = FALSE,
                      include_centre_re = TRUE,
                      age_tvc = NULL) {          
  list(
    age_df = age_df,
    bmi_spline = isTRUE(bmi_spline),
    include_mega_pool = isTRUE(include_mega_pool),
    include_cr_resid = isTRUE(include_cr_resid),
    egfr_spline_df = egfr_spline_df,
    egfr_as_quartiles = isTRUE(egfr_as_quartiles),
    include_centre_re = isTRUE(include_centre_re),
    age_tvc = age_tvc                           
  )
}




compute_cr_resid <- function(df, age_df = 3) {
  need <- c("creat", "egfr", "age", "sex")
  if (!all(need %in% names(df))) {
    warning("compute_cr_resid(): required columns missing; returning NA.")
    return(rep(NA_real_, nrow(df)))
  }
  if (!is.factor(df$sex)) df$sex <- factor(df$sex)
  
  ok <- is.finite(df$creat) & is.finite(df$egfr) & is.finite(df$age) & !is.na(df$sex)
  if (sum(ok) < 30L) {
    warning("compute_cr_resid(): < 30 complete rows; returning NA.")
    return(rep(NA_real_, nrow(df)))
  }
  dat <- df[ok, , drop = FALSE]
  if (any(dat$creat <= 0, na.rm = TRUE) || any(dat$egfr <= 0, na.rm = TRUE)) {
    warning("compute_cr_resid(): non-positive value(s); returning NA.")
    return(rep(NA_real_, nrow(df)))
  }
  fit <- try(lm(log(creat) ~ splines::ns(age, df = age_df) + sex + log(egfr), data = dat), silent = TRUE)
  if (inherits(fit, "try-error")) {
    warning("compute_cr_resid(): residual model failed; returning NA.")
    return(rep(NA_real_, nrow(df)))
  }
  r <- rep(NA_real_, nrow(df))
  r_ok <- try(as.numeric(stats::rstudent(fit)), silent = TRUE)
  if (inherits(r_ok, "try-error")) r_ok <- as.numeric(stats::residuals(fit))
  r[ok] <- r_ok
  r
}



scenario_exclude_ckd <- function(cp_list) {
  stopifnot(is.list(cp_list))
  purrr::imap(cp_list, function(dat, nm) {
    dplyr::filter(dat, is.finite(egfr) & egfr >= 60)
  })
}


scenario_exclude_covid <- function(cp_list, exclude_levels = NULL) {
  stopifnot(is.list(cp_list))
  purrr::imap(cp_list, function(dat, nm) {
    sg <- as.character(dat$study_grp)
    if (is.null(exclude_levels)) {
      keep <- !grepl("COVID", sg, ignore.case = TRUE)
      dplyr::filter(dat, keep)
    } else {
      dplyr::filter(dat, !(study_grp %in% exclude_levels))
    }
  })
}


scenario_renal_not_ordered <- function(cp_list, set_missing_to_zero = FALSE) {
  stopifnot(is.list(cp_list))
  purrr::imap(cp_list, function(dat, nm) {
    dat <- dat %>%
      mutate(
        renal_not_ordered = as.integer(is.na(egfr) | is.na(bmi)),
        egfr = if (isTRUE(set_missing_to_zero)) ifelse(is.na(egfr), 0, egfr) else egfr,
        bmi  = if (isTRUE(set_missing_to_zero)) ifelse(is.na(bmi),  0, bmi)  else bmi
      )
    dat
  })
}

scenario_add_cr_resid <- function(cp_list, age_df = 3) {
  stopifnot(is.list(cp_list))
  purrr::imap(cp_list, function(dat, nm) {
    need <- c("creat", "egfr", "age", "sex")
    if (!all(need %in% names(dat))) {
      warning("scenario_add_cr_resid(): missing inputs in ", nm, " — skipping.")
      return(dat)
    }
    dat$cr_resid <- compute_cr_resid(dat, age_df = age_df)
    if (!any(is.finite(dat$cr_resid))) {
      warning("scenario_add_cr_resid(): all-NA residuals in ", nm, " — dropping cr_resid.")
      dat$cr_resid <- NULL
    }
    dat
  })
}



scenario_creatinine_umol <- function(cp_list) {
  stopifnot(is.list(cp_list))
  purrr::imap(cp_list, function(dat, nm) {
    if ("creat" %in% names(dat)) {
      dat$creat_umol <- ifelse(is.finite(dat$creat), dat$creat * 88.4, NA_real_)
    }
    dat
  })
}

spec_egfr_spline5 <- function(base_spec = make_spec()) {
  modifyList(base_spec, list(egfr_spline_df = 5L))
}
spec_egfr_quartiles <- function(base_spec = make_spec()) {
  modifyList(base_spec, list(egfr_as_quartiles = TRUE))
}
spec_no_centre_re <- function(base_spec = make_spec()) {
  modifyList(base_spec, list(include_centre_re = FALSE))
}

spec_drop_mega_pool_covariate <- function(base_spec = make_spec()) {
  modifyList(base_spec, list(include_mega_pool = FALSE))
}


build_sensitivity_inputs <- function(cp_list,
                                     base_spec = make_spec(),
                                     covid_levels = NULL) {
  stopifnot(is.list(cp_list))
  
  out <- list(
    primary = list(
      name   = "primary",
      cp_list = cp_list,
      spec   = base_spec,
      notes  = "Primary model: log-eGFR + sex + age spline(df=3) + log-BMI; centre RE as specified."
    ),
    
    excl_ckd = list(
      name   = "exclude_ckd",
      cp_list = scenario_exclude_ckd(cp_list),
      spec   = base_spec,
      notes  = "Sensitivity: exclude baseline CKD (eGFR < 60 mL/min/1.73 m^2)."
    ),
    
    excl_covid = list(
      name   = "exclude_covid",
      cp_list = scenario_exclude_covid(cp_list, exclude_levels = covid_levels),
      spec   = base_spec,
      notes  = "Sensitivity: exclude COVID-19 sub-studies (epoch/selection effects)."
    ),
    
    renal_not_ordered = list(
      name   = "renal_not_ordered",
      cp_list = scenario_renal_not_ordered(cp_list, set_missing_to_zero = FALSE),
      spec   = base_spec,
      notes  = "Sensitivity: include binary 'renal_not_ordered' indicator."
    ),
    
    creat_umol = list(
      name   = "creatinine_umol_flag",
      cp_list = scenario_creatinine_umol(cp_list),
      spec   = base_spec,
      notes  = "Sensitivity: create creatinine µmol/L companion variable for unit-conversion checks."
    ),
    
    egfr_spline5 = list(
      name   = "egfr_spline_df5",
      cp_list = cp_list,
      spec   = spec_egfr_spline5(base_spec),
      notes  = "Sensitivity: replace log-linear eGFR with restricted cubic spline (df=5)."
    ),
    
    egfr_quartiles = list(
      name   = "egfr_quartiles",
      cp_list = cp_list,
      spec   = spec_egfr_quartiles(base_spec),
      notes  = "Sensitivity: categorize eGFR into quartiles (reference: highest quartile)."
    ),
    
    added_value_creatinine = list(
      name   = "added_value_cr_resid",
      cp_list = scenario_add_cr_resid(cp_list, age_df = base_spec$age_df %||% 3),
      spec   = modifyList(base_spec, list(include_cr_resid = TRUE)),
      notes  = "Added-value: include creatinine residual orthogonal to eGFR, age, and sex."
    ),
    excl_mega_pool = list(
      name   = "exclude_merged_small_sites",
      cp_list = scenario_exclude_mega_pool(cp_list, level = "Merged-small-sites"),
      spec   = base_spec,
      notes  = "Sensitivity: drop the pooled 'Merged-small-sites' cohort entirely."
    ),
    no_mega_pool_cov = list(
      name   = "no_mega_pool_covariate",
      cp_list = cp_list,
      spec   = spec_drop_mega_pool_covariate(base_spec),
      notes  = "Sensitivity (spec): exclude the mega_pool fixed-effect indicator."
    ),
    no_centre_re = list(
      name   = "no_centre_random_effect",
      cp_list = cp_list,
      spec   = spec_no_centre_re(base_spec),
      notes  = "Sensitivity (spec): fit models without centre-level random effects."
    )
    
  )
  
  out
}


compare_metric_bundles <- function(base_metrics, alt_metrics) {
  if (is.null(base_metrics) || is.null(alt_metrics)) {
    return(tibble(
      model = character(0),
      delta_DIC = numeric(0),
      delta_AIC = numeric(0),
      delta_BIC = numeric(0),
      delta_WAIC = numeric(0)
    ))
  }
  b <- try(base_metrics$pooled, silent = TRUE)
  a <- try(alt_metrics$pooled, silent = TRUE)
  if (inherits(b, "try-error") || inherits(a, "try-error") || is.null(b) || is.null(a)) {
    return(tibble(model=character(), delta_DIC=numeric(), delta_AIC=numeric(),
                  delta_BIC=numeric(), delta_WAIC=numeric()))
  }
  m <- dplyr::inner_join(
    dplyr::select(b, model, DIC, AIC, BIC, WAIC),
    dplyr::select(a, model, DIC, AIC, BIC, WAIC),
    by = "model",
    suffix = c(".base", ".alt")
  )
  if (!nrow(m)) {
    return(tibble(
      model = character(0),
      delta_DIC = numeric(0),
      delta_AIC = numeric(0),
      delta_BIC = numeric(0),
      delta_WAIC = numeric(0)
    ))
  }
  m %>%
    transmute(
      model,
      delta_DIC  = DIC.alt  - DIC.base,
      delta_AIC  = AIC.alt  - AIC.base,
      delta_BIC  = BIC.alt  - BIC.base,
      delta_WAIC = WAIC.alt - WAIC.base
    )
}


scenario_exclude_mega_pool <- function(cp_list, level = "Merged-small-sites") {
  stopifnot(is.list(cp_list))
  purrr::imap(cp_list, function(dat, nm) {
    dat <- dplyr::filter(dat, as.character(study_grp) != level)
    if (is.factor(dat$study_grp)) dat$study_grp <- droplevels(dat$study_grp)
    dat
  })
}

