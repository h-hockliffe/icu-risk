options(
  inla_cores        = min(parallel::detectCores(), 8L),
  inla_pc_scale     = 1,
  fg_verbose        = TRUE,
  mi_allow_parallel = TRUE
)

suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(qs2)
})

tar_source("R")

required_funs <- c(
  "clean_master", "overall_epv_check", "impute_data",
  "transform_counting_process", "epp_guardrail",
  "fit_bayesian", "fit_cox", "fit_finegray",
  "make_mcmc_dataset", "fit_jags_pwexp", "compare_inla_mcmc",
  "extract_model_metrics", "validate_models",
  "make_descriptive_tables", "plot_cumulative_incidence",
  "plot_decision_curve",
  "make_spec", "build_sensitivity_inputs",
  "compare_metric_bundles", "plot_calibration_models", "cox_ph_diagnostics",
  "make_coef_table", "compute_raw_metrics", "compute_pooled_metrics",
  "compute_group_metrics", "plot_cindex_td", "compute_cindex_td",
  "inla_mcmc_agreement"
)
check_helpers <- function(req = required_funs) {
  missing <- req[!vapply(req, exists, logical(1), mode = "function", inherits = TRUE)]
  if (length(missing)) {
    msg <- paste0(
      "Missing helper functions in R/: ",
      paste(missing, collapse = ", ")
    )
    if (identical(Sys.getenv("TAR_SKIP_GUARD"), "1")) {
      warning(msg, call. = FALSE)
    } else {
      stop(msg, call. = FALSE)
    }
  }
  invisible(TRUE)
}

check_helpers()

tar_option_set(
  seed     = 1000,
  envir    = globalenv(),
  packages = c(
    "dplyr", "tidyr", "purrr", "readxl", "lubridate",
    "mice", "INLA", "survival", "riskRegression", "timeROC",
    "splines", "gtsummary", "ggplot2", "targets", "tarchetypes",
    "qs2", "prodlim", "arrow", "frailtypack",
    "progressr", "later", "memuse", "stats",
    "rjags", "coda", "quarto",
    "cmprsk", "gridExtra", "patchwork"
  ),
  format   = "rds",
  resources = tar_resources(
    qs = tar_resources_qs(preset = "high")
  )
)

quarto_enabled <- function() identical(Sys.getenv("SKIP_QUARTO"), "1") == FALSE

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x)) return(default)
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("1","true","t","yes","y","on")
}

sens_enabled <- function() {
  opt <- getOption("enable_sens", NULL)
  if (!is.null(opt)) return(as_bool(opt))
  env <- Sys.getenv("ENABLE_SENS", unset = NA_character_)
  if (!is.na(env)) return(as_bool(env))
  TRUE
}

TIME_HORIZONS <- c(0.5, 1, 3, 7, 15, 30, 60, 90)

RAW_FILE <- {
  env <- Sys.getenv("RAW_FILE", unset = NA_character_)
  def <- "G:/Drive partagés/OI - DSI - Thesis Hockliffe/FINAL BAYESIAN_07042025 (1).xlsx"
  p   <- if (!is.na(env) && nzchar(env)) env else def
  n   <- try(normalizePath(p, winslash = "/", mustWork = FALSE), silent = TRUE)
  if (!file.exists(n)) {
    stop("\n✗ Raw workbook not found:\n   ", n,
         "\nMount the share or set RAW_FILE before running.")
  }
  n
}

list(
  
  tar_target(val_B,        200L),
  tar_target(val_horizon,  TIME_HORIZONS),
  tar_target(mi_m,         20L),
  
  tar_target(run_sens, sens_enabled()),
  tar_target(raw_file, RAW_FILE, format = "file"),
  
  tar_target(clean, clean_master(raw = raw_file), format = "rds"),
  tar_target(
    clean_cluster,
    clean_master(raw = raw_file, deduplicate = FALSE),
    format = "rds"
  ),
  
  tar_target(
    overall_epv,
    overall_epv_check(clean, event_var = "icu_death", n_par = 12)
  ),
  
  tar_target(
    imps,
    { overall_epv; impute_data(clean, m = mi_m) },
    format = "qs"
  ),
  tar_target(
    imps_cluster,
    { impute_data(clean_cluster, m = mi_m) },
    format = "qs"
  ),
  
  tar_target(count, transform_counting_process(imps), format = "qs"),
  tar_target(
    count_cluster,
    transform_counting_process(imps_cluster, keep_duplicates = TRUE),
    format = "qs"
  ),
  
  tar_target(epp_check, epp_guardrail(count, min_epp = 5, n_par = 12)),
  
  tar_target(
    cp_branch,
    {
      if (is.list(count) && length(count) == 1L && is.data.frame(count[[1]])) {
        count[[1]]
      } else if (is.list(count) &&
                 length(count) == 1L && is.list(count[[1]]) &&
                 is.data.frame(count[[1]][[1]])) {
        count[[1]][[1]]
      } else if (is.data.frame(count)) {
        count
      } else {
        stop("cp_branch: unexpected object type; wanted a data.frame")
      }
    },
    pattern   = map(count),
    iteration = "list"
  ),
  
  tar_target(
    bayes_fit,
    { epp_check; fit_bayesian(cp_branch, cause_codes = c(1L, 2L)) },
    pattern = map(cp_branch), format = "rds"
  ),
  tar_target(
    cox_fit,
    { epp_check; fit_cox(cp_branch, spec = make_spec(bmi_spline = NULL, age_tvc = "log")) },
    pattern = map(cp_branch),
    format  = "rds"
  ),
  tar_target(
    fg_fit,
    { epp_check; fit_finegray(cp_branch) },
    pattern = map(cp_branch), format = "qs"
  ),
  
  tar_target(
    cp_branch_cluster,
    {
      if (is.list(count_cluster) && length(count_cluster) == 1L &&
          is.data.frame(count_cluster[[1]])) {
        count_cluster[[1]]
      } else if (is.data.frame(count_cluster)) {
        count_cluster
      } else {
        stop("cp_branch_cluster: unexpected object shape")
      }
    },
    pattern   = map(count_cluster),
    iteration = "list"
  ),
  
  tar_target(
    bayes_fit_cluster,
    fit_bayesian(cp_branch_cluster, patient_cluster = TRUE),
    pattern = map(cp_branch_cluster), format = "rds"
  ),
  
  tar_target(
    cox_fit_cluster,
    fit_cox(cp_branch_cluster, patient_cluster = TRUE,
            spec = make_spec(bmi_spline = NULL, age_tvc = "log")),
    pattern = map(cp_branch_cluster),
    format  = "rds"
  ),
  
  tar_target(
    fg_fit_cluster,
    {
      dat1 <- cp_branch_cluster
      has_cens <- any(dplyr::distinct(dat1, patient, .keep_all = TRUE)$status == 0L)
      if (!has_cens) {
        message("FG cluster: no censoring in data — skipping (NULL).")
        NULL
      } else {
        fit_finegray(dat1, engine = "FGR", patient_cluster = TRUE)
      }
    },
    pattern = map(cp_branch_cluster), format = "qs"
  ),
  
  tar_target(
    cp_branch_nosmall,
    {
      dat <- cp_branch
      dplyr::filter(dat, study_grp != "Merged-small-sites")
    },
    pattern   = map(cp_branch),
    iteration = "list"
  ),
  
  tar_target(
    bayes_fit_nosmall,
    { epp_guardrail(list(cp_branch_nosmall), min_epp = 5, n_par = 12)
      fit_bayesian(cp_branch_nosmall, cause_codes = c(1L, 2L)) },
    pattern   = map(cp_branch_nosmall),
    iteration = "list",
    format    = "rds"
  ),
  
  tar_target(
    cox_fit_nosmall,
    { epp_guardrail(list(cp_branch_nosmall), min_epp = 5, n_par = 12)
      fit_cox(cp_branch_nosmall, spec = make_spec(bmi_spline = NULL, age_tvc = "log")) },
    pattern   = map(cp_branch_nosmall),
    iteration = "list",
    format    = "rds"
  ),
  
  tar_target(
    fg_fit_nosmall,
    { epp_guardrail(list(cp_branch_nosmall), min_epp = 5, n_par = 12)
      fit_finegray(cp_branch_nosmall) },
    pattern   = map(cp_branch_nosmall),
    iteration = "list",
    format    = "qs"
  ),
  
  tar_target(
    val_nosmall_refit,
    validate_models(
      cp_branch_nosmall,
      bayes_fit_nosmall,
      cox_fit_nosmall,
      fg_fit_nosmall,
      B = val_B,                  
      time_horizon = 30           
    ),
    format = "qs"
  ),
  
  tar_target(
    cp_list, cp_branch,
    pattern   = map(cp_branch), iteration = "list"
  ),
  tar_target(
    bayes_list, bayes_fit,
    pattern   = map(bayes_fit), iteration = "list"
  ),
  tar_target(
    cox_list, cox_fit,
    pattern   = map(cox_fit), iteration = "list"
  ),
  tar_target(
    fg_list, fg_fit,
    pattern   = map(fg_fit), iteration = "list"
  ),
  
  tar_target(coef_table, make_coef_table(bayes_list, cox_list, fg_list)),
  
  tar_target(
    cp_list_cluster, cp_branch_cluster,
    pattern = map(cp_branch_cluster), iteration = "list"
  ),
  tar_target(
    bayes_list_cluster, bayes_fit_cluster,
    pattern = map(bayes_fit_cluster), iteration = "list"
  ),
  tar_target(
    cox_list_cluster, cox_fit_cluster,
    pattern = map(cox_fit_cluster), iteration = "list"
  ),
  tar_target(
    fg_list_cluster, fg_fit_cluster,
    pattern = map(fg_fit_cluster), iteration = "list"
  ),
  
  tar_target(tables,      make_descriptive_tables(clean)),

  tar_target(
    fig_cif,
    plot_cumulative_incidence(
      bayes_list, cp_list,
      annotate_zoom_at   = 7,          
      annotate_full_at   = numeric(0), 
      show_cif_labels    = FALSE,      
      show_markers_full  = FALSE, 
    ),
    format = "rds"
  ),
  
  tar_target(
    fig_calib,
    plot_calibration_models(
      bayes_list, cox_list, fg_list, cp_list,
      horizon = val_horizon
    ),
    format = "rds"
  ),
  
  tar_target(
    fig_decision,
    setNames(
      lapply(val_horizon, function(h)
        plot_decision_curve(
          bayes_list, cox_list, fg_list, cp_list,
          time_horizon = h
        )
      ),
      paste0(val_horizon, "d")
    ),
    format = "rds"
  ),
  
  
  tar_target(
    model_metrics,
    extract_model_metrics(bayes_list, cox_list, fg_list)
  ),
  tar_target(
    raw_metrics,
    compute_raw_metrics(cp_list, bayes_list, cox_list, fg_list, horizons = val_horizon),
    format = "qs"
  ),
  tar_target(
    pooled_metrics,
    compute_pooled_metrics(raw_metrics),
    format = "qs"
  ),
  
  tar_target(cindex_td, compute_cindex_td(cp_list, bayes_list, cox_list, fg_list), format = "qs"),
  tar_target(fig_c_td,  plot_cindex_td(cindex_td)),
  
  tar_target(
    val_allh,
    validate_models(
      cp_list,
      bayes_list,
      cox_list,
      fg_list,
      B = val_B,
      time_horizon = val_horizon
    ),
    format = "qs"
  ),
  
  tar_target(
    group_metrics,
    compute_group_metrics(cp_list, bayes_list, cox_list, fg_list, horizons = val_horizon),
    format = "qs"
  ),
  
  tar_target(
    model_metrics_cluster,
    extract_model_metrics(
      bayes_list_cluster,
      cox_list_cluster,
      fg_list_cluster,
      strict = FALSE,
      context_label = "clustered branch"
    ),
    format = "qs"
  ),
  tar_target(
    val_results_cluster,
    validate_models(
      cp_list_cluster,
      bayes_list_cluster,
      cox_list_cluster,
      fg_list_cluster,
      B = val_B,
      time_horizon = val_horizon
    ),
    format = "qs"
  ),
  
  tar_target(
    val_no_small,
    validate_models(
      cp_list %>% purrr::map(~ dplyr::filter(.x, study_grp != "Merged-small-sites")),
      bayes_list, cox_list, fg_list,
      B = val_B,
      time_horizon = 30                  
    ),
    format = "qs"
  ),
  
  tar_target(cox_ph, cox_ph_diagnostics(cox_list), format = "rds"),
  
  tar_target(
    val_ckdneg,
    validate_models(
      cp_list %>% purrr::map(~ dplyr::filter(.x, egfr >= 60)),
      bayes_list, cox_list, fg_list,
      B = val_B,
      time_horizon = 30                  
    ),
    format = "qs"
  ),
  
  tar_target(
    sens_inputs,
    if (run_sens) {
      build_sensitivity_inputs(
        cp_list,
        base_spec = make_spec(
          age_df            = 3,
          bmi_spline        = TRUE,
          include_mega_pool = TRUE,
          include_cr_resid  = FALSE,
          include_centre_re = TRUE
        ),
        covid_levels = NULL
      )
    } else {
      NULL
    },
    format = "qs"
  ),
  
  tar_target(
    sens_fits,
    {
      x <- sens_inputs
      if (is.null(x)) return(NULL)
      
      purrr::imap(x, function(scn, nm) {
        cp <- scn$cp_list
        stopifnot(is.list(cp))
        sp <- scn$spec
        
        saf <- function(expr, tag) {
          tryCatch(expr, error = function(e) {
            warning(sprintf("Sensitivity '%s' – %s failed: %s", nm, tag, conditionMessage(e)))
            NULL
          })
        }
        
        bay <- purrr::map(cp, ~ saf(fit_bayesian(.x, spec = sp), "bayes"))
        cox <- purrr::map(cp, ~ saf(fit_cox(.x, spec = sp), "cox"))
        fg  <- purrr::map(cp, ~ saf({
          has_cens <- any(dplyr::distinct(.x, patient, .keep_all = TRUE)$status == 0L)
          if (!has_cens) {
            message("Sensitivity '", nm, "': no censoring — FG set to NULL.")
            NULL
          } else {
            fit_finegray(.x, spec = sp)
          }
        }, "finegray"))
        
        list(name = nm, notes = scn$notes, bay = bay, cox = cox, fg = fg)
      })
    },
    format = "qs"
  ),
  
  tar_target(
    sens_metrics,
    {
      if (is.null(sens_fits)) return(NULL)
      purrr::imap(sens_fits, function(sf, nm) {
        mb <- extract_model_metrics(sf$bay, sf$cox, sf$fg, strict = FALSE,
                                    context_label = paste0("sensitivity:", nm))
        list(name = nm, notes = sf$notes, metrics = mb)
      })
    },
    format = "qs"
  ),
  
  tar_target(
    sens_deltas,
    {
      if (is.null(sens_metrics) || is.null(model_metrics)) return(NULL)
      purrr::imap(sens_metrics, function(sm, nm) {
        del <- compare_metric_bundles(model_metrics, sm$metrics)
        list(name = nm, notes = sm$notes, deltas = del)
      })
    },
    format = "qs"
  ),
  

  tar_target(
    mcmc_agree_branch,
    inla_mcmc_agreement(
      bayes_list$death,   
      cp_list,            
      n_segments = 3,
      n_burn = 4000, n_iter = 6000, n_chains = 4,
      tau_fixed = 0.25, tau_age = 0.25,
      seed = 123L
    ),
    pattern = map(bayes_list, cp_list),
    iteration = "list",
    format = "qs"
  ),
  
  
  tar_target(
    mcmc_agree_summary,
    {
      
      br <- mcmc_agree_branch
      stopifnot(length(br) == length(bayes_list))
      br <- purrr::map2(br, bayes_list, ~{
        obj <- .x
        obj$summary <- fill_inla_cols(obj$summary, .y$death)
        obj
      })
      
      dplyr::bind_rows(purrr::imap(br, ~{
        tibble::tibble(
          branch   = .y,
          rhat_max = .x$convergence$rhat_max,
          ess_min  = .x$convergence$ess_min,
          ess_med  = .x$convergence$ess_median,
          r_total  = .x$agreement$r_total,
          r_no_age = .x$agreement$r_no_age
        )
      }))
    },
    format = "rds"
  ),
  
  tar_target(
    thesis_html,
    {
      invisible(list(
        val_B, val_horizon, mi_m,
        cp_list, bayes_list, cox_list, fg_list,
        val_allh, val_no_small, val_ckdneg, val_nosmall_refit, val_results_cluster,
        tables, fig_cif, fig_calib, fig_decision,
        pooled_metrics, cindex_td, fig_c_td, group_metrics,
        model_metrics, model_metrics_cluster, coef_table, cox_ph,
        sens_metrics, sens_deltas
      ))
      
      dir.create("_targets/logs", recursive = TRUE, showWarnings = FALSE)
      logfile <- normalizePath("_targets/logs/quarto_thesis.log", winslash = "/", mustWork = FALSE)
      if (file.exists(logfile)) file.remove(logfile)
      
      old_tmpdir <- Sys.getenv("TMPDIR", unset = NA_character_)
      on.exit({
        if (is.na(old_tmpdir)) Sys.unsetenv("TMPDIR") else Sys.setenv(TMPDIR = old_tmpdir)
      }, add = TRUE)
      Sys.unsetenv("TMPDIR")
      Sys.setenv(QUARTO_LOG_LEVEL = "DEBUG", QUARTO_PRINT_STACK = "1")
      
      qmd_in   <- "Thesis.qmd"
      out_name <- "thesis.html"
      
      zz <- file(logfile, open = "wt")
      on.exit(close(zz), add = TRUE)
      sink(zz, type = "output")
      sink(zz, type = "message", append = TRUE)
      
      ok <- TRUE
      tryCatch({
        quarto::quarto_render(
          input         = qmd_in,
          output_format = "html",
          output_file   = out_name,
          as_job        = FALSE,
          quiet         = FALSE
        )
      }, error = function(e) {
        ok <<- FALSE
        message("quarto_render() failed: ", conditionMessage(e))
      })
      
      sink(type = "message"); sink(type = "output")
      
      if (!ok || !file.exists(out_name)) {
        stop("Quarto render failed. See log: ", logfile)
      }
      
      normalizePath(out_name, winslash = "/", mustWork = TRUE)
    },
    format = "file",
    cue = tar_cue(mode = ifelse(quarto_enabled(), "thorough", "never"))
  )
  
  
)
