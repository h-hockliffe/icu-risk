suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})


ensure_nonempty <- function(obj,
                            name = deparse(substitute(obj))) {
  if (
    is.null(obj) ||
      (is.data.frame(obj) && nrow(obj) == 0) ||
      (is.vector(obj) && all(is.na(obj))) ||
      (is.list(obj) && all(vapply(obj, is.null, logical(1))))
  ) {
    stop(name, " is empty or all NA – aborting", call. = FALSE)
  }
}


overall_epv_check <- function(df,
                              event_var = "icu_death",
                              n_par,
                              min_epp = 10) {
  ensure_nonempty(df, "df")

  if (!event_var %in% names(df)) {
    stop("`", event_var, "` not found in `df`.", call. = FALSE)
  }

  n_events <- sum(df[[event_var]] == 1L, na.rm = TRUE)
  epv <- n_events / n_par

  if (epv < min_epp) {
    stop(
      sprintf(
        "Overall EPV = %.1f (< %d).  Expand data or simplify the model.",
        epv, min_epp
      ),
      call. = FALSE
    )
  }

  message(sprintf("✓ Overall EPV = %.1f (≥ %d) – OK", epv, min_epp))
  invisible(epv)
}


epp_guardrail <- function(cp_list,
                          min_epp = 5, 
                          n_par,
                          fail_hard = FALSE) { 

  ensure_nonempty(cp_list, "cp_list")

  summarise_imp <- function(dat, imp_name) {
    dat %>%
      distinct(patient, .keep_all = TRUE) %>% 
      count(status, name = "events") %>%
      filter(status == 1L) %>% 
      transmute(
        .imp = imp_name,
        events = events,
        epp = events / n_par
      )
  }

  summary_tbl <- imap_dfr(cp_list, summarise_imp)

  violate <- summary_tbl %>% filter(epp < min_epp)

  if (nrow(violate) > 0) {
    msg <- paste0(
      "Imputation-level EPV below ", min_epp, ":\n",
      paste(sprintf("  • %s: %.2f", violate$.imp, violate$epp),
        collapse = "\n"
      )
    )

    if (fail_hard) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, immediate. = TRUE, call. = FALSE)
    }
  } else {
    message("✓ All imputations ≥ ", min_epp, " EPV – OK")
  }

  invisible(structure(summary_tbl,
    class = c("tbl_df", "guardrail_summary")
  ))
}

# --------------------------------------------------------------------
# Per-cause EP(P) guard + spec simplification helpers
# --------------------------------------------------------------------

#' Compute events-per-parameter (EPP) separately for each cause
#'
#' @param cp_data counting-process data.frame with one row per admission
#'        (or per distinct patient), columns: patient, status
#' @param n_par   number of fixed-effect parameters in the linear predictor
#' @param min_epp threshold for adequate events-per-parameter (default 10)
#' @return tibble with rows for "death" (1) and "discharge" (2)
epp_guardrail_per_cause <- function(cp_data, n_par, min_epp = 10) {
  stopifnot(is.data.frame(cp_data), "patient" %in% names(cp_data), "status" %in% names(cp_data))
  dd <- dplyr::distinct(cp_data, patient, .keep_all = TRUE)
  e1 <- sum(dd$status == 1L, na.rm = TRUE)
  e2 <- sum(dd$status == 2L, na.rm = TRUE)
  tibble::tibble(
    cause  = c("death", "discharge"),
    code   = c(1L, 2L),
    events = c(e1,  e2),
    epp    = c(e1 / n_par, e2 / n_par),
    ok     = c(e1 / n_par >= min_epp, e2 / n_par >= min_epp)
  )
}

simplify_spec <- function(spec) {
  stopifnot(is.list(spec))
  spec$age_df <- spec$age_df %||% 3
  spec$bmi_spline <- isTRUE(spec$bmi_spline)
  spec$include_mega_pool <- isTRUE(spec$include_mega_pool)
  if (isTRUE(spec$bmi_spline)) {
    spec$bmi_spline <- FALSE
    return(spec)
  }
  if (isTRUE(spec$include_mega_pool)) {
    spec$include_mega_pool <- FALSE
    return(spec)
  }
  if (!is.null(spec$age_df) && spec$age_df > 2) {
    spec$age_df <- spec$age_df - 1
    return(spec)
  }
  spec
}


#' @param cp_data counting-process data.frame
#' @param base_spec initial spec list (age_df, bmi_spline, include_mega_pool, include_cr_resid)
#' @param n_par number of parameters *with the current spec* (caller supplies)
#' @param min_epp target events-per-parameter (default 10)
#' @param max_steps safety cap on simplification steps
#' @return list(spec=final_spec, steps=log_tibble, ok=TRUE/FALSE, summary=per-cause tibble)
auto_simplify_to_epp <- function(cp_data, base_spec, n_par, min_epp = 10, max_steps = 3) {
  stopifnot(is.list(base_spec))
  steps <- list()
  spec  <- base_spec
  summ  <- epp_guardrail_per_cause(cp_data, n_par = n_par, min_epp = min_epp)
  if (all(summ$ok, na.rm = TRUE)) {
    return(list(spec = spec,
                steps = tibble::tibble(step = integer(0), action = character(0), n_par = numeric(0)),
                ok = TRUE,
                summary = summ))
  }
  for (i in seq_len(max_steps)) {
    prev <- spec
    spec <- simplify_spec(spec)
    action <- if (!identical(prev$bmi_spline, spec$bmi_spline)) {
      if (!spec$bmi_spline) "drop_bmi_spline" else "add_bmi_spline"
    } else if (!identical(prev$include_mega_pool, spec$include_mega_pool)) {
      if (!spec$include_mega_pool) "drop_mega_pool" else "add_mega_pool"
    } else if (!identical(prev$age_df, spec$age_df)) {
      paste0("set_age_df_", spec$age_df)
    } else {
      "no_change"
    }
    steps[[length(steps) + 1L]] <- tibble::tibble(step = i, action = action, n_par = n_par)
    summ <- epp_guardrail_per_cause(cp_data, n_par = n_par, min_epp = min_epp)
    if (all(summ$ok, na.rm = TRUE)) {
      return(list(spec = spec, steps = dplyr::bind_rows(steps), ok = TRUE, summary = summ))
    }
  }
  list(spec = spec, steps = dplyr::bind_rows(steps), ok = FALSE, summary = summ)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_tar_load <- function(names, store = "_targets") {
  present <- names[targets::tar_exist_objects(names, store = store)]
  if (!length(present)) return(invisible(NULL))
  targets::tar_load(present, store = store)
}
