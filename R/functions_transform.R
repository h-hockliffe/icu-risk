suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(purrr)
  library(mice) 
  library(tibble)
})


choose_intervals <- function(cp_data,
                             n_max = getOption("pwexp_n_segments", 5),
                             min_events = 10,
                             max_iter = 20,
                             verbose = getOption("choose_int_verbose", TRUE)) {
  stopifnot(
    is.data.frame(cp_data),
    all(c("stop", "status") %in% names(cp_data))
  )
  
  evt <- cp_data[cp_data$status %in% c(1L, 2L), "stop", drop = TRUE]
  if (!length(evt) || sum(cp_data$status %in% c(1L, 2L)) < 2 * min_events) {
    cuts <- c(0, max(cp_data$stop, na.rm = TRUE))
    attr(cuts, "interval_counts") <- tibble(interval = factor(), death = integer(), discharge = integer(), good = logical())
    attr(cuts, "stopped_reason")  <- "degenerate"
    attr(cuts, "min_events")      <- min_events
    attr(cuts, "n_max")           <- n_max
    return(cuts)
  }
  evt <- sort(unique(evt))
  q <- stats::quantile(evt, probs = seq(0, 1, length.out = n_max + 1L), type = 2)
  cuts <- unique(as.numeric(q))
  cuts[1] <- 0
  cuts[length(cuts)] <- max(cp_data$stop, na.rm = TRUE)
  cuts <- sort(unique(cuts))
  
  count_ok <- function(cu) {
    labs <- cut(cp_data$stop, breaks = cu, include.lowest = TRUE, right = TRUE)
    tab  <- cp_data |>
      mutate(interval = labs) |>
      filter(status %in% c(1L, 2L)) |>
      count(interval, status) |>
      tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0, names_sort = TRUE) |>
      mutate(good = (`1` >= min_events) & (`2` >= min_events))
    list(ok = all(tab$good), tab = tab)
  }
  
  it <- 0L
  while (it < max_iter) {
    it <- it + 1L
    chk <- count_ok(cuts)
    if (chk$ok) {
      if (isTRUE(verbose)) message("✓ choose_intervals(): satisfied after ", it, " iteration(s).")
      out <- cuts
      attr(out, "interval_counts") <- chk$tab |>
        rename(death = `1`, discharge = `2`) |>
        mutate(interval = forcats::fct_inorder(interval))
      attr(out, "stopped_reason")  <- "ok"
      attr(out, "min_events")      <- min_events
      attr(out, "n_max")           <- n_max
      return(out)
    }
    tab <- chk$tab |> mutate(total = `1` + `2`)
    bad_idx <- which(!tab$good)
    if (!length(bad_idx)) break
    j <- bad_idx[which.min(tab$total[bad_idx])]
    left_total  <- if (j > 1)  tab$total[j - 1] else -Inf
    right_total <- if (j < nrow(tab)) tab$total[j + 1] else -Inf
    if (right_total >= left_total && j < length(cuts) - 1L) {
      cuts <- cuts[-(j + 1L)]
    } else if (j > 1L) {
      cuts <- cuts[-j]
    } else {
      break
    }
    cuts <- sort(unique(cuts))
    if (length(cuts) <= 2L) break
  }
  
  med <- stats::median(evt, na.rm = TRUE)
  cuts <- sort(unique(c(0, med, max(cp_data$stop, na.rm = TRUE))))
  chk  <- count_ok(cuts)
  out  <- cuts
  attr(out, "interval_counts") <- chk$tab |>
    rename(death = `1`, discharge = `2`) |>
    mutate(interval = forcats::fct_inorder(interval))
  attr(out, "stopped_reason")  <- "fallback_median"
  attr(out, "min_events")      <- min_events
  attr(out, "n_max")           <- n_max
  if (isTRUE(verbose)) message("choose_intervals(): used fallback split at median event time.")
  out
}



transform_counting_process <- function(
    imps,
    los_cap = 120,
    censor_date = as.POSIXct(
      "2023-12-31 23:59:59",
      tz = "UTC"
    ),
    keep_duplicates = FALSE
) {
  map(seq_len(imps$m), function(m) {
    dat <- mice::complete(imps, m)
    dat <- dat %>%
      dplyr::mutate(
        egfr = dplyr::coalesce(egfr, suppressWarnings(as.numeric(e_gfr))),
        egfr = pmin(pmax(egfr, 1), 120)   
      )

    if (all(c("creat", "egfr", "leeftijd_jaar", "sex") %in% names(dat))) {
      dd_lm <- dat %>%
        dplyr::mutate(
          log_cr   = suppressWarnings(log(creat)),
          log_egfr = suppressWarnings(log(egfr)),
          age_num  = suppressWarnings(as.numeric(leeftijd_jaar)),
          sex_fac  = if (is.factor(sex)) sex else factor(sex)
        )
      
      ok <- is.finite(dd_lm$log_cr) &
        is.finite(dd_lm$log_egfr) &
        is.finite(dd_lm$age_num) &
        !is.na(dd_lm$sex_fac)
      
      if (sum(ok) >= 30) {  
        fit <- stats::lm(log_cr ~ log_egfr + age_num + sex_fac, data = dd_lm[ok, , drop = FALSE])
        rstd <- rep(NA_real_, nrow(dat))
        rstd[ok] <- tryCatch(stats::rstudent(fit), error = function(e) NA_real_)
        dat$cr_resid <- rstd
      } else {
        dat$cr_resid <- NA_real_
      }
    } else {
      dat$cr_resid <- NA_real_
    }
    

    dat <- dat %>%
      mutate(
        status = case_when(
          icu_death == 1L ~ 1L, 
          icu_death == 0L ~ 2L, 
          TRUE ~ NA_integer_
        )
      ) %>%
      filter(!is.na(status)) 

    dat <- dat %>%
      transmute(
        .imp      = m,
        patient   = studie_patient_id,
        study_grp = study_grp,
        age       = leeftijd_jaar,
        bmi       = bmi,
        sex       = sex,
        ckd       = as.integer(egfr < 60), 
        egfr      = egfr,
        cr_resid  = cr_resid,     
        icu_admit = icu_admit,
        icu_disc  = icu_disc,
        status
      ) %>%
      mutate(
        start     = 0,
        stop_raw  = as.numeric(difftime(icu_disc, icu_admit, units = "days"))
      )

    days_to_censor <- as.numeric(difftime(censor_date,
      dat$icu_admit,
      units = "days"
    ))
    days_to_censor[days_to_censor < 0] <- 0 

    late_evt <- dat$stop_raw > days_to_censor
    dat$status[late_evt] <- 0 

    dat <- dat %>%
      mutate(
        stop = pmin(stop_raw, los_cap, days_to_censor)
      ) %>%
      filter(stop >= 0) 
    if (!keep_duplicates) {
      dat <- dplyr::distinct(dat, patient, .keep_all = TRUE)
    }
    dat
  }) %>%
    set_names(paste0("imp_", seq_len(imps$m)))
}
