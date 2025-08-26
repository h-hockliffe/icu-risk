suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(gtsummary)
  library(cmprsk)
  library(purrr)
  library(tidyr)
  library(splines)
  library(riskRegression)
  library(scales)
  library(forcats)
  library(patchwork)
  library(gridExtra)
})

DEFAULT_STATUS_MAP <- c(death = 1L, discharge = 2L)

pct_stat <- function(var) {
  if (inherits(var, "factor") && nlevels(var) > 2) "{n} 
  ({style_percent(p.row, accuracy = 0.1)})"
  else                                             
    "{n} ({style_percent(p,     accuracy = 0.1)})"
}

beta_or_zero <- function(vec, name) {
  if (name %in% names(vec)) vec[[name]] else 0
}

empty_plot <- function(msg = "Prediction unavailable") {
  ggplot() +
    annotate("text",
             x = 0.5, y = 0.5,
             label = msg, size = 5, vjust = 0.5, hjust = 0.5
    ) +
    theme_void()
}

ckd_to_int <- function(x) {
  if (is.numeric(x)) {
    as.integer(x)
  } else {
    as.integer(as.character(x))
  }
}


.round_df <- function(df, digits = 3) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) round(x, digits))
  df
}

.table_grob <- function(df, n_max = 30, digits = 3, base_size = 8) {
  df <- .round_df(df, digits)
  if (nrow(df) > n_max) {
    keep <- unique(round(seq(1, nrow(df), length.out = n_max)))
    df <- df[keep, , drop = FALSE]
  }
  gridExtra::tableGrob(
    df, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      base_size = base_size,
      core = list(fg_params = list(hjust = 1, x = 0.98))
    )
  )
}

.attach_table_below <- function(p, df, n_max = 30, digits = 3, 
                                heights = c(3, 1)) {
  tbl <- .table_grob(df, n_max = n_max, digits = digits)
  patchwork::wrap_plots(p, patchwork::wrap_elements(full = tbl,
                                                    clip = "off"), ncol = 1, heights = heights)
}


make_descriptive_tables <- function(df,
                                    by_var = "icu_death",
                                    lump_large_cats = TRUE,
                                    max_levels = 20,
                                    include_longtext = TRUE) {
  
  stopifnot(by_var %in% names(df))
  
  if (!include_longtext) {
    df <- df %>% dplyr::select(-tidyselect::matches("^opnamereden_tekst$"))
  }
  
  if (lump_large_cats) {
    df <- df %>% dplyr::mutate(
      dplyr::across(where(~ is.factor(.x) || is.character(.x)),
                    ~ {
                      nlev <- dplyr::n_distinct(.x)
                      if (nlev > max_levels)
                        forcats::fct_lump_n(as.factor(.x), n = max_levels,
                                            other_level = "Other")
                      else .x
                    }
      )
    )
  }
  
  vals <- sort(unique(na.omit(df[[by_var]])))
  if (!all(vals %in% c(0, 1))) {
    warning(sprintf("`%s` is not 0/1;
    attempting to coerce via min->0, max->1.", by_var))
    m <- min(vals); M <- max(vals)
    df[[by_var]] <- as.integer(df[[by_var]] == M)
  }
  
  df <- df %>% dplyr::mutate(
    ..outcome = factor(
      !!rlang::sym(by_var),
      levels = c(0, 1),
      labels = c("Discharged alive", "ICU death")
    )
  )
  
  opname_col <- names(df)[grepl("^opnamereden_tekst$", names(df))]
  test_list <- list(
    all_continuous()  ~ "wilcox.test",
    all_categorical() ~ "chisq.test"
  )
  if (length(opname_col) == 1) test_list[[opname_col]] <- "chisq.test"
  
  tbl1 <- tbl_summary(
    df,
    by        = ..outcome,
    type      = all_continuous() ~ "continuous2",
    percent   = "column",
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})"),
    missing   = "no"
  ) |>
    add_p(test = test_list,
          test.args = list(all_continuous() ~ list(simulate.p.value =
                                                     TRUE))) |>
    bold_labels() |>
    modify_header(label ~ "**Variable**") |>
    modify_spanning_header(all_stat_cols() ~ "**ICU outcome**")
  
  tbl2 <- tbl_summary(
    df,
    type      = all_continuous() ~ "continuous2",
    statistic = list(all_continuous() ~ "{median} ({p25}–{p75})"),
    missing   = "no"
  ) |>
    bold_labels() |>
    modify_header(label ~ "**Total cohort (n = {N})**")
  
  list(baseline_by_outcome = tbl1, overall_cohort = tbl2)
}

plot_cumulative_incidence <- function(
    inla_fit, cp_data,
    horizon = NULL,
    status_map = DEFAULT_STATUS_MAP,
    xmax_full = 120, dt = 0.5,
    annotate_zoom_at = 7,
    annotate_full_at = 30,
    show_cif_labels = TRUE,        
    show_markers_full = TRUE,      
    append_table = FALSE,          
    table_times = c(0, 1, 3, 7, 15, 30),  
    table_max = 30, table_digits = 3      
) {
  to_list <- function(x) if (!is.list(x) || is.data.frame(x)) list(x)
  else x
  inla_list <- to_list(inla_fit)
  cp_list   <- to_list(cp_data)
  stopifnot(length(cp_list) >= 1L, length(inla_list) >= 1L)
  
  t_max  <- max(dplyr::distinct(cp_list[[1]], patient, 
                                .keep_all = TRUE)$stop, na.rm = TRUE)
  t_grid <- sort(unique(c(0, seq(dt, t_max, by = dt))))
  if (!is.null(horizon)) t_grid <- t_grid[t_grid <= max(horizon)]
  if (!length(t_grid)) stop("plot_cumulative_incidence(): empty time
  grid.")
  
  m_imp <- min(length(cp_list), length(inla_list))
  
  get_cif_step <- function(dat, cause_code) {
    df <- cuminc_df(dat$stop, dat$status, cause = cause_code, cencode = 0L)
    if (!nrow(df)) return(function(tt) rep(0, length(tt)))
    sf <- stats::stepfun(df$time, c(0, df$cif), right = FALSE)
    function(tt) sf(tt)
  }
  
  obs_mat <- function(cause_code) {
    M <- matrix(NA_real_, nrow = length(t_grid), ncol = m_imp)
    for (i in seq_len(m_imp)) {
      dat  <- dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE)
      step <- get_cif_step(dat, cause_code)
      M[, i] <- step(t_grid)
    }
    M
  }
  
  orient <- function(M, tlen) {
    if (nrow(M) == tlen) return(M)
    if (ncol(M) == tlen) return(t(M))
    stop("CIF matrix has incompatible dimensions: ", paste(dim(M),
                                                           collapse = "×"),
         " vs grid length ", tlen)
  }
  obs_mat_death     <- orient(obs_mat(status_map[["death"]]), 
                              length(t_grid))
  obs_mat_discharge <- orient(obs_mat(status_map[["discharge"]]),
                              length(t_grid))
  
  rubin_ci <- function(mat) {
    qbar <- rowMeans(mat, na.rm = TRUE)
    B    <- apply(mat, 1, stats::var, na.rm = TRUE); B[is.na(B)] <- 0
    se   <- sqrt(pmax(B * (1 + 1 / m_imp), 0))
    tibble::tibble(mean = qbar, lwr = pmax(0, qbar - 1.96 * se),
                   upr = pmin(1, qbar + 1.96 * se))
  }
  
  obs_death_ci     <- rubin_ci(obs_mat_death)
  obs_discharge_ci <- rubin_ci(obs_mat_discharge)
  
  obs_df <- dplyr::bind_rows(
    tibble::tibble(time = t_grid, cif = obs_death_ci$mean,
                   lwr  = obs_death_ci$lwr, upr = obs_death_ci$upr,
                   cause = "Death",   source = "Observed"),
    tibble::tibble(time = t_grid, cif = obs_discharge_ci$mean,
                   lwr  = obs_discharge_ci$lwr, upr = obs_discharge_ci$upr,
                   cause = "Discharge", source = "Observed")
  )
  
  counts <- table(dplyr::distinct(cp_list[[1]], patient,
                                  .keep_all = TRUE)$status)
  if (isTRUE(counts[as.character(status_map[["discharge"]])] == 0))
    warning("No discharge events found — observed discharge CIF will be 0.")
  if (isTRUE(counts[as.character(status_map[["death"]])] == 0))
    warning("No death events found — observed death CIF will be 0.")
  
  pred_one_imp <- function(i) {
    dat <- dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE)
    sapply(t_grid, function(h) {
      pr <- predict_inla_cif(inla_list[[i]], dat, h)
      c(mean_death = mean(pr$death,     na.rm = TRUE),
        mean_disc  = mean(pr$discharge, na.rm = TRUE))
    })
  }
  
  pred_arr <- array(NA_real_, dim = c(2L, length(t_grid), m_imp))
  for (i in seq_len(m_imp)) pred_arr[, , i] <- pred_one_imp(i)
  
  pred_death_ci <- rubin_ci(pred_arr[1, , , drop = TRUE])
  pred_disc_ci  <- rubin_ci(pred_arr[2, , , drop = TRUE])
  
  pred_df <- dplyr::bind_rows(
    tibble::tibble(time = t_grid, cif = pred_death_ci$mean,
                   cause = "Death",     source = "Predicted"),
    tibble::tibble(time = t_grid, cif = pred_disc_ci$mean,
                   cause = "Discharge", source = "Predicted")
  )
  
  cif_tab <- dplyr::bind_rows(
    dplyr::select(obs_df, time, cause, source, cif, lwr, upr),
    dplyr::mutate(pred_df, lwr = NA_real_, upr = NA_real_)
  ) |>
    dplyr::arrange(time, cause, source)
  
  
  if (length(table_times)) {
    cif_tab <- dplyr::filter(cif_tab, time %in% 
                               intersect(table_times, unique(obs_df$time)))
  }
  
  cols <- c(Death = "#D55E00", Discharge = "#009E73")
  
  base_p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = dplyr::filter(obs_df, source == "Observed"),
      ggplot2::aes(time, ymin = lwr, ymax = upr, fill = cause),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = pred_df,
      ggplot2::aes(time, cif, color = cause, linetype = source),
      linewidth = 0.9
    ) +
    ggplot2::geom_step(
      data = dplyr::filter(obs_df, source == "Observed"),
      ggplot2::aes(time, cif, color = cause, linetype = source),
      direction = "hv", linewidth = 1.2
    ) +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_linetype_manual(values = c(Observed = "solid", 
                                              Predicted = "22"), name = NULL) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Days since ICU admission", y = "Cumulative incidence",
      title = "Cause-specific cumulative incidence (Rubin-pooled)"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
  
  lab_y <- function(val) pmin(0.98, pmax(0.02, val + 0.03))
  
  annot_at <- function(p, x, d_cif, s_cif, xlim) {
    p <- p + ggplot2::geom_vline(xintercept = x,
                                 linetype = "dotted", linewidth = 0.3)
    if (!isTRUE(show_cif_labels)) return(p)   
    x_txt <- xlim[2] - 0.03 * diff(xlim)
    p +
      ggplot2::annotate(
        "text", x = x_txt, y = lab_y(d_cif),
        label = sprintf("Death CIF@%dd \u2248 %s", x, 
                        scales::percent(d_cif, 0.1)),
        hjust = 1, size = 3
      ) +
      ggplot2::annotate(
        "text", x = x_txt, y = lab_y(s_cif),
        label = sprintf("Discharge CIF@%dd \u2248 %s", x,
                        scales::percent(s_cif, 0.1)),
        hjust = 1, size = 3
      )
  }
  
  pick_idx_at <- function(tt, target) {
    idx <- max(which(tt <= target))
    if (!is.finite(idx)) idx <- which.min(abs(tt - target))
    idx
  }
  
  p_7 <- base_p +
    ggplot2::coord_cartesian(xlim = c(0, annotate_zoom_at), expand = FALSE)
  
  if (length(annotate_zoom_at)) {
    i_z <- pick_idx_at(t_grid, annotate_zoom_at)
    d_z <- obs_death_ci$mean[i_z]; s_z <- obs_discharge_ci$mean[i_z]
    p_7 <- annot_at(p_7, annotate_zoom_at, d_z, s_z, c(0, 
                                                       annotate_zoom_at))
  }
  
  p_full <- base_p +
    ggplot2::coord_cartesian(xlim = c(0, xmax_full), expand = FALSE)
  
  if (isTRUE(show_markers_full)) {
    p_full <- p_full + ggplot2::geom_vline(xintercept = c(1, 7, 30), 
                                           linetype = "dotted", linewidth = 0.3)
  }
  
  days_full <- as.numeric(annotate_full_at)
  if (length(days_full)) {
    for (x in days_full) {
      i <- pick_idx_at(t_grid, x)
      p_full <- annot_at(p_full, x, obs_death_ci$mean[i], 
                         obs_discharge_ci$mean[i], c(0, xmax_full))
    }
  }
  
  if (isTRUE(append_table)) {
    p_7    <- .attach_table_below(p_7,    cif_tab, n_max = table_max,
                                  digits = table_digits)
    p_full <- .attach_table_below(p_full, cif_tab, n_max = table_max, 
                                  digits = table_digits)
  }
  
  list(zoom_7d = p_7, full_120d = p_full)
}

plot_decision_curve <- function(inla_fit, cox_fit, fg_fit,
                                cp_data,
                                time_horizon = c(0.5, 1, 3, 7, 15, 30, 60,
                                                 90),
                                thresholds = NULL,
                                max_thr = NULL,
                                status_map = DEFAULT_STATUS_MAP,
                                show_treat_all = TRUE,        
                                clip_to_models = TRUE,        
                                clip_quantiles = c(0.02, 0.98),  
                                verbose = TRUE,
                                append_table = FALSE,         
                                table_every = 5,              
                                table_max = 30, table_digits = 4  
) {
  
  to_list <- function(x) if (!is.list(x) || is.data.frame(x)) list(x) else
    x
  stopifnot(length(time_horizon) >= 1)
  
  cp_list   <- to_list(cp_data)
  inla_list <- to_list(inla_fit)
  cox_list  <- to_list(cox_fit)
  fg_list   <- to_list(fg_fit)
  
  if (!length(cp_list)) return(empty_plot("No data"))
  if (!(length(inla_list) == length(cox_list) &&
        length(inla_list) == length(fg_list)  &&
        length(inla_list) == length(cp_list))) {
    warning("plot_decision_curve(): list lengths differ – using common 
    indices.")
  }
  idxs <- seq_len(min(length(cp_list), length(inla_list), 
                      length(cox_list), length(fg_list)))
  
  d0 <- dplyr::distinct(cp_list[[1]], patient, .keep_all = TRUE)
  crude_inc <- mean(d0$status == status_map[["death"]], na.rm = TRUE)
  
  if (is.null(max_thr)) {
    max_thr <- pmin(pmax(0.6, 4 * crude_inc), 0.8)
  }
  if (is.null(thresholds)) thresholds <- seq(0.005, max_thr, by = 0.005)
  thresholds <- thresholds[is.finite(thresholds) & thresholds > 0 &
                             thresholds < 1]
  if (!length(thresholds)) return(empty_plot("No thresholds"))
  
  if (!(is.numeric(clip_quantiles) &&
        length(clip_quantiles) == 2L &&
        all(is.finite(clip_quantiles)) &&
        clip_quantiles[1] > 0 && clip_quantiles[2] < 1 &&
        clip_quantiles[1] < clip_quantiles[2])) {
    clip_quantiles <- c(0.02, 0.98)
  }
  
  nb_one <- function(dat, inla_mod, cox_mod, fg_mod, h) {
    nd <- dplyr::distinct(dat, patient, .keep_all = TRUE)
    
    r_inla <- suppressWarnings(try(predict_inla_cif(inla_mod, nd,
                                                    h)$death, silent = TRUE))
    r_cox  <- suppressWarnings(try(predict_cox_cif(cox_mod,  nd, h),
                                   silent = TRUE))
    r_fg   <- suppressWarnings(try(predict_fg_cif(fg_mod,   nd, h),
                                   silent = TRUE))
    
    coerce_len <- function(x, n) {
      v <- if (inherits(x, "try-error") || is.null(x)) rep(NA_real_, n)
      else as.numeric(x)
      if (length(v) != n) v <- rep(NA_real_, n)
      pmin(pmax(v, 0), 1)
    }
    n <- nrow(nd)
    r_inla <- coerce_len(r_inla, n)
    r_cox  <- coerce_len(r_cox,  n)
    r_fg   <- coerce_len(r_fg,   n)
    
    event_h <- as.integer((nd$status == status_map[["death"]])
                          & (nd$stop <= h))
    inc_h   <- mean(event_h, na.rm = TRUE)
    
    net_benefit <- function(risk, th) {
      ok <- is.finite(risk) & is.finite(event_h)
      if (!any(ok)) return(NA_real_)
      tp <- mean((risk[ok] >= th) & (event_h[ok] == 1))
      fp <- mean((risk[ok] >= th) & (event_h[ok] == 0))
      tp - fp * th / (1 - th)
    }
    
    nb <- purrr::map_dfr(thresholds, function(th) {
      tibble::tibble(th = th,
                     INLA = net_benefit(r_inla, th),
                     Cox  = net_benefit(r_cox,  th),
                     FG   = net_benefit(r_fg,   th))
    })
    
    list(
      nb = nb,
      treat_all = tibble::tibble(
        th = thresholds,
        model = "Treat all",
        nb = inc_h - (1 - inc_h) * thresholds / (1 - thresholds)
      ),
      inc_h = inc_h,
      n = n
    )
  }
  
  per_h <- lapply(time_horizon, function(h) {
    tmp <- lapply(idxs, function(i)
      nb_one(dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE),
             inla_list[[i]], cox_list[[i]], fg_list[[i]], h))
    
    nb_raw <- dplyr::bind_rows(lapply(tmp, `[[`, "nb"))
    nb_df  <- nb_raw |>
      tidyr::pivot_longer(c(INLA, Cox, FG), 
                          names_to = "model", values_to = "nb") |>
      dplyr::group_by(th, model) |>
      dplyr::summarise(nb = mean(nb, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(horizon = h)
    
    ta_df <- dplyr::bind_rows(lapply(tmp, `[[`, "treat_all")) |>
      dplyr::group_by(th, model) |>
      dplyr::summarise(nb = mean(nb, na.rm = TRUE), .groups = "drop") |>
      dplyr::mutate(horizon = h)
    
    inc_h <- mean(vapply(tmp, `[[`, numeric(1), "inc_h"))
    n_h   <- mean(vapply(tmp, `[[`, numeric(1), "n"))
    
    list(h = h, nb = dplyr::bind_rows(nb_df, ta_df),
         inc_h = inc_h, n = n_h)
  })
  
  nb_df <- dplyr::bind_rows(lapply(per_h, `[[`, "nb"))
  tn_df <- tibble::tibble(
    horizon = rep(time_horizon, each = length(thresholds)),
    th      = rep(thresholds,   times = length(time_horizon)),
    model   = "Treat none",
    nb      = 0
  )
  nb_df <- dplyr::bind_rows(nb_df, tn_df)
  
  vline_df <- tibble::tibble(
    horizon = vapply(per_h, `[[`, numeric(1), "h"),
    x       = vapply(per_h, `[[`, numeric(1), "inc_h") * 100
  )
  
  diag_tb <- tibble::tibble(
    horizon        = vapply(per_h, `[[`, numeric(1), "h"),
    crude_inc      = crude_inc,
    inc_h          = vapply(per_h, `[[`, numeric(1), "inc_h"),
    n              = vapply(per_h, `[[`, numeric(1), "n"),
    thresholds_min = min(thresholds),
    thresholds_max = max(thresholds),
    max_thr        = max_thr,
    n_imputations  = length(idxs)
  )
  
  if (isTRUE(verbose)) {
    msg <- sprintf("DCA: %d imputations | %d horizons | unique patients
    (first imp.): %d | crude incidence: %.3f",
                   length(idxs), length(time_horizon), nrow(d0), crude_inc)
    message(msg)
  }
  
  if (!isTRUE(show_treat_all)) {
    nb_df <- nb_df %>% dplyr::filter(model != "Treat all")
  }
  
  max_x <- max(thresholds) * 100
  base_title <- if (length(time_horizon) == 1L)
    sprintf("Decision-curve analysis (%s-day)", time_horizon)
  else
    "Decision-curve analysis (mean over imputations)"
  
  p <- ggplot2::ggplot(nb_df, ggplot2::aes(th * 100, nb, colour = model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(
      data = vline_df,
      mapping = ggplot2::aes(xintercept = x),
      linetype = "dashed"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, max_x, by = if (max_x <= 30) 2 else 5),
      labels = function(x) paste0(x, "%"),
      limits = c(min(thresholds) * 100, max_x)
    ) +
    ggplot2::scale_colour_manual(values = c(
      INLA = "#0072B2", Cox = "#D55E00", FG = "#009E73",
      `Treat all` = "grey40", `Treat none` = "grey20"
    )) +
    ggplot2::labs(x = "Risk threshold (%)",
                  y = "Net benefit",
                  title = base_title, colour = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position  = "right")
  
  if (length(unique(nb_df$horizon)) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ~ horizon,
      labeller = ggplot2::labeller(horizon = function(x) paste0(x, "-day"))
    )
  }
  
  if (isTRUE(clip_to_models)) {
    keep_models <- c("INLA", "Cox", "FG", "Treat none")
    rng_df <- nb_df %>%
      dplyr::filter(model %in% keep_models, is.finite(nb))
    
    if (nrow(rng_df) > 0) {
      lo <- stats::quantile(rng_df$nb, clip_quantiles[1], na.rm = TRUE)
      hi <- stats::quantile(rng_df$nb, clip_quantiles[2], na.rm = TRUE)
      if (is.finite(lo) && is.finite(hi) && hi > lo) {
        pad <- max(1e-3, (hi - lo) * 0.15)
        p <- p + ggplot2::coord_cartesian(ylim = c(lo - pad, hi + pad))
      }
    }
  }
  
  if (isTRUE(append_table)) {
    nb_tab <- nb_df
    nb_tab <- nb_tab[order(nb_tab$horizon, nb_tab$th, nb_tab$model),
                     , drop = FALSE]
    nb_tab <- dplyr::mutate(
      nb_tab,
      horizon = as.character(.data$horizon),
      threshold = scales::percent(th, accuracy = 0.1)
    ) |>
      dplyr::select(horizon, threshold, model, nb)
    
    if (is.numeric(table_every) && table_every > 1) {
      nb_tab <- nb_tab |>
        dplyr::group_by(horizon, model) |>
        dplyr::slice(seq(1, dplyr::n(), by = table_every)) |>
        dplyr::ungroup()
    }
    p <- .attach_table_below(p, nb_tab, n_max = table_max,
                             digits = table_digits)
  }
  
  structure(p, diag_tb = diag_tb, nb_df = nb_df, vline_df = vline_df)
}

plot_calibration_models <- function(
    inla_fit, cox_fit, fg_fit, cp_data,
    horizon = c(0.5, 1, 3, 7, 15, 30, 60, 90),
    status_map = DEFAULT_STATUS_MAP,
    add_smooth = FALSE, smooth_method = "loess", smooth_span = 0.75,
    append_table = FALSE, table_max = 30, table_digits = 3
) {
  horizon <- sort(unique(as.numeric(horizon)))
  
  to_list <- function(x) {
    if (is.null(x)) return(list())
    if (is.list(x) && !is.data.frame(x)) x else list(x)
  }
  
  inla_list <- to_list(inla_fit)
  cox_list  <- to_list(cox_fit)
  fg_list   <- to_list(fg_fit)
  cp_list   <- to_list(cp_data)
  
  if (length(cp_list) == 0L)   stop("cp_list is empty.",   call. = FALSE)
  if (length(inla_list) == 0L) stop("INLA list is empty.", call. = FALSE)
  if (length(cox_list) == 0L)  stop("Cox list is empty.",  call. = FALSE)
  if (length(fg_list) == 0L)   stop("FG list is empty.",   call. = FALSE)
  n_imp <- min(length(cp_list), length(inla_list),
               length(cox_list), length(fg_list))
  
  pick_col <- function(df, candidates) {
    nm <- intersect(candidates, names(df)); if (!length(nm))
      stop("Missing needed columns."); nm[[1]]
  }
  as_num_vec <- function(x, n) {
    if (inherits(x, "try-error") || is.null(x)) return(rep(NA_real_, n))
    if (is.list(x) && !is.null(x$death)) x <- x$death
    x <- as.numeric(x); if (length(x) != n) rep(NA_real_, n) else x
  }
  
  pooled <- purrr::map_dfr(seq_len(n_imp), function(i) {
    dat <- dplyr::distinct(cp_list[[i]], patient, .keep_all = TRUE)
    time_col   <- pick_col(dat, c("time", "stop", "tstop", "ftime"))
    status_col <- pick_col(dat, c("status", "status_cr", "event"))
    purrr::map_dfr(horizon, function(h) {
      n <- nrow(dat)
      tibble::tibble(
        patient = dat$patient,
        time    = dat[[time_col]],
        status  = dat[[status_col]],
        horizon = h,
        INLA    = as_num_vec(try(predict_inla_cif(inla_list[[i]],
                                                  dat, h), silent = TRUE), n),
        Cox     = as_num_vec(try(predict_cox_cif(cox_list[[i]], 
                                                 dat, h), silent = TRUE), n),
        FG      = as_num_vec(try(predict_fg_cif(fg_list[[i]], 
                                                dat, h), silent = TRUE), n)
      )
    })
  })
  
  pooled <- pooled |>
    tidyr::pivot_longer(c(INLA, Cox, FG), names_to = "model",
                        values_to = "risk") |>
    dplyr::group_by(patient, horizon, model) |>
    dplyr::summarise(
      risk   = mean(risk, na.rm = TRUE),
      status = dplyr::first(status),
      time   = dplyr::first(time),
      .groups = "drop"
    ) |>
    dplyr::mutate(event_h = as.integer((status ==
                                          status_map[["death"]]) & (time <= horizon)))
  
  make_one <- function(h) {
    cal_df <- pooled |>
      dplyr::filter(horizon == h, is.finite(risk)) |>
      dplyr::group_by(model) |>
      dplyr::filter(dplyr::n() > 0) |>
      dplyr::mutate(bin = dplyr::ntile(risk, 10)) |>
      dplyr::group_by(model, bin) |>
      dplyr::summarise(
        meanRisk = mean(risk,    na.rm = TRUE),
        obsRisk  = mean(event_h, na.rm = TRUE),
        n        = dplyr::n(),
        .groups  = "drop"
      ) |>
      dplyr::filter(n > 0, is.finite(meanRisk), is.finite(obsRisk))
    
    lim_top <- min(0.5, ceiling(
      max(stats::quantile(cal_df$meanRisk, 0.99, na.rm = TRUE),
          stats::quantile(cal_df$obsRisk,  0.99, na.rm = TRUE),
          na.rm = TRUE) * 20
    ) / 20)
    
    p <- ggplot2::ggplot(cal_df, ggplot2::aes(meanRisk, obsRisk,
                                              colour = model, group = model)) +
      ggplot2::geom_point(alpha = 0.9, size = 2) +
      ggplot2::geom_line(alpha = 0.9, linewidth = 0.7) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::scale_x_continuous("Predicted risk",
                                  labels = scales::label_percent(accuracy = 1)) +
      ggplot2::scale_y_continuous("Observed risk",
                                  labels = scales::label_percent(accuracy = 1)) +
      ggplot2::coord_cartesian(xlim = c(0, lim_top),
                               ylim = c(0, lim_top), expand = TRUE) +
      ggplot2::labs(title = paste0("Calibration (", h, "-day)")) +
      ggplot2::theme_bw() +
      ggplot2::theme(panel.grid.minor = element_blank(),
                     legend.position = "bottom",
                     legend.title = element_blank())
    if (isTRUE(add_smooth)) {
      p <- p + ggplot2::geom_smooth(method = smooth_method, se = FALSE,
                                    linewidth = 0.6, alpha = 0.6, span =
                                      smooth_span)
    }
    if (isTRUE(append_table)) {
      cal_tab <- cal_df |>
        dplyr::arrange(bin, model) |>
        dplyr::transmute(
          bin, model, n,
          meanRisk = meanRisk, obsRisk = obsRisk
        )
      p <- .attach_table_below(p, cal_tab, n_max = table_max, digits = table_digits)
    }
    p
  }
  
  plots <- setNames(lapply(horizon, make_one), paste0(horizon, "d"))
  if (length(plots) == 1L) plots[[1]] else plots
}

plot_cindex_td <- function(cindex_td, horizons = TIME_HORIZONS,
                           show_ci = TRUE) {
  stopifnot(is.list(cindex_td), !is.null(cindex_td$pooled))
  
  df <- cindex_td$pooled
  df <- dplyr::filter(df, is.finite(C_td))
  
  if (!nrow(df)) {
    return(empty_plot("No C-index values available"))
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Horizon,
                                        y = C_td, color = model)) +
    ggplot2::geom_hline(yintercept = 0.5,
                        linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::geom_line(linewidth = 0.8, alpha = 0.9) +
    ggplot2::scale_x_continuous(breaks = horizons) +
    ggplot2::scale_y_continuous(limits = c(0.5, 1), 
                                breaks = seq(0.5, 1, 0.05)) +
    ggplot2::labs(
      title = "Time-dependent C-index by horizon",
      x = "Days since ICU admission",
      y = "C-index",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom")
  
  if (isTRUE(show_ci) && !is.null(cindex_td$raw)) {
    ci_df <- cindex_td$raw |>
      dplyr::group_by(Horizon, model) |>
      dplyr::summarise(
        mean = mean(C_td, na.rm = TRUE),
        se   = stats::sd(C_td,  na.rm = TRUE) / sqrt(sum(is.finite(C_td))),
        .groups = "drop"
      ) |>
      dplyr::mutate(lwr = pmax(0.5, mean - 1.96 * se),
                    upr = pmin(1.0, mean + 1.96 * se))
    
    if (nrow(ci_df)) {
      p <- p +
        ggplot2::geom_ribbon(
          data = ci_df,
          ggplot2::aes(x = Horizon, ymin = lwr, ymax = upr,
                       fill = model, group = model),
          inherit.aes = FALSE, alpha = 0.12, linewidth = 0
        ) +
        ggplot2::guides(fill = "none")
    }
  }
  
  p
}
