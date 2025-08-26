suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(tidyr)
})

collapse_study <- function(df,
                           event_var = "icu_death",
                           min_alive = 30L,
                           min_dead = 5L,
                           warn_pct = 20,
                           new_level = "Merged-small-sites",
                           verbose = TRUE) {
  if (!event_var %in% names(df)) {
    stop("`", event_var, "` not found in `df`.", call. = FALSE)
  }
  if (!"study_id" %in% names(df)) {
    stop("`study_id` column missing in `df`.", call. = FALSE)
  }
  if (new_level %in% unique(df$study_id)) {
    stop("`new_level` must not clash with existing study_id values.",
      call. = FALSE
    )
  }

  df$study_id <- forcats::as_factor(df$study_id)

  if (verbose) {
    message(
      "Collapsing sparse study strata …  ",
      "[min_alive = ", min_alive,
      ", min_dead = ", min_dead, "]"
    )
  }

  cnt_tbl <- df %>%
    count(study_id, !!sym(event_var), name = "n") %>%   
    tidyr::pivot_wider(
      names_from  = !!sym(event_var),
      values_from = n, values_fill = 0, names_sort = TRUE
    ) %>%
    mutate(bad = (`0` < min_alive) | (`1` < min_dead))
  

  offenders <- cnt_tbl$study_id[cnt_tbl$bad]

  if (length(offenders) && verbose) {
    message(
      "   • Collapsing levels: ",
      paste(offenders, collapse = ", "),
      "  →  ", new_level
    )
  }

  df <- df %>%
    mutate(
      study_id = fct_other(
        study_id,
        keep         = setdiff(levels(study_id), offenders),
        other_level  = new_level
      )
    )

  final_cnt <- df %>%
    count(study_id, !!sym(event_var), name = "n") %>%   
    tidyr::pivot_wider(
      names_from  = !!sym(event_var),                   
      values_from = n, values_fill = 0, names_sort = TRUE
    ) %>%
    mutate(
      ok    = (`0` >= min_alive) & (`1` >= min_dead),
      total = `0` + `1`,
      pct   = 100 * total / sum(total)
    ) %>%
    arrange(desc(total))

  if (!all(final_cnt$ok)) {
    offenders <- final_cnt$study_id[!final_cnt$ok]
    stop(
      "After collapsing, the following study levels still violate ",
      sprintf(
        "≥ %d alive & ≥ %d deaths rule: %s\n",
        min_alive, min_dead,
        paste(offenders, collapse = ", ")
      ),
      "Consider raising the thresholds or revisiting the raw data.",
      call. = FALSE
    )
  }

  if (verbose) {
    message("Final strata after collapsing (alive / dead / % cohort):")
    print(final_cnt %>%
      select(study_id, alive = `0`, dead = `1`, total, pct) %>%
      as.data.frame())
  }

  big_pool <- final_cnt$study_id[1] == new_level &&
    final_cnt$pct[1] >= warn_pct

  if (big_pool) {
    warning(
      sprintf(
        "%s now holds %.1f %% of the cohort (≥ %d %%). ",
        new_level, final_cnt$pct[1], warn_pct
      ),
      "A sensitivity analysis **without** the centre random effect ",
      "is strongly recommended.",
      call. = FALSE, immediate. = TRUE
    )
  }

  df <- df %>%
    mutate(
      study_id = study_id %>%
        fct_drop() %>% 
        fct_infreq() 
    )


  attr(df, "collapse_summary") <- final_cnt
  attr(df, "needs_sensitivity_RE") <- big_pool

  invisible(df)
}
