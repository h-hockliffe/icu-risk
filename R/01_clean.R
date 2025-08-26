suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(forcats)
  library(janitor)
  library(arrow)
  library(purrr)
  library(readr)
  requireNamespace("progressr", quietly = TRUE)
})


RAW_FILE <- Sys.getenv(
  "RAW_FILE",
  unset = "G:/Drive partagés/OI - DSI - Thesis Hockliffe/FINAL BAYESIAN_07042025 (1).xlsx"
)

CENSOR <- as.POSIXct("2023-12-31 23:59:59", tz = "UTC")
MAX_REAL_DATE <- CENSOR + lubridate::days(365)
OUT_PARQ <- "data/clean.parquet"

safe_reader <- purrr::insistently(
  readxl::read_excel,
  rate = purrr::rate_backoff(
    pause_base = 2,
    pause_cap = 60,
    jitter = TRUE,
    max_times = 3
  )
)

read_with_progress <- function(path) {
  if (!requireNamespace("progressr", quietly = TRUE)) {
    return(safe_reader(path, guess_max = 1e5, progress = TRUE))
  }
  progressr::with_progress({
    p <- progressr::progressor(steps = 1)
    dat <- safe_reader(path, guess_max = 1e5, progress = TRUE)
    p()
    dat
  })
}

# ------------------------------------------------------------------ #
# 2.  Master cleaner
# ------------------------------------------------------------------ #
#' clean_master()
#' Read, clean and persist the analysis dataset.
#'
#' @param raw  Path to the raw Excel workbook
#' @param out  Output parquet file 
#' @return     Cleaned tibble 
clean_master <- function(
    raw = RAW_FILE,
    out = OUT_PARQ,
    deduplicate = TRUE
) {
  raw_norm <- tryCatch(
    normalizePath(raw,
                  winslash = "/",
                  mustWork = FALSE
    ),
    error = function(e) raw
  )
  if (!file.exists(raw)) {
    stop("Raw workbook not found at:\n   ", raw_norm,
         "\nIs the shared drive mounted?",
         call. = FALSE
    )
  }
  
  message("Reading raw workbook …  (up to 3 retries on failure)")
  df <- read_with_progress(raw_norm)
  
  if (!is.data.frame(df)) {
    stop("Failed to read raw workbook after 3 attempts.", call. = FALSE)
  }
  
  df <- df %>%
    janitor::clean_names() %>%
    rename(
      icu_death = tidyselect::matches("^overlijden_op_icu_0_nee_1_ja_?$"),
      mort1m    = tidyselect::matches("^x1_month_mortality_0_nee_1_ja$"),
      mort3m    = tidyselect::matches("^x3_month_mortality_0_nee_1_ja$"),
      study_id  = tidyselect::matches("^herkomst$"),
      sex       = tidyselect::matches("^geslacht_1_m_2_v$"),
      icu_admit = tidyselect::matches("^opnamedatum_ite$"),
      icu_disc  = tidyselect::matches("^ontslagdatum_icu"),
      hosp_disc = tidyselect::matches("^ontslagdatum_ziekenhuis"),
      bmi_raw   = tidyselect::matches("^bmi_kg_m2$")
    )
  
  df <- df %>%
    mutate(
      icu_death = as.integer(icu_death == 1L),
      mort1m    = as.integer(mort1m == 1L),
      mort3m    = as.integer(mort3m == 1L)
    )
  
  df <- df %>%
    filter(
      leeftijd_jaar >= 18,
      icu_admit <= icu_disc,
      icu_disc  <= hosp_disc,
      hosp_disc <= MAX_REAL_DATE,
      year(icu_admit) <= 2023
    )
  
  to_numeric_checked <- function(vec, varname) {
    vec_chr <- as.character(vec)
    num_vec <- suppressWarnings(as.numeric(vec_chr))
    probs <- which(is.na(num_vec) & !is.na(vec_chr) & nzchar(vec_chr))
    if (length(probs)) {
      warning(sprintf("%s: %d value(s) could not be coerced to numeric and were set to NA.", varname, length(probs)),
              call. = FALSE, immediate. = TRUE
      )
    }
    num_vec
  }
  
  df <- df %>%
    mutate(
      creat_raw = to_numeric_checked(creatinine_mg_dl, "creatinine_mg_dl"),
      creat_conv = case_when(
        is.na(creat_raw) ~ NA_real_,
        creat_raw > 30 & creat_raw < 3000 ~ creat_raw / 88.4,
        TRUE ~ creat_raw
      ),
      creat = case_when(
        is.na(creat_conv) ~ NA_real_,
        creat_conv > 17 ~ NA_real_,
        creat_conv < 0.23 ~ NA_real_,
        TRUE ~ creat_conv
      ),
      egfr_raw = to_numeric_checked(e_gfr, "e_gfr"),
      egfr = case_when(
        is.na(egfr_raw) ~ NA_real_,
        egfr_raw < 1 | egfr_raw > 250 ~ NA_real_,
        TRUE ~ egfr_raw
      ),
      egfr_capped = as.integer(!is.na(egfr) & egfr > 120),
      bmi = to_numeric_checked(bmi_raw, "bmi_kg_m2")
    )
  
  n_conv <- sum(!is.na(df$creat_raw) &
                  df$creat_raw > 30 & df$creat_raw < 3000 &
                  !is.na(df$creat), na.rm = TRUE)
  n_high <- sum(df$creat_raw > 17, na.rm = TRUE)
  
  message(sprintf(
    "✓ Converted %d creatinine values from µmol/L → mg/dl.",
    n_conv
  ))
  if (n_high > 0) {
    warning(
      sprintf(
        "%d creatinine values > 17 mg/dl coerced to NA.",
        n_high
      ),
      call. = FALSE, immediate. = TRUE
    )
  }
  
  df <- df %>%
    mutate(
      sex      = factor(sex, levels = c(1, 2), labels = c("Man", "Vrouw")),
      study_id = factor(study_id)
    )
  
  df <- collapse_study(
    df,
    event_var = "icu_death",
    min_alive = CFG$collapse_min_alive,
    min_dead  = CFG$collapse_min_dead,
    new_level = "Merged-small-sites"
  ) %>%
    dplyr::rename(study_grp = study_id)
  
  df <- df %>%
    mutate(
      los_icu    = as.numeric(difftime(icu_disc, icu_admit, units = "days")),
      los_after  = as.numeric(difftime(hosp_disc, icu_disc, units = "days")),
      los_total  = as.numeric(difftime(hosp_disc, icu_admit, units = "days"))
    ) %>%
    mutate(across(
      c(los_icu, los_after, los_total),
      ~ if_else(.x < 0, NA_real_, .x)
    ))
  
  dup_ids <- df %>%
    count(studie_patient_id, name = "row_n") %>%
    filter(row_n > 1L)
  
  if (deduplicate) {
    if (nrow(dup_ids)) {
      warning(sprintf(
        "clean_master(): %d Study IDs have >1 row; keeping earliest ICU_admit.",
        nrow(dup_ids)
      ), call. = FALSE, immediate. = TRUE)
      df <- df %>%
        arrange(studie_patient_id, icu_admit) %>%
        distinct(studie_patient_id, .keep_all = TRUE)
    } else {
      message("✓ No duplicate Study IDs detected.")
    }
  } else {
    if (nrow(dup_ids)) {
      message(sprintf(
        "✓ Retaining duplicate Study IDs for patient‑level clustering (%d duplicates).",
        nrow(dup_ids)
      ))
    }
  }
  

  clean_index <- df %>%
    group_by(studie_patient_id) %>%    
    slice_max(icu_disc, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  check_index_terminal <- function(dat) {
    stopifnot(
      !anyDuplicated(dat$studie_patient_id),      
      all(dat$los_icu >= 0),
      all(dat$icu_death %in% c(0, 1)),
      all(dat$icu_disc >= dat$icu_admit)
    )
    icu_death_n <- sum(dat$icu_death == 1, na.rm = TRUE)
    icu_alive_n <- sum(dat$icu_death == 0, na.rm = TRUE)
    if ((icu_death_n + icu_alive_n) == 0)
      stop("No terminal events found for index stays.")
    invisible(dat)
  }
  clean_index <- check_index_terminal(clean_index)
  

  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(clean_index, out)
  if (exists("coercion_summary", mode = "function")) {
    cs <- try(coercion_summary(df), silent = TRUE)  
    if (!inherits(cs, "try-error") && is.data.frame(cs)) {
      dir.create("data", recursive = TRUE, showWarnings = FALSE)
      readr::write_csv(cs, file = file.path("data", "coercion_summary.csv"))
      message("Wrote coercion summary: data/coercion_summary.csv")
    }
  }
  message(
    "Wrote cleaned data: ", out,
    "  (n = ", nrow(clean_index), ",  p = ", ncol(clean_index), ")"
  )
  
  invisible(clean_index)
}

