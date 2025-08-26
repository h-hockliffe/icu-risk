suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
})


read_clean_data <- function(path = "data/clean.parquet") {
  if (!file.exists(path)) {
    stop("Clean parquet not found – run 01_clean.R or check the path.")
  }
  arrow::read_parquet(path) %>% dplyr::as_tibble()
}

#' @param df A data.frame/tibble with the columns above.
#' @return A 1-row tibble with counts.
coercion_summary <- function(df) {
  stopifnot(is.data.frame(df))
  
  n_bmi_to_na <- {
    x <- as.character(df$bmi_raw)
    num <- suppressWarnings(as.numeric(x))
    sum(!is.na(x) & nzchar(x) & is.na(num), na.rm = TRUE)
  }
  
  tibble::tibble(
    n_creat_converted = sum(
      !is.na(df$creat_raw) & df$creat_raw > 30 & df$creat_raw < 3000 &
        !is.na(df$creat), na.rm = TRUE
    ),
    n_creat_high_na  = sum(df$creat_raw > 17, na.rm = TRUE),
    n_egfr_out_of_rng = sum(
      !is.na(df$egfr_raw) & (df$egfr_raw < 1 | df$egfr_raw > 250),
      na.rm = TRUE
    ),
    n_bmi_to_na       = n_bmi_to_na
  )
}

#'
#' @param path Path to the CSV written by 01_clean.R
#' @return A tibble or NULL if the file is missing/unreadable
read_coercion_summary <- function(path = "data/coercion_summary.csv") {
  if (!file.exists(path)) return(NULL)
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
}
