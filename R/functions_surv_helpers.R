suppressPackageStartupMessages({
  library(cmprsk)
  library(tibble)
  library(dplyr)
})

cuminc_df <- function(time, status, cause = 1L, cencode = 0L) {
  stopifnot(length(time) == length(status))
  ok <- is.finite(time) & !is.na(status)
  if (sum(ok) < 2L) return(tibble::tibble(time = numeric(0), cif = numeric(0)))
  
  fstatus <- as.integer(as.character(status[ok]))
  
  ci <- try(cmprsk::cuminc(ftime = time[ok], fstatus = fstatus, cencode = cencode),
            silent = TRUE)
  if (inherits(ci, "try-error")) return(tibble::tibble(time = numeric(0), cif = numeric(0)))
  
  want <- as.character(cause)
  nm   <- names(ci)
  nm_t <- trimws(nm)
  
  last_tok <- vapply(strsplit(nm_t, "\\s+"), function(z) utils::tail(z, 1), character(1))
  idx <- which(last_tok == want)
  
  if (!length(idx)) idx <- which(nm_t == want)
  if (!length(idx)) idx <- which(nm == paste0(" ", want))
  
  if (!length(idx)) return(tibble::tibble(time = numeric(0), cif = numeric(0)))
  
  comp <- ci[[idx[1]]]
  if (is.null(comp) || length(comp$time) == 0L)
    return(tibble::tibble(time = numeric(0), cif = numeric(0)))
  
  out <- tibble::tibble(
    time = c(0, as.numeric(comp$time)),
    cif  = c(0, as.numeric(comp$est))
  )
  out <- out[order(out$time), , drop = FALSE]
  out[!duplicated(out$time), , drop = FALSE]
}
