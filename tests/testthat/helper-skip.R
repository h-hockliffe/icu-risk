skip_if_no_models <- function() {
  if (!exists("bayes_list", .GlobalEnv))
    testthat::skip("Pipeline artefacts not built yet")
}

skip_if_no_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    skip(paste0("Package '", pkg, "' not available"))
  }
}

skip_if_not_mac_or_linux <- function() {
  if (.Platform$OS.type != "unix") skip("Non-Unix OS (Windows) — skipping")
}
