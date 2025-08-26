test_that("overall_epv_check stops when EPV < threshold", {
  df <- data.frame(icu_death = c(rep(1L, 5), rep(0L, 95)))
  expect_error(overall_epv_check(df, event_var = "icu_death", n_par = 10, min_epp = 10),
               "Overall EPV")
})

test_that("overall_epv_check passes when EPV >= threshold", {
  df <- data.frame(icu_death = c(rep(1L, 50), rep(0L, 50)))
  expect_silent(overall_epv_check(df, event_var = "icu_death", n_par = 5, min_epp = 10))
})

test_that("epp_guardrail warns (not stop) when some imputations < 5 EPV", {
  cp1 <- make_cp_toy(60, seed = 1)
  cp2 <- make_cp_toy(60, seed = 2)
  cp2$status[cp2$status == 1L] <- 0L  
  expect_warning(
    epp_guardrail(list(imp_1 = cp1, imp_2 = cp2), min_epp = 5, n_par = 12, fail_hard = FALSE),
    "Imputation-level EPV below"
  )
})
