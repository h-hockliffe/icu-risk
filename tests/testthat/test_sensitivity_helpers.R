test_that("drop_merged_small_sites removes that stratum", {
  df <- make_cp_toy(50)
  kept <- drop_merged_small_sites(df, level = "Merged-small-sites")
  expect_false("Merged-small-sites" %in% levels(kept$study_grp))
})

test_that("add_creatinine_residual adds cr_resid orthogonal to log-eGFR", {
  df <- make_cp_toy(80)
  df$creat_umol <- pmax(30, 1000 / pmax(df$egfr, 1) + rnorm(nrow(df), 0, 20))
  out <- add_creatinine_residual(df, creat_umol_col = "creat_umol", egfr_col = "egfr")
  expect_true("cr_resid" %in% names(out))
  expect_lt(abs(cor(out$cr_resid, log(pmax(out$egfr, 1)))), 0.1)
})

test_that("limit_to_ckd_negative filters to egfr >= 60", {
  df <- make_cp_toy(60)
  out <- limit_to_ckd_negative(df, egfr_col = "egfr", threshold = 60)
  expect_true(all(out$egfr >= 60, na.rm = TRUE))
})
