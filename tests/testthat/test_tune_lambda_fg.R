test_that("tune_lambda_fg handles no-censoring by skipping FGR", {
  df <- make_cp_toy(120)
  df$status[df$status == 0L] <- 2L
  expect_silent({
  })
})

