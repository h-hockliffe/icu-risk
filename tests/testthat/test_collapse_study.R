testthat::test_that("collapse_study errors on level clash", {
  df <- data.frame(
    study_id = factor("Merged-small-sites"),
    icu_death = 0
  )
  testthat::expect_error(
    collapse_study(df, new_level = "Merged-small-sites"),
    "must not clash"
  )
})

