testthat::test_that("validate_models pooling keeps original names", {
  fake <- tibble::tibble(
    model = c("INLA","Cox","FG"),
    c_index = c(0.75, 0.72, 0.70),
    intercept = c(0.0, 0.0, 0.0),
    slope = c(1.0, 1.0, 1.0),
    c_index_se = c(0.02, 0.02, 0.03),
    intercept_se = c(NA, NA, NA),
    slope_se = c(0.1, 0.1, 0.1)
  )
  
  pooled <- fake |>
    dplyr::group_by(model) |>
    dplyr::summarise(
      dplyr::across(
        c(c_index, intercept, slope, c_index_se, intercept_se, slope_se),
        ~ list(.), .names = "{.col}"
      ),
      .groups = "drop"
    )
  
  expect_true(all(c("c_index","c_index_se","intercept","intercept_se","slope","slope_se") %in% names(pooled)))
})