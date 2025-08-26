testthat::test_that("core pipeline artefacts exist", {
  skip_if_no_models()
  
  testthat::expect_true(inherits(bayes_list, "list"))
  testthat::expect_gt(nrow(model_metrics), 0)
})
