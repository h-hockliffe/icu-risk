testthat::test_that("extract_model_metrics handles cox list(death, discharge)", {
  skip_if_not(exists("cox_list_cluster", .GlobalEnv))
  mm <- extract_model_metrics(bayes_list_cluster, cox_list_cluster, fg_list_cluster)
  expect_true(is.data.frame(mm))
  expect_true(any(grepl("Cox frailty", mm$model)))
})