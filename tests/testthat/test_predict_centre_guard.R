testthat::test_that("predict_inla_cif stops on unseen study_grp", {
  skip_if_no_models()
  
  newdata <- toy_cp
  newdata$study_grp <- factor("ZZ")     
  testthat::expect_error(
    predict_inla_cif(bayes_list[[1]], newdata),
    "new study_grp"
  )
})
