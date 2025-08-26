test_that("predict_inla_cif returns finite vector on toy data", {
  skip_if_not(exists("bayes_list", .GlobalEnv),
              "Run pipeline first to create bayes_list and cp_list")
  
  toy <- cp_list[[1]][1:2, ]              # two patients
  risk <- predict_inla_cif(bayes_list[[1]], toy, 30)$death
  expect_length(risk, 2)
  expect_true(all(is.finite(risk)))
})
