test_that("compare_inla_mcmc maps parameters and computes diffs without INLA draws", {
  inla_death <- fake_inla_death_fit()
  mcmc_samp  <- fake_mcmc_samp(1500)
  
  tab <- compare_inla_mcmc(inla_death, mcmc_samp, use_inla_draws = FALSE)
  expect_true(all(c("par","mean_inla","mean_mcmc","diff","width_ratio") %in% names(tab)))
  expect_true(all(c("log_egfr","log_bmi","sexVrouw","age_spline1","age_spline2","age_spline3") %in% tab$par))
  expect_true(all(is.finite(tab$mean_inla)))
  expect_true(all(is.finite(tab$mean_mcmc)))
})
