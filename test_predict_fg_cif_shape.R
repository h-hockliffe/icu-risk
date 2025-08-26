test_that("predict_fg_cif returns numeric vector (crr engine)", {
  skip_if_no_pkg("cmprsk")
  
  df <- make_cp_toy(180)
  df1 <- df[!duplicated(df$patient), ]
  X <- model.matrix(~ log(egfr) + log(bmi) + sex, data = transform(df1, sex = factor(sex)))[, -1]
  fit <- try(
    cmprsk::crr(ftime = df1$stop, fstatus = df1$status, cov1 = X, failcode = 1, cencode = 0),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) skip("crr fit failed on toy data")
  
  res <- predict_fg_cif(fit, df1, horizon = 30)
  expect_type(res, "double")
  expect_length(res, nrow(df1))
  expect_gte(min(res, na.rm = TRUE), 0)
  expect_lte(max(res, na.rm = TRUE), 1)
})
