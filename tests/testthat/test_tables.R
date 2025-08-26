test_that("baseline table builds without error", {
  clean <- arrow::read_parquet("data/clean.parquet")
  expect_error(make_descriptive_tables(clean), NA)
})
