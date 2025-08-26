test_that("choose_intervals enforces min_events per cause", {
  set.seed(1)
  n <- 400
  df <- data.frame(
    stop = rexp(n, 0.05),
    status = sample(c(0L,1L,2L), n, replace = TRUE, prob = c(0.4, 0.3, 0.3))
  )
  cu <- choose_intervals(df, n_max = 5, min_events = 10, max_iter = 50, verbose = FALSE)
  labs <- cut(df$stop, breaks = cu, include.lowest = TRUE, right = TRUE)
  tab  <- df[df$status %in% c(1,2), ] |>
    dplyr::mutate(interval = labs) |>
    dplyr::count(interval, status) |>
    tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)
  expect_true(all(tab$`1` >= 10))
  expect_true(all(tab$`2` >= 10))
})



