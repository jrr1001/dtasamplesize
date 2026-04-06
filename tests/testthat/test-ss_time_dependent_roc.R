test_that("ss_time_dependent_roc runs without error", {
  skip_if_not_installed("timeROC")
  result <- ss_time_dependent_roc(B = 50, seed = 42,
    N_range = seq(100, 300, by = 50),
    censoring_rates = 0.20)
  expect_s3_class(result, "dtasamplesize")
})

test_that("N increases with censoring rate", {
  skip_if_not_installed("timeROC")
  result <- ss_time_dependent_roc(B = 100, seed = 2026,
    N_range = seq(100, 500, by = 50),
    censoring_rates = c(0.10, 0.30))
  ns <- result$results$N_required
  expect_true(ns[2] >= ns[1])
})
