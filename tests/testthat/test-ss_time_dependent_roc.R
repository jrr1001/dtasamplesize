
test_that("ss_time_dependent_roc runs without error", {
  skip_if_not_installed("timeROC")
  expect_no_warning(
    result <- ss_time_dependent_roc(B = 50, seed = 42,
      N_range = seq(100, 600, by = 100),
      censoring_rates = 0.20)
  )
  expect_s3_class(result, "dtasamplesize")
  # Converged strictly inside the grid, not pinned to its ceiling.
  expect_lt(result$results$N_required, 600)
  expect_gte(result$results$prob_achieved, 0.80)
})

test_that("N increases with censoring rate", {
  skip_if_not_installed("timeROC")
  result <- ss_time_dependent_roc(B = 100, seed = 2026,
    N_range = seq(100, 600, by = 50),
    censoring_rates = c(0.10, 0.30))
  ns <- result$results$N_required
  expect_true(ns[2] >= ns[1])
})

test_that("A short N_range still warns and does not silently converge", {
  skip_if_not_installed("timeROC")
  expect_warning(
    ss_time_dependent_roc(B = 30, seed = 42, N_range = c(100, 150),
                          censoring_rates = 0.20),
    "target precision not"
  )
})
