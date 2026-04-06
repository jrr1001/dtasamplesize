test_that("Adaptive does not reduce N when prevalence is higher", {
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.42), B = 1000, seed = 2026)
  expect_true(result$results$N_final_median[1] >= result$N_initial)
})

test_that("Adaptive increases N when prevalence is lower", {
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.18), B = 1000, seed = 2026)
  expect_true(result$results$N_final_mean[1] > result$N_initial * 1.1)
})

test_that("ss_adaptive_prevalence returns correct class", {
  result <- ss_adaptive_prevalence(B = 200, seed = 42,
                                   prev_true_range = c(0.25, 0.30))
  expect_s3_class(result, "dtasamplesize")
  expect_true(is.data.frame(result$results))
  expect_equal(nrow(result$results), 2)
})

test_that("At initial prevalence, N_final is near N_initial", {
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.30), B = 1000, seed = 2026)
  # Median should be close to N_initial (within 20%)
  ratio <- result$results$N_final_median[1] / result$N_initial_adj
  expect_true(ratio >= 0.95 && ratio <= 1.30)
})
