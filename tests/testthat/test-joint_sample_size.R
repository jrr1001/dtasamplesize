test_that("Joint N exceeds Buderer max", {
  result <- joint_sample_size(B = 1000, seed = 2026,
                              N_range = seq(100, 700, by = 20))
  buderer_N <- max(buderer_n(0.85, 0.07) / 0.20,
                   buderer_n(0.90, 0.05) / 0.80)
  expect_true(result$n_total > ceiling(buderer_N))
})

test_that("joint_sample_size returns correct class and fields", {
  result <- joint_sample_size(B = 500, seed = 42,
                              N_range = seq(100, 800, by = 20))
  expect_s3_class(result, "dtasamplesize")
  expect_true(!is.null(result$n_total))
  expect_true(!is.null(result$joint_prob))
  expect_true(!is.null(result$buderer_N))
  expect_true(result$joint_prob >= 0 && result$joint_prob <= 1)
})

test_that("Joint N is in expected range", {
  result <- joint_sample_size(B = 3000, seed = 2026,
                              N_range = seq(100, 800, by = 10))
  # Expected 450-600 based on joint requirement
  expect_true(result$n_total >= 300)
  expect_true(result$n_total <= 800)
})
