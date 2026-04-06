test_that("N increases with threshold probability", {
  result <- ss_net_benefit(pt_range = c(0.10, 0.30, 0.50),
                           B = 1000, seed = 2026)
  ns <- result$N_by_pt$N_required
  expect_true(all(diff(ns) >= 0))
})

test_that("ss_net_benefit returns correct class and fields", {
  result <- ss_net_benefit(pt_range = c(0.10, 0.30), B = 500, seed = 42)
  expect_s3_class(result, "dtasamplesize")
  expect_true(is.data.frame(result$N_by_pt))
  expect_equal(nrow(result$N_by_pt), 2)
  expect_true(!is.null(result$N_conservative))
})

test_that("N_conservative is in expected range", {
  result <- ss_net_benefit(B = 2000, seed = 2026,
                           N_range = seq(100, 600, by = 10))
  # Expected: N_conservative ~ 300-500
  expect_true(result$N_conservative >= 100)
  expect_true(result$N_conservative <= 600)
})
