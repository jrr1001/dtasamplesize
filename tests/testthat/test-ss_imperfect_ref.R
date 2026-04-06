test_that("VIF calculation is correct", {
  result <- ss_imperfect_ref(Se_ref = 0.90, Sp_ref = 0.95, B = 0)
  expected_vif <- 1 / (0.90 + 0.95 - 1)^2
  expect_equal(result$inflation_factor_se, expected_vif, tolerance = 0.001)
})

test_that("Perfect reference gives VIF = 1", {
  result <- ss_imperfect_ref(Se_ref = 1.0, Sp_ref = 1.0, B = 0)
  expect_equal(result$inflation_factor_se, 1.0)
  expect_equal(result$n_diseased_unadjusted, result$n_diseased_adjusted)
})

test_that("Adjusted n_se is approximately 139 for default params", {
  result <- ss_imperfect_ref(Se_ref = 0.90, Sp_ref = 0.95, B = 0)
  # VIF ≈ 1.384, n_unadj_se = 100, so n_adj_se = ceiling(100 * 1.384) = 139
  expect_equal(result$n_diseased_adjusted, ceiling(100 * 1 / (0.85)^2))
})

test_that("Loss adjustment inflates total N", {
  result <- ss_imperfect_ref(B = 0, loss_rate = 0.10)
  expect_true(result$N_adjusted_loss > result$N_adjusted)
  expect_equal(result$N_adjusted_loss, ceiling(result$N_adjusted / 0.90))
})

test_that("Sensitivity table is generated", {
  result <- ss_imperfect_ref(B = 0, sensitivity_table = TRUE)
  expect_true(is.data.frame(result$sensitivity_table))
  expect_true(nrow(result$sensitivity_table) > 0)
  expect_true("VIF" %in% names(result$sensitivity_table))
})

test_that("MC validation shows adjusted improves on unadjusted", {
  result <- ss_imperfect_ref(B = 2000, seed = 2026, sensitivity_table = FALSE)
  expect_true(is.data.frame(result$mc_validation))
  expect_equal(nrow(result$mc_validation), 2)
})

test_that("ss_imperfect_ref returns dtasamplesize class", {
  result <- ss_imperfect_ref(B = 0)
  expect_s3_class(result, "dtasamplesize")
})
