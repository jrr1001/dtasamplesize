test_that("Unified N exceeds Buderer", {
  result <- ss_unified(B = 1000, seed = 2026,
    N_range = seq(200, 1000, by = 50),
    delta_auc = 0, check_nb = FALSE)
  buderer_N <- ceiling(buderer_n(0.85, 0.07) / (5 / 20))
  expect_true(result$N_effective > buderer_N)
})

test_that("Unified runs with AUC enabled", {
  result <- ss_unified(B = 500, seed = 2026,
    N_range = seq(300, 1000, by = 50),
    delta_auc = 0.06, check_nb = FALSE)
  expect_s3_class(result, "dtasamplesize")
  expect_true(result$joint_assurance >= 0.80)
})

test_that("Unified runs with NB enabled and informative Sp prior", {
  result <- ss_unified(B = 500, seed = 2026,
    prior_sp = c(18, 2),
    N_range = seq(300, 1000, by = 50),
    delta_auc = 0, check_nb = TRUE)
  expect_s3_class(result, "dtasamplesize")
})

test_that("Unified returns comparison table", {
  result <- ss_unified(B = 500, seed = 42,
    N_range = seq(200, 1000, by = 50),
    delta_auc = 0, check_nb = FALSE)
  expect_true(is.data.frame(result$comparison))
  expect_equal(nrow(result$comparison), 3)
})
