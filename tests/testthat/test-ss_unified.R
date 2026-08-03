
test_that("Unified N exceeds Buderer", {
  expect_no_warning(
    result <- ss_unified(B = 1000, seed = 2026,
                         N_range = seq(200, 1200, by = 50),
                         delta_auc = 0, check_nb = FALSE)
  )
  buderer_N <- ceiling(buderer_n(0.85, 0.07) / (5 / 20))
  expect_true(result$N_effective > buderer_N)
})

test_that("Unified runs with AUC enabled", {
  expect_no_warning(
    result <- ss_unified(B = 500, seed = 2026,
                         N_range = seq(300, 1200, by = 50),
                         delta_auc = 0.06, check_nb = FALSE)
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(result$joint_assurance >= 0.80)
})

test_that("Unified runs with NB enabled and informative Sp prior", {
  # The CI-based NB criterion is genuinely demanding: N_range must reach ~1600.
  expect_no_warning(
    result <- ss_unified(B = 500, seed = 2026,
                         prior_sp = c(18, 2),
                         N_range = seq(400, 2200, by = 100),
                         delta_auc = 0, check_nb = TRUE)
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(result$joint_assurance >= 0.80)
})

test_that("M-1: the CI-based NB criterion demands MORE N than no NB check", {
  # Under the old point-estimate criterion, check_nb = TRUE barely moved N,
  # because "NB_hat > 0 and NB_hat > NB_all" is satisfied with near-certainty
  # whenever the test is useful at all. With the CI-based criterion it bites.
  common <- list(B = 400, seed = 2026, prior_sp = c(18, 2),
                 delta_auc = 0, N_range = seq(400, 2400, by = 100))
  without <- do.call(ss_unified, c(common, list(check_nb = FALSE)))
  with_nb <- do.call(ss_unified, c(common, list(check_nb = TRUE)))
  expect_gt(with_nb$N_effective, without$N_effective)
})

test_that("M-2: the assurance denominator is B (degenerate draws are failures)", {
  # A tiny N with a low-prevalence prior produces many degenerate replications
  # (n_d < 5). Those used to be dropped from the denominator, yielding an
  # assurance CONDITIONAL on the study not being degenerate. Now they are
  # failures, so the assurance at a small N must be very low, not undefined
  # or inflated.
  res <- suppressWarnings(
    ss_unified(B = 500, seed = 2026, prior_prev = c(1, 40),
               N_range = c(30, 40), delta_auc = 0, check_nb = FALSE)
  )
  expect_true(res$joint_assurance >= 0)
  expect_true(res$joint_assurance < 0.10)
  expect_false(is.na(res$joint_assurance))
})

test_that("Unified returns comparison table", {
  expect_no_warning(
    result <- ss_unified(B = 500, seed = 42,
                         N_range = seq(200, 1200, by = 50),
                         delta_auc = 0, check_nb = FALSE)
  )
  expect_true(is.data.frame(result$comparison))
  expect_equal(nrow(result$comparison), 3)
  expect_match(result$comparison$method[2], "Rogan-Gladen")
})
