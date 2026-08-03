
test_that("Joint N exceeds Buderer max", {
  expect_no_warning(
    result <- joint_sample_size(B = 1000, seed = 2026,
                                N_range = seq(100, 800, by = 20))
  )
  buderer_N <- max(buderer_n(0.85, 0.07) / 0.20,
                   buderer_n(0.90, 0.05) / 0.80)
  expect_true(result$n_total > ceiling(buderer_N))
})

test_that("joint_sample_size returns correct class and fields", {
  expect_no_warning(
    result <- joint_sample_size(B = 500, seed = 42,
                                N_range = seq(100, 800, by = 20))
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(!is.null(result$n_total))
  expect_true(!is.null(result$joint_prob_se_sp))
  expect_true(!is.null(result$buderer_N))
  expect_true(result$joint_prob_se_sp >= 0 && result$joint_prob_se_sp <= 1)
  expect_true(result$auc_gate_passed)
})

test_that("joint_prob is a deprecated alias of joint_prob_se_sp", {
  result <- joint_sample_size(B = 500, seed = 42,
                              N_range = seq(100, 800, by = 20))
  expect_identical(result$joint_prob, result$joint_prob_se_sp)
})

test_that("Joint N is in expected range", {
  expect_no_warning(
    result <- joint_sample_size(B = 3000, seed = 2026,
                                N_range = seq(100, 800, by = 10))
  )
  expect_true(result$n_total >= 300)
  expect_true(result$n_total <= 800)
})

# --- (a): geometric coherence of AUC with Se and Sp ---------------------

test_that("An AUC below the concave-ROC minimum is refused", {
  # Se = 0.85, Sp = 0.90 => AUC_min = 0.5*(1-Sp)*Se + 0.5*Sp*(Se+1) = 0.875.
  # The old default of AUC = 0.80 was therefore geometrically impossible, and
  # the package happily returned a sample size for it.
  expect_error(
    joint_sample_size(Se = 0.85, Sp = 0.90, AUC = 0.80, B = 100),
    "geometrically impossible"
  )
  expect_error(
    joint_sample_size(Se = 0.85, Sp = 0.90, AUC = 0.80, B = 100),
    "0\\.8750"
  )
})

test_that("AUC_min matches the closed form and the default AUC clears it", {
  result <- joint_sample_size(B = 200, seed = 1,
                              N_range = seq(100, 800, by = 50))
  Se <- 0.85
  Sp <- 0.90
  expect_equal(result$AUC_min, 0.5 * (1 - Sp) * Se + 0.5 * Sp * (Se + 1))
  expect_equal(result$AUC_min, 0.875)
  expect_gte(result$AUC, result$AUC_min)   # default AUC = 0.90
})

test_that("An AUC exactly at the geometric minimum is accepted", {
  expect_error(
    suppressWarnings(
      joint_sample_size(Se = 0.85, Sp = 0.90, AUC = 0.875, B = 100,
                        N_range = seq(100, 800, by = 50))
    ),
    NA
  )
})

# --- (c): the AUC gate never passing must give NA, not 0 ----------------

test_that("joint_prob_se_sp is NA (not 0) when the AUC gate never passes", {
  # delta_auc so tight that no N in range clears it: the Se/Sp Monte Carlo is
  # never run, so its probability was never computed. Reporting 0 (the old
  # initialisation value) would be indistinguishable from a genuine 0.
  expect_warning(
    result <- joint_sample_size(delta_auc = 0.001, B = 200, seed = 1,
                                N_range = seq(100, 300, by = 50)),
    "AUC precision target"
  )
  expect_false(result$auc_gate_passed)
  expect_true(is.na(result$joint_prob_se_sp))
  expect_true(is.na(result$joint_prob))   # alias agrees
})
