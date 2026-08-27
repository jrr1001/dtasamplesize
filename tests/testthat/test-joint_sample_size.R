
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

# --- design = "cohort" (default) vs "fixed" (legacy, versions <= 0.4.0) --

test_that("design defaults to \"cohort\" and is echoed in the result", {
  expect_equal(formals(joint_sample_size)[["design"]],
               quote(c("cohort", "fixed")))
  result <- joint_sample_size(B = 500, seed = 42,
                              N_range = seq(100, 800, by = 20))
  expect_equal(result$design, "cohort")
})

test_that("an invalid design is rejected by match.arg", {
  expect_error(
    joint_sample_size(design = "bogus", B = 100),
    "'arg' should be one of"
  )
})

test_that("design = \"fixed\" reproduces the pre-0.5.0 fixed-margin behaviour", {
  # Adversarial operating point verified independently at B = 300000:
  # fixed-margin joint probability at N = 650 is 0.8071. B is kept moderate
  # here to run fast in the suite; the tolerance covers the resulting MC
  # error (MCSE ~= sqrt(0.81*0.19/8000) ~= 0.0044, so 6*MCSE ~= 0.026).
  # N_range is a single point below the search's own target_prob, so the
  # "target not reached, consider expanding N_range" warning is expected
  # and not the object of this test.
  result <- suppressWarnings(
    joint_sample_size(Se = 0.70, Sp = 0.80, prev = 0.20,
                      delta_se = 0.08, delta_sp = 0.06,
                      N_range = c(650), B = 8000, seed = 2026,
                      design = "fixed")
  )
  expect_equal(result$design, "fixed")
  expect_equal(result$joint_prob_se_sp, 0.8071, tolerance = 0.03)
})

test_that("design = \"cohort\" gives a LOWER joint probability than \"fixed\"", {
  # Conditioning on the expected number of diseased subjects (design =
  # "fixed") ignores the sampling variability of that count in a
  # prospective cohort, and so overstates the assurance. At this
  # adversarial operating point the gap is large (~0.09) and B = 8000 is
  # far more than enough to resolve it (MCSE ~= 0.004-0.005 per arm).
  common <- list(Se = 0.70, Sp = 0.80, prev = 0.20, delta_se = 0.08,
                 delta_sp = 0.06, N_range = c(650), B = 8000, seed = 2026)
  fixed <- suppressWarnings(
    do.call(joint_sample_size, c(common, list(design = "fixed")))
  )
  cohort <- suppressWarnings(
    do.call(joint_sample_size, c(common, list(design = "cohort")))
  )
  expect_lt(cohort$joint_prob_se_sp, fixed$joint_prob_se_sp)
})

test_that("design = \"cohort\" is reproducible given the same seed", {
  common <- list(Se = 0.70, Sp = 0.80, prev = 0.20, delta_se = 0.08,
                 delta_sp = 0.06, N_range = c(650), B = 2000, seed = 2026,
                 design = "cohort")
  r1 <- suppressWarnings(do.call(joint_sample_size, common))
  r2 <- suppressWarnings(do.call(joint_sample_size, common))
  expect_identical(r1$joint_prob_se_sp, r2$joint_prob_se_sp)
})

test_that("degenerate cohort replicates count as failures (denominator B)", {
  # With loose CI-width targets almost every NON-degenerate replicate
  # passes, so the joint probability is bounded above by the population
  # probability that a replicate is NOT degenerate under n_d ~ Bin(N, prev).
  # If degenerate replicates were instead dropped from the denominator (the
  # bug already fixed in ss_unified for versions <= 0.2.0), the joint
  # probability could sit arbitrarily close to 1 regardless of this bound.
  Se <- 0.85; Sp <- 0.90; prev <- 0.10; N <- 30
  B <- 5000
  res <- suppressWarnings(
    joint_sample_size(Se = Se, Sp = Sp, prev = prev,
                      delta_se = 0.30, delta_sp = 0.30, delta_auc = 0.30,
                      N_range = c(N), B = B, seed = 1, design = "cohort")
  )
  p_degenerate <- stats::pbinom(1, N, prev) +
    (1 - stats::pbinom(N - 2, N, prev))
  upper_bound <- 1 - p_degenerate
  expect_false(is.na(res$joint_prob_se_sp))
  # small MC slack around the population bound (worst-case MCSE at p = 0.5)
  expect_lte(res$joint_prob_se_sp, upper_bound + 4 * sqrt(0.25 / B))
})

test_that("the manuscript scenario under the cohort default reaches ~0.80 at N = 580", {
  # Independently verified at B = 300000: joint_prob_se_sp at N = 580 is
  # 0.8080 under the cohort default (vs. 0.8558 previously reported under
  # the fixed-margin design). This locks in that the cohort default no
  # longer overstates the assurance for the published scenario.
  result <- suppressWarnings(joint_sample_size(B = 2000, seed = 2026,
                                                prev = 0.20,
                                                N_range = seq(100, 900, by = 20)))
  expect_equal(result$design, "cohort")
  expect_equal(result$n_total, 580)
  expect_gte(result$joint_prob_se_sp, 0.80)
})
