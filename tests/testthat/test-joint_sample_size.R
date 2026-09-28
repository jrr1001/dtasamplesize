
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
  # Degenerate iff n_d = 0 or n_nd = 0 (n_d = N); an arm of size 1 is NOT
  # degenerate (see ?joint_sample_size, @details).
  p_degenerate <- stats::dbinom(0, N, prev) + stats::dbinom(N, N, prev)
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

# --- H-06: unified degeneracy rule (n = 0 is degenerate; n = 1 is not) -----
# These pin the SAME threshold used by bam_sample_size(): an arm is
# degenerate, and always scored as a failure, only when it received ZERO
# subjects (n_d = 0 or n_nd = 0), the point at which wilson_width() is
# literally undefined (division by n). An arm of size ONE is NOT degenerate:
# the Wilson interval is well-defined there (see ?joint_sample_size,
# @details) and is scored on its actual width.

test_that("design = \"fixed\": an arm of size 1 is scored on its Wilson width, not auto-failed", {
  # With n_d_exp = 1, wilson_width(x, 1) is EXACTLY 0.7934507 regardless of
  # whether x = 0 or x = 1 (the interval is symmetric at n = 1), so the Se
  # arm's pass/fail is deterministic across every replicate, independent of
  # Se and of the random draw. This isolates the n = 1 rule with no Monte
  # Carlo noise on that arm. N = 4, prev = 0.25 gives n_d_exp = 1, n_nd_exp
  # = 3 exactly (floor(4 * 0.25) = 1).
  # delta_auc = 0.5 (target AUC width 1.0) is needed to clear the
  # deterministic AUC gate at these small expected margins (n_d_exp = 1,
  # n_nd_exp = 3): hanley_mcneil_var(0.90, 1, 3) gives an AUC width of
  # 0.9475, verified separately.
  target_width_n1 <- 2 * 0.4 # 0.80 > 0.7934507 -> the n = 1 arm always PASSES
  res_pass <- suppressWarnings(joint_sample_size(
    Se = 0.85, Sp = 0.90, prev = 0.25,
    delta_se = 0.40, delta_sp = 0.45, delta_auc = 0.5,
    N_range = c(4), B = 1000, seed = 11, design = "fixed"
  ))
  expect_equal(res_pass$n_diseased, 1)
  expect_equal(res_pass$joint_prob_se_sp, 1)

  target_width_fail <- 2 * 0.39 # 0.78 < 0.7934507 -> the n = 1 arm always FAILS
  res_fail <- suppressWarnings(joint_sample_size(
    Se = 0.85, Sp = 0.90, prev = 0.25,
    delta_se = 0.39, delta_sp = 0.45, delta_auc = 0.5,
    N_range = c(4), B = 1000, seed = 11, design = "fixed"
  ))
  expect_equal(res_fail$n_diseased, 1)
  expect_equal(res_fail$joint_prob_se_sp, 0)
})

test_that("design = \"cohort\": only n_d = 0 or n_nd = 0 replicates are forced failures", {
  # Targets loose enough (full width 1.0, i.e. the theoretical maximum
  # Wilson width) that EVERY non-degenerate replicate (n_d in 1..N-1)
  # passes on both arms; only n_d = 0 or n_d = N (n_nd = 0) can fail. The
  # joint probability should then match the analytic non-degeneracy
  # probability under n_d ~ Binomial(N, prev), not the old (n < 2) rule's
  # analytic bound.
  N <- 10; prev <- 0.30; B <- 40000
  res <- suppressWarnings(joint_sample_size(
    Se = 0.85, Sp = 0.90, prev = prev,
    delta_se = 0.5, delta_sp = 0.5, delta_auc = 0.5,
    N_range = c(N), B = B, seed = 7, design = "cohort"
  ))
  p_non_degenerate <- 1 - stats::dbinom(0, N, prev) - stats::dbinom(N, N, prev)
  expect_false(is.na(res$joint_prob_se_sp))
  expect_equal(res$joint_prob_se_sp, p_non_degenerate, tolerance = 4 * sqrt(0.25 / B))
})

test_that("an expected margin of exactly 0 is skipped as a candidate N (n = 0, not n = 1)", {
  # N = 2 with prev = 0.30 gives n_d_exp = floor(0.6) = 0: this candidate N
  # must be skipped (as before), so the AUC gate is never even evaluated
  # for it and auc_gate_passed stays FALSE -- hence the "AUC precision
  # target" warning below (the same code path an unmet AUC gate takes),
  # not the (unrelated) "No N in N_range achieved" warning. N = 4 with the
  # same prev gives n_d_exp = 1, which must NOT be skipped now that the
  # threshold is n = 0 (checked below).
  expect_warning(
    res_skipped <- joint_sample_size(prev = 0.30, delta_se = 0.001, delta_sp = 0.001,
                       N_range = c(2), B = 100, seed = 1),
    "AUC precision target"
  )
  expect_equal(res_skipped$n_diseased, 0) # confirms N = 2 was never evaluated
  expect_no_warning(
    res <- joint_sample_size(Se = 0.85, Sp = 0.90, prev = 0.30,
                              delta_se = 0.45, delta_sp = 0.45, delta_auc = 0.5,
                              N_range = c(4), B = 100, seed = 1, design = "fixed")
  )
  expect_equal(res$n_diseased, 1)
})
