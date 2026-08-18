
test_that("ss_net_benefit returns correct class and fields", {
  result <- ss_net_benefit(pt_range = c(0.10, 0.30), B = 1000, seed = 42)
  expect_s3_class(result, "dtasamplesize")
  expect_true(is.data.frame(result$N_by_pt))
  expect_equal(nrow(result$N_by_pt), 2)
  expect_true(all(c("pt", "feasible", "N_required", "prob_achieved",
                    "NB_true", "NB_treat_all") %in% names(result$N_by_pt)))
  expect_true(!is.null(result$N_conservative))
  expect_identical(result$design, "cohort")   # new default
})

test_that("design is cohort by default and 'fixed' is still available", {
  coh <- ss_net_benefit(pt_range = 0.30, B = 1000, seed = 1)
  fix <- ss_net_benefit(pt_range = 0.30, B = 1000, seed = 1, design = "fixed")
  expect_identical(coh$design, "cohort")
  expect_identical(fix$design, "fixed")
  expect_error(ss_net_benefit(design = "bogus", B = 10), "should be one of")
})

test_that("the cohort design requires MORE subjects than the fixed design", {
  # This is the whole point of C-2. The fixed-margin variance conditions away
  # the sampling variability of the prevalence, so it understates the true
  # variance of a prospective cohort and therefore understates N.
  coh <- ss_net_benefit(pt_range = 0.30, B = 4000, seed = 2026)
  fix <- ss_net_benefit(pt_range = 0.30, B = 4000, seed = 2026,
                        design = "fixed")
  expect_gt(coh$N_by_pt$N_required[1], fix$N_by_pt$N_required[1])
})

test_that("Useful thresholds are feasible and return a sample size", {
  expect_no_warning(result <- ss_net_benefit(B = 2000, seed = 2026))
  # With Se=0.85, Sp=0.90, prev=0.20 the test is useful at every threshold
  expect_true(all(result$N_by_pt$feasible))
  expect_true(all(!is.na(result$N_by_pt$N_required)))
  expect_true(is.finite(result$N_conservative))
})

test_that("Infeasible thresholds are flagged, not silently sized", {
  # Se = Sp = 0.55, prev = 0.10: at pt = 0.05 treat-all dominates the test
  result <- ss_net_benefit(Se = 0.55, Sp = 0.55, prev = 0.10,
                           pt_range = c(0.05), B = 1000, seed = 1)
  expect_false(result$N_by_pt$feasible[1])
  expect_true(is.na(result$N_by_pt$N_required[1]))
})

test_that("Criterion is non-trivial: extreme thresholds need more N", {
  # CI-based criterion -> sample size is U-shaped in pt (largest where the
  # test is hardest to distinguish from a default strategy), not constant.
  result <- ss_net_benefit(pt_range = c(0.10, 0.25, 0.50),
                           B = 3000, seed = 2026)
  ns <- result$N_by_pt$N_required
  names(ns) <- as.character(result$N_by_pt$pt)
  # Guard: all three thresholds must be feasible and met, otherwise the
  # shape comparisons below would compare against NA.
  expect_false(anyNA(ns))
  expect_true(ns["0.5"] >= ns["0.25"])
  expect_true(ns["0.1"] >= ns["0.25"])
  expect_true(ns["0.5"] > min(ns))  # not all pinned to the search floor
})

test_that("target_assurance is an equivalent alias of target_prob", {
  by_prob <- ss_net_benefit(pt_range = 0.30, target_prob = 0.75,
                            B = 1000, seed = 1)
  by_assurance <- ss_net_benefit(pt_range = 0.30, target_assurance = 0.75,
                                 B = 1000, seed = 1)
  expect_identical(by_prob$N_by_pt, by_assurance$N_by_pt)
})

test_that("target_assurance takes precedence over target_prob when both are supplied", {
  result <- ss_net_benefit(pt_range = 0.30, target_prob = 0.99,
                           target_assurance = 0.75, B = 1000, seed = 1)
  expected <- ss_net_benefit(pt_range = 0.30, target_prob = 0.75,
                             B = 1000, seed = 1)
  expect_identical(result$N_by_pt, expected$N_by_pt)
})

test_that("Required N grows as assurance target increases", {
  lo <- ss_net_benefit(pt_range = 0.50, target_prob = 0.70,
                       B = 3000, seed = 2026)
  hi <- ss_net_benefit(pt_range = 0.50, target_prob = 0.90,
                       B = 3000, seed = 2026)
  expect_true(hi$N_by_pt$N_required[1] >= lo$N_by_pt$N_required[1])
})

# --- Independent re-simulation of the cohort DGP ---------------------------
# Helper: measure the TRUE assurance at a given N by simulating the cohort
# data-generating process from scratch (different seed, criterion re-derived
# here rather than reused from the package).
cohort_assurance <- function(N, Se, Sp, prev, pt, alpha = 0.05,
                             M = 40000, seed = 99) {
  w <- pt / (1 - pt)
  z <- stats::qnorm(1 - alpha / 2)
  set.seed(seed)
  n_d  <- stats::rbinom(M, N, prev)
  n_nd <- N - n_d
  TP <- stats::rbinom(M, n_d, Se)
  FP <- stats::rbinom(M, n_nd, 1 - Sp)
  FN <- n_d - TP
  TN <- n_nd - FP

  p1 <- TP / N
  p2 <- FP / N
  nb  <- p1 - w * p2
  vnb <- (p1 * (1 - p1) + w^2 * p2 * (1 - p2) + 2 * w * p1 * p2) / N

  q1 <- FN / N
  q2 <- TN / N
  dd  <- w * q2 - q1
  vd  <- (q1 * (1 - q1) + w^2 * q2 * (1 - q2) + 2 * w * q1 * q2) / N

  mean((nb - z * sqrt(pmax(vnb, 0)) > 0) & (dd - z * sqrt(pmax(vd, 0)) > 0))
}

test_that("the DECLARED cohort assurance matches an INDEPENDENT simulation", {
  res <- ss_net_benefit(pt_range = 0.30, target_prob = 0.80,
                        B = 4000, seed = 2026)
  N <- res$N_by_pt$N_required[1]
  expect_false(is.na(N))

  real <- cohort_assurance(N, Se = 0.85, Sp = 0.90, prev = 0.20, pt = 0.30,
                           seed = 99)
  # The recommended N must genuinely deliver the target, and the reported
  # assurance must agree with the independently measured one.
  expect_gt(real, 0.78)
  expect_lt(abs(real - res$N_by_pt$prob_achieved[1]), 0.03)
})

test_that("the fixed-margin design OVERSTATES the assurance of a cohort", {
  # Take the N the legacy fixed design recommends, then measure what a real
  # prospective cohort would actually achieve at that N.
  fix <- ss_net_benefit(pt_range = 0.30, target_prob = 0.80,
                        B = 4000, seed = 2026, design = "fixed")
  N <- fix$N_by_pt$N_required[1]

  real <- cohort_assurance(N, Se = 0.85, Sp = 0.90, prev = 0.20, pt = 0.30,
                           seed = 123)
  # The fixed design claims >= 0.80 at this N; a real cohort falls short.
  expect_gte(fix$N_by_pt$prob_achieved[1], 0.80)
  expect_lt(real, 0.80)
})
