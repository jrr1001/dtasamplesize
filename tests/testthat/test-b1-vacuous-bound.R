
# Regression tests for the B = 1 vacuous assurance_lower defect: all three
# decision rules of ss_unified() computed the lower confidence bound on the
# Monte Carlo assurance with the ordinary normal-approximation (Wald)
# formula, joint_assurance - z * sqrt(joint_assurance * (1 - joint_assurance)
# / B). That standard error is EXACTLY 0 whenever joint_assurance is 0 or 1,
# regardless of B: at B = 1, a single replicate that happened to pass every
# active target gave joint_assurance = 1 with a Wald SE of
# sqrt(1 * 0 / 1) = 0, so assurance_lower came out as exactly 1.000 --
# manufactured certainty from the smallest possible amount of evidence -- and
# every decision rule accepted the first N tried on that basis.
#
# assurance_lower is now a one-sided Wilson score bound (wilson_lower(),
# R/helpers.R), which stays strictly inside (0, 1) for any finite effective
# sample size. These tests fail against the pre-fix Wald formula: at
# phat = 1, n = 1, z = qnorm(0.95), the Wald lower bound is EXACTLY 1 (the
# assertions below require it to be well below 1).

wilson_lower <- get("wilson_lower", envir = asNamespace("dtasamplesize"))

test_that("wilson_lower() does not collapse to 1 at phat = 1, unlike the Wald bound it replaced", {
  z <- stats::qnorm(0.95)
  wald_lower_at_1_1 <- 1 - z * sqrt(1 * 0 / 1)  # the OLD formula, verbatim
  expect_equal(wald_lower_at_1_1, 1)  # sanity check: the old formula really is vacuous here

  expect_lt(wilson_lower(1, 1, z), 0.9)
  expect_gt(wilson_lower(1, 1, z), 0)
})

test_that("wilson_lower() does not collapse to 0 at phat = 0, unlike the Wald bound it replaced", {
  z <- stats::qnorm(0.95)
  expect_equal(0 - z * sqrt(0 * 1 / 1), 0)  # the OLD formula: also exactly 0 here
  # The Wilson bound at phat = 0 is 0 too (there is genuinely no evidence of
  # ANY success), so this is not "wrong" in the same way as the phat = 1
  # case -- included for completeness/symmetry of the formula, not because
  # the old value was itself a defect at this end.
  expect_equal(wilson_lower(0, 1, z), 0)
})

test_that("wilson_lower() converges to the Wald bound as n grows, so realistic-B results are unaffected", {
  z <- stats::qnorm(0.95)
  phat <- 0.80
  wald <- phat - z * sqrt(phat * (1 - phat) / 20000)
  wilson <- wilson_lower(phat, 20000, z)
  expect_equal(wilson, wald, tolerance = 1e-3)
})

test_that("ss_unified(B = 1, decision = 'lower_bound') never reports a vacuous assurance_lower = 1.000", {
  # A generous scenario (tight, near-perfect priors; low target) so a SINGLE
  # B = 1 replicate at a small N plausibly clears every active check --
  # exactly the condition that fabricated joint_assurance = 1 pre-fix.
  result <- suppressWarnings(ss_unified(
    prior_se = c(40, 2), prior_sp = c(40, 2), prior_prev = c(10, 10),
    Se_ref = 1, Sp_ref = 1, delta_se = 0.30, delta_sp = 0.30, delta_auc = 0,
    check_nb = FALSE, target_assurance = 0.5, decision = "lower_bound",
    N_range = seq(50, 400, by = 10), B = 1, seed = 2026
  ))
  # This scenario is deliberately generous enough that the search accepts
  # on the very first replicate tried (joint_assurance == 1); if it did not,
  # the test below would not exercise the vacuous-bound condition at all.
  expect_equal(result$joint_assurance, 1)
  expect_lt(result$assurance_lower, 0.9)
})

test_that("ss_unified(B = 1, decision = 'point') also reports a non-vacuous assurance_lower even though 'point' does not use it to decide", {
  # decision = "point" accepts on joint_assurance alone, so assurance_lower
  # does not drive its acceptance -- but the field is still populated and
  # reported, and must not mislead a caller who reads it anyway.
  result <- suppressWarnings(ss_unified(
    prior_se = c(40, 2), prior_sp = c(40, 2), prior_prev = c(10, 10),
    Se_ref = 1, Sp_ref = 1, delta_se = 0.30, delta_sp = 0.30, delta_auc = 0,
    check_nb = FALSE, target_assurance = 0.5, decision = "point",
    N_range = seq(50, 400, by = 10), B = 1, seed = 2026
  ))
  expect_equal(result$joint_assurance, 1)
  expect_lt(result$assurance_lower, 0.9)
})

test_that("ss_unified(B = 1, decision = 'isotonic') (the default) also avoids a vacuous margin-adjusted bound (defect 5)", {
  # Isotonic pooling raises the effective sample size of a flat block above
  # B = 1 (several consecutive N's sharing the same fitted value), so a
  # margin far below 1.000 does not by itself prove the pooling is
  # trustworthy at this B: a WIDE pooled block built from mostly-coincidental
  # ties (see the next test) could still report a deceptively strong bound
  # -- e.g. close to 0.9 -- while still clearing the old, near-vacuous
  # 0.999999 threshold. `assert_lt(..., 0.999999)` alone therefore passed
  # even with the pre-fix Wald formula's exact 1.0 minus a 1e-6 slack, and
  # continued to pass with a merely LARGE (not vacuous) Wilson bound after
  # the Wald-to-Wilson fix -- neither is what "not vacuous" should mean.
  # The block-pooling multiplier is now additionally capped at
  # min(block_len, B) (defect 5): at B = 1 no block's effective sample
  # size can exceed 1 * 1 = 1, so assurance_lower can never clear
  # wilson_lower(1, 1, qnorm(0.95)) =~ 0.270, regardless of how many grid
  # points end up pooled into the accepting block.
  result <- suppressWarnings(ss_unified(
    prior_se = c(40, 2), prior_sp = c(40, 2), prior_prev = c(10, 10),
    Se_ref = 1, Sp_ref = 1, delta_se = 0.30, delta_sp = 0.30, delta_auc = 0,
    check_nb = FALSE, target_assurance = 0.5, decision = "isotonic",
    N_range = seq(50, 400, by = 10), B = 1, seed = 2026
  ))
  expect_equal(result$joint_assurance, 1)
  cap_bound <- wilson_lower(1, 1, stats::qnorm(0.95))
  expect_lt(result$assurance_lower, cap_bound + 1e-8)
})

test_that("a minimal B cannot manufacture a strong isotonic margin from a wide, coarse curve (defect 5)", {
  # The REAL hole defect 5 closes: rle(fitted) cannot tell a block that
  # isoreg() genuinely pooled to enforce monotonicity apart from a run of
  # already-DISTINCT, non-duplicated N whose RAW assurance simply happened
  # to coincide -- both look identical in `fitted`. At B = 1, raw
  # assurance can only take 2 values (0, 1), so long flat runs of distinct
  # N arise routinely by chance, not because the curve is genuinely flat.
  # Independently reproduced against the pre-defect-5 code (Wald-to-Wilson
  # fix already in place, block_len cap not yet added) on this exact
  # scenario: the search silently returned a converged N_effective deep
  # inside N_range with assurance_lower above 0.80 -- built from ONE bit
  # of information per contributing N -- and NO non-convergence warning.
  # With the pooling multiplier capped at min(block_len, B), a B = 1
  # search can no longer manufacture a bound anywhere near
  # target_assurance = 0.80 from this data, so it now correctly falls back
  # to the generic "did not converge" warning instead.
  grid <- round(seq(50, 3050, length.out = 101))
  expect_warning(
    result <- ss_unified(
      prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
      Se_ref = 0.90, Sp_ref = 0.95, delta_se = 0.07, delta_sp = 0.05,
      delta_auc = 0, check_nb = FALSE, target_assurance = 0.80,
      decision = "isotonic", N_range = grid, B = 1, seed = 2026,
      full_grid = TRUE
    ),
    "Consider expanding N_range"
  )
  # Whatever is reported, it must not resemble a strong bound: at B = 1 the
  # cap limits ANY block's effective sample size to 1.
  expect_lt(result$assurance_lower, 0.5)
})

test_that("B = 1's assurance_lower is not vacuous, but B = 20000 (the manuscript's B) is essentially unaffected by the formula change", {
  # Same scenario, realistic B: the Wilson and Wald bounds must be close
  # enough that this package's own published results (computed at
  # B = 20000) do not move.
  z <- stats::qnorm(0.95)
  phat <- 0.85  # a representative joint_assurance value, not from a search
  wald <- phat - z * sqrt(phat * (1 - phat) / 20000)
  wilson <- wilson_lower(phat, 20000, z)
  expect_equal(wilson, wald, tolerance = 1e-4)
})
