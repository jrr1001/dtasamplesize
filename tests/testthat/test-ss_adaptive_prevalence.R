
test_that("Adaptive does not reduce N when prevalence is higher", {
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.42), B = 1000, seed = 2026)
  expect_true(result$results$N_final_median[1] >= result$N_initial)
})

test_that("Adaptive increases N when prevalence is lower", {
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.18), B = 1000, seed = 2026)
  expect_true(result$results$N_final_mean[1] > result$N_initial * 1.1)
})

test_that("ss_adaptive_prevalence returns correct class", {
  result <- ss_adaptive_prevalence(B = 200, seed = 42,
                                   prev_true_range = c(0.25, 0.30))
  expect_s3_class(result, "dtasamplesize")
  expect_true(is.data.frame(result$results))
  expect_equal(nrow(result$results), 2)
  expect_true("N_analysed_median" %in% names(result$results))
})

test_that("At initial prevalence, N_final is near N_initial", {
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.30), B = 1000, seed = 2026)
  # Median should be close to N_initial (within 20%)
  ratio <- result$results$N_final_median[1] / result$N_initial_adj
  expect_true(ratio >= 0.95 && ratio <= 1.30)
})

# --- losses must COST precision, never buy it -------------------------

test_that("precision_achieved does NOT increase with loss_rate", {
  # The bug: the inflated sample was analysed in full, so the inflation was
  # free extra data and precision ROSE with the loss rate. Now the losses
  # actually happen, so inflating merely COMPENSATES them and precision is
  # flat in loss_rate (up to Monte Carlo noise).
  lrs <- c(0, 0.10, 0.20, 0.30)
  prec <- vapply(lrs, function(lr) {
    ss_adaptive_prevalence(B = 4000, seed = 2026, loss_rate = lr,
                           prev_true_range = 0.30)$results$precision_achieved
  }, numeric(1))

  # MC standard error at B = 4000 and p ~ 0.63 is about 0.0076; allow 4 SE.
  tol <- 0.03
  expect_lt(prec[length(prec)] - prec[1], tol)   # no upward drift end-to-end
  expect_lt(max(prec) - min(prec), 2 * tol)      # essentially flat
  # And emphatically not the old behaviour (0.70 -> 0.98).
  expect_lt(prec[length(prec)] - prec[1], 0.10)
})

test_that("recruited N grows with loss_rate but ANALYSED N does not", {
  r0 <- ss_adaptive_prevalence(B = 2000, seed = 2026, loss_rate = 0,
                               prev_true_range = 0.30)$results
  r3 <- ss_adaptive_prevalence(B = 2000, seed = 2026, loss_rate = 0.30,
                               prev_true_range = 0.30)$results
  # You must enrol substantially more...
  expect_gt(r3$N_final_mean, r0$N_final_mean * 1.2)
  # ...to end up analysing about the same number of people.
  expect_lt(abs(r3$N_analysed_median - r0$N_analysed_median) /
              r0$N_analysed_median, 0.10)
})

# --- the interim prevalence estimate must be bounded ABOVE too --------

test_that("prev_hat is truncated above, so results stay finite", {
  # With prev_true near 1 the stage-1 estimate hits 1, and the old code
  # computed n_sp / (1 - 1) = Inf, poisoning the whole results table.
  result <- suppressWarnings(ss_adaptive_prevalence(B = 300, seed = 1,
                                   prev_initial = 0.90, prev_true_range = 0.97))
  expect_true(all(is.finite(as.matrix(result$results))))
  expect_true(is.finite(result$results$N_final_mean))
  expect_true(is.finite(result$results$N_final_P75))
})

test_that("the truncation caps the re-estimated N at a known bound", {
  # prev_hat in [0.05, 0.95] => N_new <= max(n_se / 0.05, n_sp / 0.05)
  result <- suppressWarnings(
    ss_adaptive_prevalence(B = 300, seed = 3, prev_true_range = 0.02))
  cap <- ceiling(max(result$n_se / 0.05, result$n_sp / 0.05) / (1 - 0.10))
  expect_lte(result$results$N_final_P75, cap)
})

test_that("a truncated interim prevalence estimate emits a warning", {
  # At very low true prevalence the stage-1 estimate frequently falls below the
  # 0.05 floor, so the truncation binds and the caller must be warned that the
  # diseased arm may be under-sized (silent truncation would hide this).
  expect_warning(
    ss_adaptive_prevalence(B = 300, seed = 3, prev_true_range = 0.02),
    "truncated")
  # At a central prevalence the truncation does not bind: no such warning.
  expect_no_warning(
    ss_adaptive_prevalence(B = 300, seed = 2026, prev_true_range = 0.30))
})

# --- stage 1 must be a genuine INTERNAL pilot, not a discarded one -----
#
# Stark & Zapf (2020) describe an "internal pilot study": the n_stage1
# subjects recruited to re-estimate the prevalence are retained and are
# part of the final analysed sample. Before this fix, stage 1 was
# simulated and then thrown away: stage 2 redrew a brand-new sample of
# size N_final_adj from scratch, so n_stage1 subjects were recruited in
# real life but never counted in the reported N, and the interim
# prevalence estimate (prev_hat) was statistically independent of the
# analysed sample -- the opposite of what an internal pilot is.

test_that("reported N reflects the loss-exempt internal pilot, not a fully re-inflated cohort", {
  # Mathematical guarantee under the fix: n_stage1 pilot subjects are
  # loss-exempt (they are already fully observed -- that is what makes
  # D_stage1 usable to re-estimate the prevalence; see @details "Losses to
  # follow-up"), so only the stage-2 top-up is inflated for losses. When no
  # upward adaptation is triggered (prev_hat >= prev_initial, the majority
  # of replicates here since prev_true > prev_initial), the reported
  # N_final drops BELOW N_initial_adj = ceiling(N_initial / (1 -
  # loss_rate)), the bound that would result from inflating the WHOLE
  # cohort for losses.
  #
  # Under the pre-fix code this can never happen: stage 1 was discarded and
  # stage 2 alone was inflated as N_final_adj = ceiling(max(buderer_total_N(
  # prev_hat), N_initial) / (1 - loss_rate)), which is >= N_initial_adj in
  # EVERY replicate (N_final >= N_initial always), so the median can never
  # fall below N_initial_adj. This test fails against that code and passes
  # against the fix.
  result <- ss_adaptive_prevalence(
    prev_true_range = 0.35, prev_initial = 0.30, B = 3000, seed = 2026
  )
  expect_lt(result$results$N_final_median, result$N_initial_adj)
})

test_that("the analysed:recruited ratio exceeds (1 - loss_rate), because stage-1 subjects are not re-subjected to loss", {
  # If the ENTIRE reported N_final were subjected to loss (the pre-fix
  # behaviour, since stage 1 was discarded and stage 2 alone -- treated as
  # the whole recruited cohort -- absorbed the full loss draw), the
  # analysed:recruited ratio would sit at essentially exactly (1 -
  # loss_rate) = 0.90. Under the fix, the n_stage1 subjects already
  # contributing to the analysed sample (D_stage1 diseased subjects
  # included) are loss-exempt, so the ratio is systematically HIGHER than
  # 0.90.
  loss_rate <- 0.10
  result <- ss_adaptive_prevalence(
    prev_true_range = 0.30, B = 5000, seed = 2026, loss_rate = loss_rate
  )
  ratio <- result$results$N_analysed_median / result$results$N_final_median
  expect_gt(ratio, (1 - loss_rate) + 0.02)
})

test_that("stage-1 diseased subjects (D_stage1) are folded into the analysed sample", {
  # A minimal, deterministic check that the analysed sample can never be
  # smaller than the stage-1 pilot itself: the n_stage1 pilot subjects are
  # always part of the analysed sample (N_analysis = n_stage1 +
  # n_stage2_analysis >= n_stage1), never discarded outright.
  result <- ss_adaptive_prevalence(
    prev_true_range = c(0.18, 0.42), B = 1000, seed = 2026
  )
  expect_true(all(result$results$N_analysed_median >= result$n_stage1))
})

test_that("the fraction_stage1 -> 1 edge case (N_final close to n_stage1) does not error", {
  # Cuidá el caso borde: when the stage-1 pilot is nearly as large as
  # N_initial, the re-estimated N_final can equal n_stage1, so the stage-2
  # top-up (N_final - n_stage1) must be clamped at 0 rather than going
  # negative.
  result <- suppressWarnings(ss_adaptive_prevalence(
    fraction_stage1 = 0.99, prev_true_range = c(0.30, 0.42),
    B = 300, seed = 2026
  ))
  expect_true(all(is.finite(as.matrix(result$results))))
  expect_true(all(result$results$N_final_median >= result$n_stage1))
})
