
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
