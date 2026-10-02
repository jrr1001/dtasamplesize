
test_that("ss_time_dependent_roc runs without error", {
  skip_if_not_installed("timeROC")
  expect_no_warning(
    result <- ss_time_dependent_roc(B = 50, seed = 42,
      N_range = seq(100, 600, by = 100),
      censoring_rates = 0.20)
  )
  expect_s3_class(result, "dtasamplesize")
  # Converged strictly inside the grid, not pinned to its ceiling.
  expect_lt(result$results$N_required, 600)
  expect_gte(result$results$prob_achieved, 0.80)
})

test_that("N increases with censoring rate", {
  skip_if_not_installed("timeROC")
  result <- ss_time_dependent_roc(B = 100, seed = 2026,
    N_range = seq(100, 600, by = 50),
    censoring_rates = c(0.10, 0.30))
  ns <- result$results$N_required
  expect_true(ns[2] >= ns[1])
})

test_that("A short N_range still warns and does not silently converge", {
  skip_if_not_installed("timeROC")
  w <- testthat::capture_warnings(
    result <- ss_time_dependent_roc(B = 30, seed = 42, N_range = c(100, 150),
                          censoring_rates = 0.20)
  )
  expect_true(any(grepl("target precision not", w)))
  # Non-crossing contract (bam_sample_size()/joint_sample_size()/
  # ss_net_benefit(), NEWS 0.6.6/0.6.7): a non-crossing search must
  # report NA, never max(N_range) as though it were a validated design.
  # Versions <= 0.6.6 set N_required/n_total to max(N_range) (= 150) here.
  expect_false(result$target_reached)
  expect_false(result$results$target_reached)
  expect_true(is.na(result$results$N_required))
  expect_true(is.na(result$results$prob_achieved))
  expect_true(is.na(result$n_total))
  expect_true(is.na(result$n_diseased))
  expect_equal(result$results$N_at_max_prob, 150L)
})

# --- regression: non-crossing must never report max(N_range) as a real
# solution (confirmed against the installed 0.6.6) ------

test_that("non-crossing across ALL censoring_rates: n_total/n_diseased are NA, never max(N_range)", {
  skip_if_not_installed("timeROC")
  # ss_time_dependent_roc() raises multiple warnings here (one per
  # censoring_rate that fails to converge, plus one top-level warning);
  # capture all of them and require the top-level one.
  w <- testthat::capture_warnings(
    result <- ss_time_dependent_roc(
      delta_auc = 0.001, target_prob = 0.999,
      N_range = seq(100, 140, by = 20), B = 20, seed = 2026
    )
  )
  expect_true(any(grepl("n_total/n_diseased are NA", w)))
  expect_false(result$target_reached)
  expect_true(is.na(result$n_total))
  expect_true(is.na(result$n_diseased))
  expect_true(all(!result$results$target_reached))
  expect_true(all(is.na(result$results$N_required)))
  expect_true(all(is.na(result$results$prob_achieved)))
  # The pre-fix (<= 0.6.6) defect: n_total silently became max(N_range).
  expect_false(isTRUE(result$n_total == max(seq(100, 140, by = 20))))
})

test_that("non-crossing for only SOME censoring_rates: n_total is still NA overall, not the smaller converged max", {
  skip_if_not_installed("timeROC")
  # One easy rate (converges well inside the grid) and one essentially
  # impossible rate (delta_auc too tight for any N here): the overall
  # n_total must be NA because the true worst case is unknown, not the
  # (necessarily smaller and misleading) max over the rate that DID
  # converge.
  w <- testthat::capture_warnings(
    result <- ss_time_dependent_roc(
      delta_auc = 0.06, target_prob = 0.80,
      censoring_rates = c(0.001, 0.9999),
      N_range = seq(100, 200, by = 50), B = 40, seed = 2026
    )
  )
  expect_true(length(w) >= 1)
  expect_false(result$target_reached)
  expect_true(is.na(result$n_total))
  expect_true(is.na(result$n_diseased))
  expect_false(all(result$results$target_reached)) # at least one failed
})

# --- regression: n_total must not depend on the caller's normal.kind
# (defect 9) -----------------------------------------------------------

test_that("ss_time_dependent_roc gives the same n_total under different normal.kind settings (defect 9)", {
  skip_if_not_installed("timeROC")
  # set.seed(seed, kind = "Mersenne-Twister") (added in 0.6.1) fixes the
  # UNIFORM generator explicitly, but this function is the one place in
  # the package that also draws from the NORMAL generator
  # (stats::rnorm(), for the simulated biomarker below), which
  # `normal.kind` controls independently of `kind`. Versions <= 0.6.1
  # left `normal.kind` inherited from the caller, so identical arguments
  # and seed gave n_total = 420 under the caller's normal.kind =
  # "Inversion" (R's default) but 400 under "Box-Muller". The caller's
  # own RNG state (kind, normal.kind, seed) is saved and restored here via
  # on.exit(), so this test does not leak its RNGkind() changes.
  old_kind <- RNGkind()
  on.exit(suppressWarnings(RNGkind(
    kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3]
  )), add = TRUE)

  n_by_kind <- vapply(
    c("Inversion", "Box-Muller", "Kinderman-Ramage", "Ahrens-Dieter"),
    function(nk) {
      RNGkind("Mersenne-Twister", nk)
      set.seed(1)
      suppressWarnings(ss_time_dependent_roc(B = 50, seed = 2026))$n_total
    },
    numeric(1)
  )

  expect_equal(unname(n_by_kind), rep(unname(n_by_kind[1]), length(n_by_kind)))
})
