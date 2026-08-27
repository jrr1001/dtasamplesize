
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
  expect_warning(
    ss_time_dependent_roc(B = 30, seed = 42, N_range = c(100, 150),
                          censoring_rates = 0.20),
    "target precision not"
  )
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
