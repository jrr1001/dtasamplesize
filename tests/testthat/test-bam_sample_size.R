
test_that("BAM n_diseased is near Buderer with informative prior", {
  # n_range must reach n_sp (= 373 here), otherwise the Sp search falls back to
  # max(n_range) with a warning and this test would pass on that branch.
  expect_no_warning(
    result <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14,
                              B = 5000, seed = 2026, n_range = 20:600)
  )
  # With prior Beta(17,3) and delta_se=0.14, BAM gives ~107 diseased
  expect_true(result$n_diseased > 90 & result$n_diseased < 130)
  # Converged strictly inside the range (not pinned to the ceiling)
  expect_lt(result$n_diseased, 600L)
  expect_lt(result$n_non_diseased, 600L)
})

test_that("BAM assurance meets target", {
  expect_no_warning(
    result <- bam_sample_size(B = 5000, seed = 2026, n_range = 20:600)
  )
  expect_true(result$assurance_se >= 0.80)
  expect_true(result$assurance_sp >= 0.80)
})

test_that("BAM returns dtasamplesize class with expected fields", {
  expect_no_warning(
    result <- bam_sample_size(B = 1000, seed = 42, n_range = 20:600)
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(!is.null(result$n_diseased))
  expect_true(!is.null(result$n_non_diseased))
  expect_true(!is.null(result$N_total_median))
  expect_true(!is.null(result$assurance_se))
  expect_true(!is.null(result$assurance_sp))
})

test_that("BAM vague prior requires more subjects than informative", {
  expect_no_warning(
    res_informative <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14,
                                       B = 3000, seed = 2026, n_range = 20:600)
  )
  expect_no_warning(
    res_vague <- bam_sample_size(prior_se = c(2, 2), delta_se = 0.14,
                                 B = 3000, seed = 2026, n_range = 20:600)
  )
  expect_true(res_vague$n_diseased > res_informative$n_diseased)
})

test_that("BAM total N covers BOTH arms (Se and Sp requirements)", {
  result <- bam_sample_size(B = 2000, seed = 2026, n_range = 20:600)
  # The total must supply n_se diseased AND n_sp non-diseased. With the vague
  # default prior_sp the specificity arm is the binding one.
  expect_gte(result$N_total_median,
             ceiling(result$n_non_diseased / (1 - 0.30)))
  expect_true(is.finite(result$N_total_P90))
})

test_that("BAM still warns when n_range is genuinely too short", {
  # The fallback branch must keep announcing itself. n_range = 20:150 lets the
  # Se search converge (n_se ~ 107) but not the Sp one (n_sp ~ 373), so exactly
  # one warning is raised.
  expect_warning(
    res <- bam_sample_size(B = 500, seed = 2026, n_range = 20:150),
    "achieved target assurance for Sp"
  )
  expect_equal(res$n_non_diseased, 150L)
  expect_lt(res$n_diseased, 150L)   # the Se arm did converge
})
