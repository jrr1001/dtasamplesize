test_that("BAM n_diseased is near Buderer with informative prior", {
  result <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14,
                            B = 5000, seed = 2026, n_range = 20:200)
  # With prior Beta(17,3) and delta_se=0.14, BAM gives ~107 diseased
  expect_true(result$n_diseased > 90 & result$n_diseased < 130)
})

test_that("BAM assurance meets target", {
  result <- bam_sample_size(B = 5000, seed = 2026)
  expect_true(result$assurance_se >= 0.80)
  expect_true(result$assurance_sp >= 0.80)
})

test_that("BAM returns dtasamplesize class with expected fields", {
  result <- bam_sample_size(B = 1000, seed = 42)
  expect_s3_class(result, "dtasamplesize")
  expect_true(!is.null(result$n_diseased))
  expect_true(!is.null(result$n_non_diseased))
  expect_true(!is.null(result$N_total_median))
  expect_true(!is.null(result$assurance_se))
  expect_true(!is.null(result$assurance_sp))
})

test_that("BAM vague prior requires more subjects than informative", {
  res_informative <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14,
                                     B = 3000, seed = 2026, n_range = 20:300)
  res_vague <- bam_sample_size(prior_se = c(2, 2), delta_se = 0.14,
                               B = 3000, seed = 2026, n_range = 20:500)
  expect_true(res_vague$n_diseased > res_informative$n_diseased)
})
