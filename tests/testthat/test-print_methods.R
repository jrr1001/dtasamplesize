test_that("print works for mc_validate_buderer", {
  result <- mc_validate_buderer(B = 100, seed = 1)
  expect_output(print(result), "Monte Carlo Validation")
  expect_output(print(result), "n_diseased")
})

test_that("print works for bam_sample_size", {
  # n_range widened: n_se ~ 109 and n_sp ~ 372 here, so 50:100 exercised
  # the non-converged fallback branch.
  result <- bam_sample_size(B = 100, seed = 1, n_range = 20:600)
  expect_output(print(result), "Bayesian Assurance")
  expect_output(print(result), "n_diseased")
})

test_that("print works for joint_sample_size", {
  # N_range widened: the real N is ~570, so seq(100, 300, 50) exercised
  # the non-converged fallback branch.
  result <- joint_sample_size(B = 100, seed = 1,
                              N_range = seq(100, 800, by = 50))
  expect_output(print(result), "Joint Sample Size")
})

test_that("print works for ss_imperfect_ref", {
  result <- ss_imperfect_ref(B = 0)
  expect_output(print(result), "Imperfect Reference")
})

test_that("print works for ss_adaptive_prevalence", {
  result <- ss_adaptive_prevalence(B = 50, seed = 1,
                                   prev_true_range = 0.30)
  expect_output(print(result), "Adaptive")
})

test_that("print works for ss_net_benefit", {
  result <- ss_net_benefit(pt_range = 0.20, B = 100, seed = 1)
  expect_output(print(result), "Net Benefit")
})

test_that("print works for ss_unified", {
  # N_range widened: the real N is ~700, so seq(200, 400, 100) exercised
  # the non-converged fallback branch.
  result <- ss_unified(B = 50, seed = 1,
                       N_range = seq(200, 1000, by = 100),
                       delta_auc = 0, check_nb = FALSE)
  expect_output(print(result), "Unified")
})

test_that("print returns invisibly", {
  result <- mc_validate_buderer(B = 100, seed = 1)
  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
})
