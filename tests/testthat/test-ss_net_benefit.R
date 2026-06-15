test_that("ss_net_benefit returns correct class and fields", {
  result <- ss_net_benefit(pt_range = c(0.10, 0.30), B = 1000, seed = 42)
  expect_s3_class(result, "dtasamplesize")
  expect_true(is.data.frame(result$N_by_pt))
  expect_equal(nrow(result$N_by_pt), 2)
  expect_true(all(c("pt", "feasible", "N_required", "prob_achieved",
                    "NB_true", "NB_treat_all") %in% names(result$N_by_pt)))
  expect_true(!is.null(result$N_conservative))
})

test_that("Useful thresholds are feasible and return a sample size", {
  result <- ss_net_benefit(B = 2000, seed = 2026)
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

test_that("Required N grows as assurance target increases", {
  lo <- ss_net_benefit(pt_range = 0.50, target_prob = 0.70,
                       B = 3000, seed = 2026)
  hi <- ss_net_benefit(pt_range = 0.50, target_prob = 0.90,
                       B = 3000, seed = 2026)
  expect_true(hi$N_by_pt$N_required[1] >= lo$N_by_pt$N_required[1])
})
