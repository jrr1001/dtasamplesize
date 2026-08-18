
test_that("Unified N exceeds Buderer", {
  expect_no_warning(
    result <- ss_unified(B = 1000, seed = 2026,
                         N_range = seq(200, 1200, by = 50),
                         delta_auc = 0, check_nb = FALSE)
  )
  buderer_N <- ceiling(buderer_n(0.85, 0.07) / (5 / 20))
  expect_true(result$N_effective > buderer_N)
})

test_that("Unified runs with AUC enabled", {
  expect_no_warning(
    result <- ss_unified(B = 500, seed = 2026,
                         N_range = seq(300, 1200, by = 50),
                         delta_auc = 0.06, check_nb = FALSE)
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(result$joint_assurance >= 0.80)
})

test_that("Unified runs with NB enabled and informative Sp prior", {
  # The CI-based NB criterion is genuinely demanding: N_range must reach ~1600.
  expect_no_warning(
    result <- ss_unified(B = 500, seed = 2026,
                         prior_sp = c(18, 2),
                         N_range = seq(400, 2200, by = 100),
                         delta_auc = 0, check_nb = TRUE)
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(result$joint_assurance >= 0.80)
})

test_that("M-1: the CI-based NB criterion demands MORE N than no NB check", {
  # Under the old point-estimate criterion, check_nb = TRUE barely moved N,
  # because "NB_hat > 0 and NB_hat > NB_all" is satisfied with near-certainty
  # whenever the test is useful at all. With the CI-based criterion it bites.
  common <- list(B = 400, seed = 2026, prior_sp = c(18, 2),
                 delta_auc = 0, N_range = seq(400, 2400, by = 100))
  without <- do.call(ss_unified, c(common, list(check_nb = FALSE)))
  with_nb <- do.call(ss_unified, c(common, list(check_nb = TRUE)))
  expect_gt(with_nb$N_effective, without$N_effective)
})

test_that("M-2: the assurance denominator is B (degenerate draws are failures)", {
  # A tiny N with a low-prevalence prior produces many degenerate replications
  # (n_d < 5). Those used to be dropped from the denominator, yielding an
  # assurance CONDITIONAL on the study not being degenerate. Now they are
  # failures, so the assurance at a small N must be very low, not undefined
  # or inflated.
  res <- suppressWarnings(
    ss_unified(B = 500, seed = 2026, prior_prev = c(1, 40),
               N_range = c(30, 40), delta_auc = 0, check_nb = FALSE)
  )
  expect_true(res$joint_assurance >= 0)
  expect_true(res$joint_assurance < 0.10)
  expect_false(is.na(res$joint_assurance))
})

test_that("Unified returns comparison table", {
  expect_no_warning(
    result <- ss_unified(B = 500, seed = 42,
                         N_range = seq(200, 1200, by = 50),
                         delta_auc = 0, check_nb = FALSE)
  )
  expect_true(is.data.frame(result$comparison))
  expect_equal(nrow(result$comparison), 3)
  expect_match(result$comparison$method[2], "Rogan-Gladen")
})

# --- grid_results / full_grid -----------------------------------------

test_that("grid_results is a data frame with N and assurance columns", {
  result <- ss_unified(B = 300, seed = 2026,
                       N_range = seq(300, 900, by = 100),
                       delta_auc = 0, check_nb = FALSE)
  expect_true(is.data.frame(result$grid_results))
  expect_true(all(c("N", "assurance") %in% names(result$grid_results)))
  expect_true(nrow(result$grid_results) >= 1)
  expect_true(all(result$grid_results$assurance >= 0 &
                    result$grid_results$assurance <= 1))
})

test_that("full_grid = TRUE evaluates the whole N_range without changing N_effective", {
  common <- list(B = 300, seed = 2026, N_range = seq(300, 900, by = 100),
                 delta_auc = 0, check_nb = FALSE)
  truncated <- do.call(ss_unified, c(common, list(full_grid = FALSE)))
  full <- do.call(ss_unified, c(common, list(full_grid = TRUE)))

  # The default (full_grid = FALSE) stops the grid at the optimum, so it has
  # STRICTLY FEWER rows than the untruncated search over the same N_range.
  expect_gt(nrow(full$grid_results), nrow(truncated$grid_results))
  expect_equal(nrow(full$grid_results), length(common$N_range))

  # optimal_N / joint_assurance / n_total must be identical either way: they
  # always refer to the FIRST N to reach the target, never the last.
  expect_equal(full$N_effective, truncated$N_effective)
  expect_equal(full$n_total, truncated$n_total)
  expect_equal(full$joint_assurance, truncated$joint_assurance)
})

# --- N_buderer / N_imperfect / seed / target_assurance -----------------

test_that("N_buderer, N_imperfect, seed and target_assurance are present and coherent", {
  result <- ss_unified(B = 300, seed = 4242,
                       N_range = seq(300, 900, by = 100),
                       delta_auc = 0, check_nb = FALSE,
                       target_assurance = 0.75)

  expect_false(is.null(result$N_buderer))
  expect_false(is.null(result$N_imperfect))
  expect_equal(result$seed, 4242)
  expect_equal(result$target_assurance, 0.75)

  # These must agree exactly with the corresponding rows of $comparison,
  # since they are the same quantities computed once.
  expect_equal(
    result$N_buderer,
    result$comparison$N[result$comparison$method == "Buderer (classical)"]
  )
  expect_equal(
    result$N_imperfect,
    result$comparison$N[grepl("Imperfect ref", result$comparison$method)]
  )
})

test_that("target_assurance defaults to the formal default when not overridden", {
  result <- ss_unified(B = 300, seed = 2026,
                       N_range = seq(300, 900, by = 100),
                       delta_auc = 0, check_nb = FALSE)
  expect_equal(result$target_assurance, formals(ss_unified)[["target_assurance"]])
})
