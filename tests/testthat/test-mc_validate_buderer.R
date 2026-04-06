test_that("mc_validate_buderer runs without error", {
  result <- mc_validate_buderer(Se = 0.85, d = 0.07, B = 500, seed = 42)
  expect_s3_class(result, "dtasamplesize")
  expect_equal(result$n_diseased, 100)
})

test_that("Wald coverage is below nominal", {
  result <- mc_validate_buderer(Se = 0.85, d = 0.07, B = 5000, seed = 2026)
  wald_row <- result$results[result$results$ci_method == "wald", ]
  # Wald coverage should be approximately 92-94% (below 95%)
  expect_true(wald_row$coverage < 0.95)
  expect_true(wald_row$coverage > 0.90)
  # P(width <= target) should be approximately 55-62%
  expect_true(wald_row$P_width_target > 0.50)
  expect_true(wald_row$P_width_target < 0.65)
})

test_that("Wilson has tighter P90 width than Wald", {
  result <- mc_validate_buderer(Se = 0.85, d = 0.07, B = 5000, seed = 2026)
  wald <- result$results[result$results$ci_method == "wald", ]
  wilson <- result$results[result$results$ci_method == "wilson", ]
  # Wilson produces tighter upper-tail widths
  expect_true(wilson$P90_width < wald$P90_width)
})

test_that("Wilson outperforms Wald near boundary", {
  # At Se=0.95, Wilson coverage advantage is clear
  result <- mc_validate_buderer(Se = 0.95, d = 0.05, B = 5000, seed = 2026)
  wald <- result$results[result$results$ci_method == "wald", ]
  wilson <- result$results[result$results$ci_method == "wilson", ]
  expect_true(wilson$coverage > wald$coverage)
})

test_that("ci_method parameter works", {
  res_wald <- mc_validate_buderer(B = 500, seed = 1, ci_method = "wald")
  expect_equal(nrow(res_wald$results), 1)
  expect_equal(res_wald$results$ci_method, "wald")

  res_wilson <- mc_validate_buderer(B = 500, seed = 1, ci_method = "wilson")
  expect_equal(nrow(res_wilson$results), 1)
  expect_equal(res_wilson$results$ci_method, "wilson")

  res_both <- mc_validate_buderer(B = 500, seed = 1, ci_method = "both")
  expect_equal(nrow(res_both$results), 2)
})
