
test_that("VIF calculation is correct", {
  result <- ss_imperfect_ref(Se_ref = 0.90, Sp_ref = 0.95, B = 0)
  expected_vif <- 1 / (0.90 + 0.95 - 1)^2
  expect_equal(result$inflation_factor_se, expected_vif, tolerance = 0.001)
})

# --- a legal-but-useless reference standard must fail informatively ---

test_that("a near-coin-flip reference standard is refused, not exploded", {
  # Se_ref = Sp_ref = 0.51 passes the old check (Se_ref + Sp_ref > 1) but gives
  # Youden = 0.02 and VIF = 2500, which used to try to allocate ~9 Gb.
  expect_error(
    ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51),
    "Reference standard too weak"
  )
  # The message must name the Youden index and the VIF.
  expect_error(
    ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51),
    "Youden"
  )
  # It fails fast, i.e. before any Monte Carlo allocation, so B = 0 also errors.
  expect_error(
    ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51, B = 0),
    "Reference standard too weak"
  )
})

test_that("min_youden can be lowered deliberately, and the MC guard holds", {
  # Lowering min_youden lets the (huge) VIF through, but then the memory guard
  # on B * N_adjusted must catch the Monte Carlo allocation.
  expect_error(
    ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51, min_youden = 0.01,
                     B = 5000),
    "max_mc_cells"
  )
  # With B = 0 there is nothing to allocate, so it succeeds and reports the VIF.
  res <- ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51, min_youden = 0.01,
                          B = 0, sensitivity_table = FALSE)
  expect_equal(res$youden_ref, 0.02, tolerance = 1e-9)
  expect_equal(res$VIF, 1 / 0.02^2, tolerance = 1e-6)
})

test_that("the memory guard also catches an extreme prevalence", {
  # A fine reference standard (Youden 0.85) but prev = 0.001 blows N_adjusted
  # up on its own; the Youden bound alone would not catch this.
  expect_error(
    ss_imperfect_ref(prev = 0.001, B = 5000),
    "max_mc_cells"
  )
})

test_that("the method is named Rogan-Gladen, not Staquet", {
  result <- ss_imperfect_ref(B = 0)
  expect_match(result$method, "Rogan-Gladen")
  expect_false(grepl("Staquet", result$method, fixed = TRUE))
})

test_that("Default parameters stay comfortably inside both guards", {
  expect_error(ss_imperfect_ref(B = 0), NA)
  expect_error(ss_imperfect_ref(B = 5000, sensitivity_table = FALSE), NA)
})

test_that("Perfect reference gives VIF = 1", {
  result <- ss_imperfect_ref(Se_ref = 1.0, Sp_ref = 1.0, B = 0)
  expect_equal(result$inflation_factor_se, 1.0)
  expect_equal(result$n_diseased_unadjusted, result$n_diseased_adjusted)
})

test_that("Adjusted n_se is approximately 139 for default params", {
  result <- ss_imperfect_ref(Se_ref = 0.90, Sp_ref = 0.95, B = 0)
  # VIF ≈ 1.384, n_unadj_se = 100, so n_adj_se = ceiling(100 * 1.384) = 139
  expect_equal(result$n_diseased_adjusted, ceiling(100 * 1 / (0.85)^2))
})

test_that("Loss adjustment inflates total N", {
  result <- ss_imperfect_ref(B = 0, loss_rate = 0.10)
  expect_true(result$N_adjusted_loss > result$N_adjusted)
  expect_equal(result$N_adjusted_loss, ceiling(result$N_adjusted / 0.90))
})

test_that("Sensitivity table is generated", {
  result <- ss_imperfect_ref(B = 0, sensitivity_table = TRUE)
  expect_true(is.data.frame(result$sensitivity_table))
  expect_true(nrow(result$sensitivity_table) > 0)
  expect_true("VIF" %in% names(result$sensitivity_table))
})

test_that("MC validation shows adjusted improves precision on unadjusted", {
  result <- ss_imperfect_ref(B = 3000, seed = 2026, sensitivity_table = FALSE)
  mc <- result$mc_validation
  expect_true(is.data.frame(mc))
  expect_equal(nrow(mc), 2)
  # The inflated (adjusted) sample size achieves the target CI width more
  # often than the unadjusted one.
  adj <- mc[mc$scenario == "adjusted", ]
  unadj <- mc[mc$scenario == "unadjusted", ]
  expect_true(adj$P_width_target > unadj$P_width_target)
})

test_that("MC validation exposes the imperfect-reference bias", {
  result <- ss_imperfect_ref(B = 3000, seed = 2026, sensitivity_table = FALSE)
  mc <- result$mc_validation
  expect_true(all(c("se_apparent", "se_true", "bias") %in% names(mc)))
  # Apparent Se against an imperfect reference is biased below true Se.
  expect_true(all(mc$bias < 0))
})

test_that("ss_imperfect_ref returns dtasamplesize class", {
  result <- ss_imperfect_ref(B = 0)
  expect_s3_class(result, "dtasamplesize")
})
