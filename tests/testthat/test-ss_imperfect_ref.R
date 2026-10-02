
# --- estimand argument: default, validation, and "always both" ---------

test_that("estimand defaults to \"apparent\"", {
  result <- ss_imperfect_ref(B = 0)
  expect_equal(result$estimand, "apparent")
  expect_equal(result$N_adjusted, result$N_apparent)
  expect_equal(result$N_adjusted_loss, result$N_apparent_loss)
})

test_that("estimand = \"corrected\" selects the corrected numbers", {
  result <- ss_imperfect_ref(B = 0, estimand = "corrected")
  expect_equal(result$estimand, "corrected")
  expect_equal(result$N_adjusted, result$N_corrected)
  expect_equal(result$N_adjusted_loss, result$N_corrected_loss)
})

test_that("an invalid estimand is rejected", {
  expect_error(ss_imperfect_ref(B = 0, estimand = "bogus"))
})

test_that("both N_apparent and N_corrected are always computed, regardless of estimand", {
  by_apparent <- ss_imperfect_ref(B = 0, estimand = "apparent")
  by_corrected <- ss_imperfect_ref(B = 0, estimand = "corrected")
  # Same underlying numbers either way -- estimand only picks which pair
  # populates the generic N_adjusted/n_total aliases.
  expect_equal(by_apparent$N_apparent, by_corrected$N_apparent)
  expect_equal(by_apparent$N_corrected, by_corrected$N_corrected)
  expect_false(is.null(by_apparent$N_corrected))
  expect_false(is.null(by_corrected$N_apparent))
  # The two are materially different numbers (not the same quantity twice).
  expect_true(by_apparent$N_apparent != by_apparent$N_corrected)
  # The results table echoes both, with a flag for which was selected.
  expect_equal(sort(by_apparent$results$estimand), c("apparent", "corrected"))
  expect_true(by_apparent$results$selected[by_apparent$results$estimand == "apparent"])
  expect_true(by_corrected$results$selected[by_corrected$results$estimand == "corrected"])
})

# --- REGRESSION: the fixed defect ---------------------------------------
#
# Versions <= 0.5.0 multiplied the Buderer Se/Sp sample size by
# VIF = 1/(Se_ref+Sp_ref-1)^2, calling the result the misclassification-
# corrected Se/Sp sample size. VIF is the Rogan-Gladen PREVALENCE
# variance-inflation factor; applying it to Se/Sp substitutes both the
# estimand and the baseline. The targets below (N*Var(Se_corrected) and
# the required variance multiplier relative to Buderer) were verified
# independently by delta method, Fisher information, and Monte Carlo
# (40000 replicates), which agreed to within 1e-4. This test would fail
# under the old VIF-based code (the fields either don't exist, or -- for
# the multiplier -- old code implicitly forces
# multiplier_se_corrected == VIF).

test_that("REGRESSION: corrected-estimand variance matches the independently-verified reference values, and is NOT the retired VIF", {
  scenarios <- list(
    list(Se = .85, Sp = .90, prev = .30, Se_ref = .90, Sp_ref = .95,
         nvar_se = 0.8116, mult_se = 1.910),
    list(Se = .85, Sp = .90, prev = .20, Se_ref = .90, Sp_ref = .95,
         nvar_se = 1.5749, mult_se = 2.470),
    list(Se = .95, Sp = .95, prev = .05, Se_ref = .80, Sp_ref = .85,
         nvar_se = 99.643, mult_se = 104.9)
  )
  for (s in scenarios) {
    res <- ss_imperfect_ref(Se = s$Se, Sp = s$Sp, prev = s$prev,
                            Se_ref = s$Se_ref, Sp_ref = s$Sp_ref, B = 0)
    expect_equal(res$NVar_se_corrected, s$nvar_se, tolerance = 1e-3)
    expect_equal(res$multiplier_se_corrected, s$mult_se, tolerance = 2e-2)
    # the defect: old code implicitly claimed multiplier_se_corrected == VIF
    expect_true(abs(res$multiplier_se_corrected - res$VIF) > 0.05)
  }
})

test_that("REGRESSION: manuscript configuration (prev = 0.20) reproduces the audited N's", {
  res_app <- ss_imperfect_ref(prev = 0.20, B = 0)  # default estimand = apparent
  res_cor <- ss_imperfect_ref(prev = 0.20, B = 0, estimand = "corrected")
  expect_equal(res_app$N_apparent, 732)
  expect_equal(res_cor$N_corrected, 1235)

  # The manuscript's now-superseded number came from VIF-multiplying the
  # Buderer arms directly; reconstruct it here (without calling any
  # ss_imperfect_ref internals) purely to show it differs from both of the
  # correct estimand-specific numbers above.
  old_defective_N <- ceiling(max(
    ceiling(buderer_n(0.85, 0.07) * res_app$VIF) / 0.20,
    ceiling(buderer_n(0.90, 0.05) * res_app$VIF) / 0.80
  ))
  expect_equal(old_defective_N, 695)
  expect_false(res_app$N_apparent == old_defective_N)
  expect_false(res_cor$N_corrected == old_defective_N)
})

test_that("VIF is retained as an informative quantity but no longer sizes anything", {
  result <- ss_imperfect_ref(Se_ref = 0.90, Sp_ref = 0.95, B = 0)
  expected_vif <- 1 / (0.90 + 0.95 - 1)^2
  expect_equal(result$VIF, expected_vif, tolerance = 0.001)
  # N_apparent is the exact closed form, not n_unadjusted * VIF.
  expect_false(isTRUE(all.equal(
    result$n_refpos_apparent,
    ceiling(buderer_n(0.85, 0.07) * result$VIF)
  )))
  # N_corrected is the delta-method total N, not n_unadjusted * VIF either.
  expect_false(isTRUE(all.equal(
    result$n_se_corrected,
    ceiling(buderer_n(0.85, 0.07) * result$VIF)
  )))
})

# --- a legal-but-useless reference standard must fail informatively ---

test_that("a near-coin-flip reference standard is refused, not exploded", {
  # Se_ref = Sp_ref = 0.51 passes the old check (Se_ref + Sp_ref > 1) but gives
  # Youden = 0.02.
  expect_error(
    ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51),
    "Reference standard too weak"
  )
  # The message must name the Youden index.
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
  # Lowering min_youden lets the (huge) corrected-estimand N through, but
  # then the memory guard on B * N_adjusted must catch the Monte Carlo
  # allocation.
  expect_error(
    ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51, min_youden = 0.01,
                     estimand = "corrected", B = 5000),
    "max_mc_cells"
  )
  # With B = 0 there is nothing to allocate, so it succeeds and reports the
  # (informative-only) VIF and Youden index.
  res <- ss_imperfect_ref(Se_ref = 0.51, Sp_ref = 0.51, min_youden = 0.01,
                          B = 0, sensitivity_table = FALSE)
  expect_equal(res$youden_ref, 0.02, tolerance = 1e-9)
  expect_equal(res$VIF, 1 / 0.02^2, tolerance = 1e-6)
})

test_that("the memory guard catches the corrected estimand at low prevalence, even with a fine reference standard", {
  # Youden_ref = 0.85 comfortably clears the default min_youden = 0.5, so
  # the guard above does not fire; but for estimand = "corrected" a low
  # prevalence alone drives N_adjusted into the hundreds of thousands (see
  # "apparent stays bounded ... " test below), which the memory guard must
  # still catch.
  expect_error(
    ss_imperfect_ref(prev = 0.01, estimand = "corrected", B = 100),
    "max_mc_cells"
  )
})

test_that("apparent stays bounded at low prevalence; corrected does not", {
  # This is the qualitative signature of the fixed defect: the apparent
  # estimand's variance is driven by P(ref+)/P(ref-), which stay bounded
  # away from 0 by the reference standard's own false-positive/negative
  # floor, while the corrected estimand's variance scales like
  # 1/(prev*Youden)^2 and explodes as prevalence shrinks.
  r_low <- ss_imperfect_ref(prev = 0.001, B = 0)
  expect_true(r_low$N_apparent < 5000)
  expect_true(r_low$N_corrected > 1e6)
  expect_true(r_low$N_corrected > 1000 * r_low$N_apparent)
})

test_that("the method no longer claims Staquet, and is not silent about the correction", {
  result <- ss_imperfect_ref(B = 0)
  expect_false(grepl("Staquet", result$method, fixed = TRUE))
})

test_that("Default parameters stay comfortably inside both guards", {
  expect_error(ss_imperfect_ref(B = 0), NA)
  expect_error(ss_imperfect_ref(B = 5000, sensitivity_table = FALSE), NA)
  expect_error(ss_imperfect_ref(B = 5000, sensitivity_table = FALSE,
                                estimand = "corrected"), NA)
})

test_that("Perfect reference: both estimands reduce exactly to Buderer", {
  # Se = 0.85, Sp = 0.90, prev = 0.30 (package defaults, not overridden here).
  result <- ss_imperfect_ref(Se_ref = 1.0, Sp_ref = 1.0, B = 0)
  expect_equal(result$VIF, 1.0)
  expect_equal(result$Se_apparent, 0.85, tolerance = 1e-9)
  expect_equal(result$n_refpos_apparent, result$n_diseased_unadjusted)
  expect_equal(result$n_refneg_apparent, result$n_nondiseased_unadjusted)
  expect_equal(result$N_apparent, result$N_unadjusted)
  # The corrected estimand also reduces to Buderer's variance exactly when
  # the reference is perfect (Se_ref = Sp_ref = 1 makes Se_hat/Sp_hat the
  # naive proportions).
  expect_equal(result$NVar_se_corrected, 0.85 * 0.15 / 0.30, tolerance = 1e-9)
  expect_equal(result$NVar_sp_corrected, 0.90 * 0.10 / 0.70, tolerance = 1e-9)
})

test_that("Loss adjustment inflates total N, for both estimands", {
  result <- ss_imperfect_ref(B = 0, loss_rate = 0.10)
  expect_true(result$N_adjusted_loss > result$N_adjusted)
  expect_equal(result$N_adjusted_loss, ceiling(result$N_adjusted / 0.90))
  expect_equal(result$N_apparent_loss, ceiling(result$N_apparent / 0.90))
  expect_equal(result$N_corrected_loss, ceiling(result$N_corrected / 0.90))
})

test_that("Sensitivity table is generated with both estimands' N", {
  result <- ss_imperfect_ref(B = 0, sensitivity_table = TRUE)
  tbl <- result$sensitivity_table
  expect_true(is.data.frame(tbl))
  expect_true(nrow(tbl) > 0)
  expect_true(all(c("VIF", "N_apparent", "N_corrected", "N_adj") %in% names(tbl)))
  # N_adj tracks the selected estimand (apparent, here).
  expect_equal(tbl$N_adj, tbl$N_apparent)
})

test_that("MC validation shows the selected estimand's N improves apparent-Se precision over Buderer's", {
  result <- ss_imperfect_ref(B = 3000, seed = 2026, sensitivity_table = FALSE)
  mc <- result$mc_validation
  expect_true(is.data.frame(mc))
  expect_equal(nrow(mc), 2)
  adj <- mc[mc$scenario == "adjusted", ]
  unadj <- mc[mc$scenario == "unadjusted", ]
  expect_true(adj$P_width_target > unadj$P_width_target)
})

test_that("MC validation exposes the imperfect-reference bias under apparent, and its (near) absence under the corrected estimator", {
  result <- ss_imperfect_ref(B = 5000, seed = 2026, sensitivity_table = FALSE)
  mc <- result$mc_validation
  expect_true(all(c("se_apparent", "se_true", "bias",
                     "se_corrected", "bias_corrected") %in% names(mc)))
  # Apparent Se against an imperfect reference is biased below true Se.
  expect_true(all(mc$bias < 0))
  # The misclassification-corrected estimator has (at most) a small-sample
  # bias, an order of magnitude smaller than the apparent-Se bias.
  expect_true(all(abs(mc$bias_corrected) < abs(mc$bias) / 5))
})

test_that("ss_imperfect_ref returns dtasamplesize class", {
  result <- ss_imperfect_ref(B = 0)
  expect_s3_class(result, "dtasamplesize")
})

# --- delta_se / delta_sp aliases of d_se / d_sp -------------------------

test_that("delta_se/delta_sp are equivalent aliases of d_se/d_sp", {
  by_d <- ss_imperfect_ref(d_se = 0.07, d_sp = 0.05, B = 0)
  by_delta <- ss_imperfect_ref(delta_se = 0.07, delta_sp = 0.05, B = 0)
  expect_identical(by_d$n_diseased_adjusted, by_delta$n_diseased_adjusted)
  expect_identical(by_d$n_nondiseased_adjusted, by_delta$n_nondiseased_adjusted)
  expect_identical(by_d$N_adjusted, by_delta$N_adjusted)
})

test_that("delta_se/delta_sp take precedence over d_se/d_sp when both are supplied", {
  result <- ss_imperfect_ref(d_se = 0.20, d_sp = 0.20,
                             delta_se = 0.07, delta_sp = 0.05, B = 0)
  expected <- ss_imperfect_ref(d_se = 0.07, d_sp = 0.05, B = 0)
  expect_identical(result$n_diseased_adjusted, expected$n_diseased_adjusted)
  expect_identical(result$n_nondiseased_adjusted, expected$n_nondiseased_adjusted)
})

# --- regression: mc_validation must not depend on the caller's
# normal.kind/sample.kind (defect 9) -------------------------------------

test_that("ss_imperfect_ref's mc_validation is unaffected by the caller's normal.kind (defect 9)", {
  # set.seed(seed) (no `kind`) left this function's Monte Carlo bias check
  # dependent on the caller's ACTIVE generator entirely -- not just
  # normal.kind, but the uniform kind too, since no kind was named at all.
  # Reproduced directly with the exact call used in
  # validation/reproduce_manuscript.R: the manuscript's published apparent
  # sensitivity of 0.764 (bias -0.086) reproduced under Mersenne-Twister
  # but rounded to 0.765 under Wichmann-Hill. The caller's own RNG state
  # is saved and restored here via on.exit(), so this test does not leak
  # its RNGkind() changes.
  old_kind <- RNGkind()
  on.exit(suppressWarnings(RNGkind(
    kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3]
  )), add = TRUE)

  gens <- c("Mersenne-Twister", "Wichmann-Hill", "Marsaglia-Multicarry",
            "Super-Duper", "Knuth-TAOCP-2002", "L'Ecuyer-CMRG")
  se_by_gen <- vapply(gens, function(g) {
    suppressWarnings(RNGkind(g))
    set.seed(1)
    ir <- ss_imperfect_ref(B = 6000, seed = 2026, sensitivity_table = FALSE)
    ir$mc_validation$se_apparent[ir$mc_validation$scenario == "adjusted"]
  }, numeric(1))

  expect_equal(unname(se_by_gen), rep(unname(se_by_gen[1]), length(se_by_gen)))
  # Published-manuscript figure (R <= 4.5 sampler): within +/- 3 MCSE at
  # B = 6000, robust to the RNG-stream change on R-devel/4.6.x (see NEWS
  # 0.6.7); exact rounded check kept as a guarded
  # sanity check on R < 4.6.0.
  p_hat <- unname(se_by_gen[1])
  mcse <- sqrt(p_hat * (1 - p_hat) / 6000)
  expect_lt(abs(p_hat - 0.764), 3 * mcse)
  if (getRversion() < "4.6.0") {
    expect_equal(unname(round(se_by_gen[1], 3)), 0.764)
  }
})
