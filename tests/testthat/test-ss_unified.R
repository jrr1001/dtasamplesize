
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
  # At this B, the first N to clear the default lower-bound decision rule
  # lands close enough to its own Monte Carlo error to trip the new
  # insufficient-B warning (see the "close to its own Monte Carlo noise"
  # tests below); that warning is expected here and is not what this test
  # is about, so it is suppressed rather than asserted away.
  result <- suppressWarnings(
    ss_unified(B = 500, seed = 2026,
              prior_sp = c(18, 2),
              N_range = seq(400, 2200, by = 100),
              delta_auc = 0, check_nb = TRUE, nb_B_ceiling = 50000)
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(result$joint_assurance >= 0.80)
})

test_that("M-1: the CI-based NB criterion demands MORE N than no NB check", {
  # Under the old point-estimate criterion, check_nb = TRUE barely moved N,
  # because "NB_hat > 0 and NB_hat > NB_all" is satisfied with near-certainty
  # whenever the test is useful at all. With the CI-based criterion it bites.
  # (B = 400 also lands close enough to the decision boundary to trip the
  # insufficient-B warning; suppressed here since this test is about
  # N_effective, not about that warning.)
  common <- list(B = 400, seed = 2026, prior_sp = c(18, 2),
                 delta_auc = 0, N_range = seq(400, 2400, by = 100),
                 nb_B_ceiling = 50000)
  without <- suppressWarnings(
    do.call(ss_unified, c(common, list(check_nb = FALSE)))
  )
  with_nb <- suppressWarnings(
    do.call(ss_unified, c(common, list(check_nb = TRUE)))
  )
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
  # The imperfect-reference correction contributes two rows, not one: the
  # apparent estimand (exact closed form) and the corrected estimand
  # (delta-method misclassification correction) cannot be collapsed into a
  # single "inflation factor" row (see ss_imperfect_ref()'s Details).
  expect_equal(nrow(result$comparison), 4)
  expect_match(result$comparison$method[2], "apparent")
  expect_match(result$comparison$method[3], "corrected")
  expect_false(any(grepl("Rogan-Gladen", result$comparison$method)))
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
  # decision = "lower_bound": this test is about full_grid's early-stopping
  # mechanics for a first-crossing rule. decision = "isotonic" (the default
  # since 0.5.1) always evaluates the complete N_range regardless of
  # full_grid, which would make truncated/full identical here and defeat
  # the point of the test.
  common <- list(B = 300, seed = 2026, N_range = seq(300, 900, by = 100),
                 delta_auc = 0, check_nb = FALSE, decision = "lower_bound")
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
  # since they are the same quantities computed once. N_imperfect is
  # specifically the APPARENT-estimand row (see ?ss_unified); the
  # corrected-estimand row is a separate, larger N, also in $comparison.
  expect_equal(
    result$N_buderer,
    result$comparison$N[result$comparison$method == "Buderer (classical)"]
  )
  expect_equal(
    result$N_imperfect,
    result$comparison$N[grepl("Imperfect ref.*apparent", result$comparison$method)]
  )
})

test_that("target_assurance defaults to the formal default when not overridden", {
  result <- ss_unified(B = 300, seed = 2026,
                       N_range = seq(300, 900, by = 100),
                       delta_auc = 0, check_nb = FALSE)
  expect_equal(result$target_assurance, formals(ss_unified)[["target_assurance"]])
})

# --- decision = "isotonic" (default) vs "lower_bound"/"point" (legacy) -

test_that("decision defaults to \"isotonic\" and is echoed in the result", {
  # isotonic became the default in 0.5.1: measured against a B = 1e7
  # reference N* (the step-5 scenario of Figure 4), it has far lower
  # seed-to-seed sd, bias and RMSE than lower_bound's first-crossing rule,
  # and (unlike lower_bound) actually delivers its guarantee across seeds.
  # See ss_unified()'s @details for the measured numbers.
  expect_equal(formals(ss_unified)[["decision"]],
               quote(c("isotonic", "lower_bound", "point")))
  result <- ss_unified(B = 300, seed = 2026,
                       N_range = seq(300, 900, by = 100),
                       delta_auc = 0, check_nb = FALSE)
  expect_equal(result$decision, "isotonic")
})

test_that("an invalid decision is rejected by match.arg", {
  expect_error(
    ss_unified(decision = "bogus", B = 100, N_range = seq(300, 900, by = 200)),
    "'arg' should be one of"
  )
})

test_that("decision = \"point\" reproduces the pre-0.5.0 behaviour (N_effective = 920)", {
  # Independently verified: at B = 1200, seed = 2026, this exact search
  # (prior_prev = c(4, 16), delta_auc = 0, check_nb = FALSE) selects
  # N_effective = 920 under the point-estimate stopping rule. This is NOT
  # the N = 900 reported by versions <= 0.5.0: 0.5.1 replaced the per-N
  # set.seed(seed) pattern (only PARTIAL common random numbers -- rbinom's
  # uniform consumption depends on N, so streams diverged from the second
  # replication on) with genuinely independent L'Ecuyer-CMRG streams per N
  # (see @details), which changes the simulated draws at every N even
  # though `seed` is unchanged. The near-threshold warning is expected
  # here (the point estimate at N = 920 is close to its own MCSE, which is
  # exactly the situation decision = "lower_bound" exists to guard
  # against) and is not the object of this test.
  result <- suppressWarnings(ss_unified(
    decision = "point", B = 1200, seed = 2026, prior_prev = c(4, 16),
    N_range = seq(200, 1600, by = 20), delta_auc = 0, check_nb = FALSE
  ))
  expect_equal(result$decision, "point")
  expect_equal(result$N_effective, 920)
})

test_that("assurance_mcse and assurance_lower are present and internally consistent", {
  # decision = "lower_bound": assurance_mcse/assurance_lower are the raw
  # single-N Monte Carlo diagnostics this test checks the formula for; they
  # are NA under decision = "isotonic" by design (see @details), so this
  # test must not rely on the default.
  #
  # assurance_mcse remains the plain (Wald) standard error,
  # sqrt(p(1-p)/B) -- a descriptive quantity, unaffected by defect 3's fix.
  # assurance_lower is now a WILSON score lower bound, not
  # joint_assurance - z * assurance_mcse (see @details, "The lower
  # confidence bound is Wilson, not Wald"): the two formulas agree closely
  # at this B but are not algebraically identical, which is exactly the
  # point of the fix (the Wald version collapses to a vacuous 1.000 at the
  # extremes regardless of B; see test-b1-vacuous-bound.R).
  result <- suppressWarnings(ss_unified(
    decision = "lower_bound",
    B = 1200, seed = 2026, prior_prev = c(4, 16),
    N_range = seq(200, 1600, by = 20), delta_auc = 0, check_nb = FALSE
  ))
  expect_false(is.null(result$assurance_mcse))
  expect_false(is.null(result$assurance_lower))
  expect_equal(
    result$assurance_mcse,
    sqrt(result$joint_assurance * (1 - result$joint_assurance) / result$B)
  )

  z <- stats::qnorm(0.95)
  phat <- result$joint_assurance
  n <- result$B
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half_width <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  expect_equal(result$assurance_lower, max(0, center - half_width))

  # The Wilson and Wald bounds must be close (not identical) at this B.
  wald_lower <- result$joint_assurance - z * result$assurance_mcse
  expect_true(abs(result$assurance_lower - wald_lower) < 0.01)
  expect_false(isTRUE(all.equal(result$assurance_lower, wald_lower)))
})

test_that("decision = \"lower_bound\" never accepts an N below its own lower bound", {
  # By construction, decision = "lower_bound" only stops the search once
  # wilson_lower(joint_assurance, B, qnorm(0.95)) >= target_assurance, so
  # assurance_lower at the returned N must never fall short of the target
  # (up to floating-point slack). Explicit decision: assurance_lower is NA
  # under the default decision = "isotonic" (see @details), so this
  # property is specific to "lower_bound" and must not rely on the default.
  for (s in c(1, 42, 2026)) {
    result <- suppressWarnings(ss_unified(
      decision = "lower_bound",
      B = 800, seed = s, prior_prev = c(4, 16),
      N_range = seq(200, 1600, by = 40), delta_auc = 0, check_nb = FALSE
    ))
    expect_gte(result$assurance_lower, result$target_assurance - 1e-8)
  }
})

test_that("decision = \"lower_bound\" never selects an EARLIER N than \"point\", for the same draws", {
  # lower_bound's stopping condition (wilson_lower(assurance, B, z) >=
  # target) implies point's (assurance >= target), because the Wilson
  # lower bound is never above its own point estimate (a valid confidence
  # interval always contains the estimate it is centred on -- see
  # wilson_lower()'s tests in test-helpers.R). So for identical random
  # draws (same seed, same B, same N_range) the first N accepted under
  # lower_bound can never come before the first N accepted under point.
  # This is a structural guarantee, not a statistical tendency, so it
  # holds even at the small B used here for speed.
  for (s in c(1, 7, 42, 99)) {
    common <- list(B = 500, seed = s, prior_prev = c(4, 16),
                   N_range = seq(200, 1600, by = 40),
                   delta_auc = 0, check_nb = FALSE)
    r_point <- suppressWarnings(
      do.call(ss_unified, c(common, list(decision = "point")))
    )
    r_lower <- suppressWarnings(
      do.call(ss_unified, c(common, list(decision = "lower_bound")))
    )
    expect_gte(r_lower$N_effective, r_point$N_effective)
  }
})

test_that("a warning suggests increasing B when the selected N is close to its own MC noise", {
  # This is distinct from warn_small_B(): it fires based on the margin
  # between the achieved assurance and target_assurance relative to the
  # MCSE at the SELECTED N, not on B alone, and its message is
  # distinguishable from the small-B warning's text. Not meaningful under
  # decision = "isotonic" (assurance_mcse is NA there, see @details), so
  # this test pins decision = "lower_bound" explicitly rather than relying
  # on the (now isotonic) default.
  # seed/N_range re-verified against the 0.5.1 independent-stream RNG (see
  # @details): the pre-0.5.1 scenario (seed = 2026, N_range up to 2200) no
  # longer lands close to the boundary under the new streams, so a scenario
  # that does was re-selected rather than forcing the old one to fit.
  expect_warning(
    result <- ss_unified(decision = "lower_bound",
                         B = 500, seed = 9, prior_sp = c(18, 2),
                         N_range = seq(400, 3000, by = 50),
                         delta_auc = 0, check_nb = TRUE,
                         nb_B_ceiling = 50000),
    "Monte Carlo noise"
  )
  expect_false(is.null(result$assurance_mcse))
})

# --- nb_ceiling: the net-benefit criterion's achievable ceiling --------

test_that("nb_ceiling is NA when check_nb = FALSE and numeric when check_nb = TRUE", {
  # N_range = seq(300, 500, by = 100) is deliberately small (this test only
  # cares about nb_ceiling, not about convergence), so the generic
  # "did not converge" warning is expected and suppressed here.
  no_nb <- suppressWarnings(ss_unified(B = 200, seed = 2026,
                      N_range = seq(300, 500, by = 100),
                      delta_auc = 0, check_nb = FALSE))
  expect_true("nb_ceiling" %in% names(no_nb))
  expect_true(is.na(no_nb$nb_ceiling))

  # N_range = c(300, 400) is far too small to reach target_assurance with
  # check_nb = TRUE, so the generic "did not converge" warning is expected
  # here and is not what this test is about; nb_ceiling is computed before
  # the search regardless of whether the search itself converges.
  # nb_B_ceiling is set small here: this test only checks presence/range,
  # not precision (see "nb_ceiling reproduces..." below for the precision
  # check, which needs the full default B_ceiling).
  with_nb <- suppressWarnings(ss_unified(
    B = 100, seed = 2026, prior_sp = c(18, 2), prior_prev = c(4, 16),
    Se_ref = 0.90, Sp_ref = 0.95, N_range = c(300, 400),
    delta_auc = 0, check_nb = TRUE, nb_B_ceiling = 50000
  ))
  expect_false(is.na(with_nb$nb_ceiling))
  expect_true(with_nb$nb_ceiling > 0 && with_nb$nb_ceiling < 1)
})

test_that("nb_ceiling reproduces the independently-verified value for the default pt_range (defect 4)", {
  # 0.879606 was derived INDEPENDENTLY of nb_assurance_ceiling() entirely,
  # by tensor-product Gauss-Legendre quadrature over the (prev, Se, Sp)
  # priors (n = 300 and n = 600 nodes per dimension agree to six decimals),
  # cross-checked against a 5e6-draw Monte Carlo run. This uses the
  # function's real default nb_B_ceiling (200000), which -- since 0.6.3 --
  # is a deterministic quadrature evaluation budget, not a Monte Carlo
  # sample size (see ?dtasamplesize:::nb_assurance_ceiling): versions
  # <= 0.6.0 used a Monte Carlo B_ceiling = 200000 and got 0.87864 here,
  # off by ~0.00097 (more than the tolerance below allows), and versions
  # 0.6.1-0.6.2 fixed precision by raising Monte Carlo B_ceiling to 2e8,
  # which reproduced this value but exhausted memory on ordinary hardware.
  result <- suppressWarnings(ss_unified(
    B = 100, seed = 2026, prior_se = c(17, 3), prior_sp = c(18, 2),
    prior_prev = c(4, 16), Se_ref = 0.90, Sp_ref = 0.95,
    N_range = c(300, 400), delta_auc = 0, check_nb = TRUE
  ))
  expect_equal(result$nb_ceiling, 0.879606, tolerance = 2e-4)
})

test_that("nb_assurance_ceiling() is a deterministic quadrature: identical across seed and RNGkind", {
  # Versions <= 0.6.2 estimated the ceiling by Monte Carlo, so its value
  # depended on `seed` (and, before an earlier fix, even on the caller's
  # ambient RNGkind -- see test-rng_state.R for that regression test).
  # Since 0.6.3 the ceiling is computed by closed-form integration over
  # prevalence and Gauss-Legendre quadrature over Se and Sp: it touches no
  # random number generator at all, so `seed` -- kept only for backward
  # compatibility -- can no longer move the result even a single bit.
  nb_assurance_ceiling <- get("nb_assurance_ceiling", envir = asNamespace("dtasamplesize"))
  prior_se <- c(17, 3); prior_sp <- c(18, 2); prior_prev <- c(4, 16)
  Se_ref <- 0.90; Sp_ref <- 0.95; pt_range <- c(0.15, 0.40)

  val_seed_a <- nb_assurance_ceiling(prior_se, prior_sp, prior_prev, Se_ref, Sp_ref,
                                      pt_range, seed = 1, B_ceiling = 200000L)
  val_seed_b <- nb_assurance_ceiling(prior_se, prior_sp, prior_prev, Se_ref, Sp_ref,
                                      pt_range, seed = 999999, B_ceiling = 200000L)
  expect_identical(val_seed_a, val_seed_b)
  expect_equal(val_seed_a, 0.879606, tolerance = 1e-5)

  # A second pt_range, called immediately after the first with the SAME
  # seed, must NOT reuse whatever the first call computed (the historical
  # bug this test used to target, back when the underlying draws were
  # random): it is a different, independently-checked ceiling.
  val_other_range <- nb_assurance_ceiling(prior_se, prior_sp, prior_prev, Se_ref, Sp_ref,
                                           c(0.20, 0.30), seed = 1, B_ceiling = 200000L)
  expect_false(isTRUE(all.equal(val_seed_a, val_other_range)))
  expect_equal(val_other_range, 0.968998, tolerance = 1e-5)
})

test_that("nb_B_ceiling above the maximum budget errors with a clear message instead of running", {
  # See NB_B_CEILING_MAX (R/helpers.R) and A1 in the audit this guards
  # against: earlier versions' Monte Carlo default (2e8 draws) allocated
  # several ~4.8 GB double vectors and could exhaust memory outright. The
  # cap below is checked before any computation is attempted.
  expect_error(
    ss_unified(B = 10, seed = 1, N_range = 300, delta_auc = 0,
               check_nb = TRUE, nb_B_ceiling = 50000000),
    "exceeds the maximum allowed budget"
  )
})

test_that("an unreachable check_nb ceiling REPLACES, not joins, the generic 'expand N_range' warning", {
  # pt_range = c(0.05, 0.50) drives the achievable ceiling to EXACTLY 0 for
  # these priors/Se_ref/Sp_ref -- a structural fact (at pt = 0.05, w =
  # 0.05/0.95 = (1 - Sp_ref)/Sp_ref, which makes the treat-all inequality's
  # Sp coefficient exactly 0, so it can never hold), not a Monte Carlo
  # estimate, so it is exact at any B_ceiling. target_assurance = 0.80 can
  # therefore never be reached at ANY N. This is a regression test for the
  # pre-search ceiling check: it FAILS against the pre-0.5.1 code, which
  # has no $nb_ceiling element at all and instead runs the full (here,
  # futile) grid search before emitting the generic "No N in N_range
  # achieved target assurance. Consider expanding N_range." message --
  # exactly the misleading suggestion this check replaces, since no
  # N_range could ever help here. nb_B_ceiling is kept small: the ceiling
  # here is exact, not approximate, so precision is not what this test
  # checks.
  w <- tryCatch({
    ss_unified(B = 100, seed = 2026, prior_sp = c(18, 2),
               prior_prev = c(4, 16), Se_ref = 0.90, Sp_ref = 0.95,
               pt_range = c(0.05, 0.50), N_range = c(300, 400),
               delta_auc = 0, check_nb = TRUE, nb_B_ceiling = 50000)
    NULL
  }, warning = function(w) w)

  expect_false(is.null(w))
  msg <- conditionMessage(w)
  expect_match(msg, "cannot reach target_assurance")
  expect_match(msg, "ceiling")
  expect_false(grepl("Consider expanding N_range", msg, fixed = TRUE))

  result <- suppressWarnings(
    ss_unified(B = 100, seed = 2026, prior_sp = c(18, 2),
               prior_prev = c(4, 16), Se_ref = 0.90, Sp_ref = 0.95,
               pt_range = c(0.05, 0.50), N_range = c(300, 400),
               delta_auc = 0, check_nb = TRUE, nb_B_ceiling = 50000)
  )
  expect_false(is.null(result$nb_ceiling))
  expect_equal(result$nb_ceiling, 0)
  expect_true(is.na(result$joint_assurance))
  # A2 fix: N_effective/n_total must be NA here too, not max(N_range) --
  # otherwise a caller reading only n_total would adopt a "sample size"
  # with a target that is structurally unreachable at any N.
  expect_true(is.na(result$N_effective))
  expect_true(is.na(result$n_total))
  expect_identical(result$status, "unreachable")
})

# --- decision = "isotonic" -----------------------------------------------

test_that("decision = \"isotonic\" is accepted and is the formal default", {
  # match.arg() would reject decision = "isotonic" against the pre-0.5.1
  # signature (only "lower_bound"/"point" were valid), so this line alone
  # is a regression check against that code.
  result <- ss_unified(decision = "isotonic", B = 300, seed = 2026,
                       N_range = seq(300, 900, by = 100),
                       delta_auc = 0, check_nb = FALSE)
  expect_equal(result$decision, "isotonic")
})

test_that("decision = \"isotonic\" always evaluates the complete N_range, even with full_grid = FALSE", {
  n_range <- seq(300, 900, by = 100)
  result <- ss_unified(decision = "isotonic", B = 300, seed = 2026,
                       N_range = n_range, delta_auc = 0, check_nb = FALSE,
                       full_grid = FALSE)
  expect_equal(nrow(result$grid_results), length(n_range))
})

test_that("decision = \"isotonic\" reports NA assurance_mcse, a non-NA margin-adjusted assurance_lower, and a fitted joint_assurance at/above target when found", {
  # assurance_mcse has no meaning for a pooled, possibly off-grid estimate
  # and stays NA, but assurance_lower is no longer NA: it is now the
  # margin-adjusted curve isotonic actually inverts to choose N_effective
  # (see @details for the margin), so it is populated and, by
  # construction, both it and the unadjusted joint_assurance clear
  # target_assurance whenever N_effective is found.
  result <- suppressWarnings(ss_unified(
    decision = "isotonic", B = 500, seed = 2026, prior_prev = c(4, 16),
    N_range = seq(200, 1600, by = 20), delta_auc = 0, check_nb = FALSE
  ))
  expect_true(is.na(result$assurance_mcse))
  expect_false(is.na(result$assurance_lower))
  expect_false(is.na(result$N_effective))
  expect_gte(result$assurance_lower, result$target_assurance - 1e-8)
  expect_gte(result$joint_assurance, result$assurance_lower - 1e-8)
  expect_gte(result$joint_assurance, result$target_assurance - 1e-8)
})

test_that("decision = \"isotonic\" falls back to the generic warning when the fitted curve never reaches target_assurance", {
  result <- suppressWarnings(ss_unified(
    decision = "isotonic", B = 100, seed = 2026,
    N_range = c(100, 150), delta_auc = 0, check_nb = FALSE
  ))
  expect_warning(
    ss_unified(decision = "isotonic", B = 100, seed = 2026,
               N_range = c(100, 150), delta_auc = 0, check_nb = FALSE),
    "Consider expanding N_range"
  )
  # N_effective/n_total are NA when the search does not converge (A2 fix):
  # a caller who reads only n_total must not see a number that looks like a
  # validated design. status distinguishes this ("grid_exhausted") from the
  # check_nb ceiling case ("unreachable").
  expect_true(is.na(result$N_effective))
  expect_true(is.na(result$n_total))
  expect_identical(result$status, "grid_exhausted")
})

test_that("decision = \"isotonic\" is exactly reproducible for a fixed seed", {
  args <- list(decision = "isotonic", B = 300, seed = 4242,
               N_range = seq(300, 1200, by = 50), delta_auc = 0,
               check_nb = FALSE)
  r1 <- do.call(ss_unified, args)
  r2 <- do.call(ss_unified, args)
  expect_identical(r1$N_effective, r2$N_effective)
  expect_identical(r1$joint_assurance, r2$joint_assurance)
  expect_identical(r1$grid_results, r2$grid_results)
})

# --- regression: single-usable-grid-point crash under decision = "isotonic" ---

test_that("decision = \"isotonic\" does not error on a single-point N_range that reaches target", {
  # Before this fix, inverting the isotonic fit called stats::approx(x, y,
  # xout, rule = 2) unconditionally, which requires >= 2 x values. With a
  # single-N N_range (e.g. checking one specific, already-decided N, the
  # most basic use case of the function), grid_results/Nseq/fitted all had
  # length 1. If that one N did NOT reach target_assurance, the code never
  # reached the approx() call and was fine; if it DID reach it -- exactly
  # the case that matters, since that is what "verify this N works" means
  # -- approx() crashed with "need at least two non-NA values to
  # interpolate". N = 1600 here clears target_assurance = 0.80 with a wide
  # margin (joint_assurance ~= 0.97), so this reproduces the crash
  # regardless of the conservative margin added alongside this fix.
  expect_no_error(
    result <- suppressWarnings(ss_unified(
      B = 800, seed = 2026, prior_prev = c(4, 16),
      N_range = 1600:1600, delta_auc = 0, check_nb = FALSE, full_grid = TRUE
    ))
  )
  expect_equal(result$N_effective, 1600L)
  expect_equal(nrow(result$grid_results), 1)
  expect_false(is.na(result$joint_assurance))
  expect_false(is.na(result$assurance_lower))
})

test_that("decision = \"isotonic\" does not error on the reported single-N reproduction (N_range = 2216:2216)", {
  # The literal reproduction from the bug report: prior_prev = c(4, 16),
  # Se_ref = 0.90, Sp_ref = 0.95, check_nb = TRUE, target_assurance = 0.80,
  # a single N = 2216 whose raw assurance sits almost exactly on target
  # (independently verified: joint_assurance = 0.800112 at B = 1e6,
  # seed = 21). N_range = 2110:2110 never crashed pre-fix, because the
  # target was not reached there and the crashing approx() call was never
  # reached -- i.e. the bug fired exactly when the single N DID reach the
  # target, which is the case exercised here.
  expect_no_error(
    result <- suppressWarnings(ss_unified(
      prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
      Se_ref = 0.90, Sp_ref = 0.95, loss_rate = 0,
      delta_se = 0.07, delta_sp = 0.05, delta_auc = 0.06,
      check_nb = TRUE, target_assurance = 0.80,
      N_range = 2216:2216, B = 1e6, seed = 21, full_grid = TRUE,
      nb_B_ceiling = 50000
    ))
  )
  # The point estimate (0.800112) clears the target, but the conservative
  # Monte Carlo margin puts the lower bound below it, so the rule declines
  # this N rather than accepting a size whose true assurance sits on the
  # boundary. An independent reference curve places the truth at N = 2216
  # at 0.80000, so declining is the correct, conservative outcome.
  expect_equal(result$joint_assurance, 0.800112, tolerance = 1e-4)
  expect_lt(result$assurance_lower, 0.80)
  expect_true(is.na(result$N_effective))
  expect_identical(result$status, "grid_exhausted")
})

# --- regression: decision = "isotonic" must keep a safety margin ----------

test_that("decision = \"isotonic\"'s conservative margin guarantees assurance_lower/joint_assurance at/above target whenever N_effective is found", {
  # Before this fix, decision = "isotonic" inverted the smoothed curve
  # with NO safety margin: it accepted the first N whose SMOOTHED POINT
  # ESTIMATE crossed target_assurance, and a point estimate lands above
  # its own true value roughly half the time. Measured against a
  # B = 1e7-per-point reference curve on the step-5 scenario (see
  # @details), that left a substantial fraction of seeds selecting an N
  # whose TRUE assurance fell short of target_assurance = 0.80, despite
  # sd/bias/RMSE all far better than the pre-isotonic first-crossing rule.
  # The fix subtracts a conservative per-N Monte Carlo margin from the
  # fitted curve before inverting it, which makes it algebraically
  # guaranteed (not merely likely) that whenever N_effective is found,
  # both assurance_lower (the margin-adjusted curve actually inverted) and
  # joint_assurance (always >= assurance_lower, since the margin is >= 0
  # pointwise) are >= target_assurance. That guarantee is a structural
  # property of the search rule and can be checked cheaply here; the
  # calibration this buys back against ground truth (a high fraction of
  # seeds landing at/above the true target) is measured at larger scale in
  # @details rather than in this fast test.
  for (s in c(1, 7, 42, 99, 2026)) {
    result <- suppressWarnings(ss_unified(
      decision = "isotonic", B = 400, seed = s, prior_prev = c(4, 16),
      N_range = seq(200, 1600, by = 20), delta_auc = 0, check_nb = FALSE
    ))
    if (is.na(result$N_effective)) next  # did-not-converge fallback: N/A here
    expect_gte(result$assurance_lower, result$target_assurance - 1e-8)
    expect_gte(result$joint_assurance, result$assurance_lower - 1e-8)
  }
})

# --- regression: decision = "isotonic" must be invariant to N_range's order --

test_that("decision = \"isotonic\" gives identical results for ascending, descending and shuffled N_range", {
  # stats::isoreg() returns $x in the ORIGINAL (input) order and $yf in the
  # SORTED order -- pairing them positionally, as the pre-fix code did,
  # silently mismatches each N against a fitted value belonging to a
  # DIFFERENT N whenever N_range is not itself already increasing, with no
  # warning of any kind. The RNG stream assigned to each N was also keyed
  # by N's POSITION in N_range rather than by its value, so reordering
  # N_range legitimately handed a different candidate N to a different
  # random stream on top of the isoreg mismatch. Independently verified
  # against the pre-fix code for this exact scenario (same N values, same
  # B, same seed, only the order of N_range changes): ascending selected
  # N_effective = 729 (joint_assurance = 0.8104), descending selected
  # N_effective = 677 (joint_assurance = 0.8101), and a shuffled order
  # selected N_effective = 415 with joint_assurance = 0.4156 and
  # assurance_lower = 0.4041 -- BELOW target_assurance = 0.80, and with no
  # warning at all. This test fails against that code and must pass here.
  base <- seq(200, 1200, by = 100)
  ascending <- base
  descending <- rev(base)
  shuffled <- c(600, 200, 1000, 400, 1200, 300, 800, 500, 700, 900, 1100)
  stopifnot(setequal(shuffled, base))  # same SET of N values, different order

  common <- list(B = 4000, delta_auc = 0, seed = 2026)
  r_asc <- suppressWarnings(do.call(ss_unified, c(list(N_range = ascending), common)))
  r_desc <- suppressWarnings(do.call(ss_unified, c(list(N_range = descending), common)))
  r_shuf <- suppressWarnings(do.call(ss_unified, c(list(N_range = shuffled), common)))

  expect_identical(r_asc$N_effective, r_desc$N_effective)
  expect_identical(r_asc$N_effective, r_shuf$N_effective)
  expect_identical(r_asc$joint_assurance, r_desc$joint_assurance)
  expect_identical(r_asc$joint_assurance, r_shuf$joint_assurance)
  expect_identical(r_asc$assurance_lower, r_desc$assurance_lower)
  expect_identical(r_asc$assurance_lower, r_shuf$assurance_lower)

  # The documented guarantee (assurance_lower >= target_assurance whenever
  # N_effective is found) must hold regardless of N_range's order -- it is
  # exactly what the pre-fix shuffled case violated (0.4041 < 0.80).
  expect_gte(r_shuf$assurance_lower, r_shuf$target_assurance - 1e-8)
})

test_that("decision %in% c(\"point\", \"lower_bound\") select the SAME (smallest) N_effective regardless of N_range's order (defect 1)", {
  # Versions <= 0.6.1 accepted whichever N reached target_assurance FIRST
  # IN N_range's OWN ORDER -- correct only when N_range happens to already
  # be ascending. Independently verified against that code for this exact
  # scenario (identical priors, B, seed and SET of 17 candidate N
  # throughout): decision = "point" with full_grid = TRUE selected
  # N_effective = 850 (the true minimum) given N_range ascending, but 1200
  # (+41%) given descending, and a third value given shuffled -- with no
  # warning of any kind.
  G <- seq(400, 1200, by = 50)
  descending <- rev(G)
  set.seed(99)
  shuffled <- sample(G)
  stopifnot(setequal(shuffled, G))

  common <- list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
                 Se_ref = 1, Sp_ref = 1, delta_se = 0.07, delta_sp = 0.05,
                 delta_auc = 0, check_nb = FALSE, target_assurance = 0.80,
                 B = 1000, seed = 2026, full_grid = TRUE)

  for (dec in c("point", "lower_bound")) {
    r_asc  <- suppressWarnings(do.call(ss_unified,
      c(list(N_range = G, decision = dec), common)))
    r_desc <- suppressWarnings(do.call(ss_unified,
      c(list(N_range = descending, decision = dec), common)))
    r_shuf <- suppressWarnings(do.call(ss_unified,
      c(list(N_range = shuffled, decision = dec), common)))

    expect_identical(r_asc$N_effective, r_desc$N_effective)
    expect_identical(r_asc$N_effective, r_shuf$N_effective)
    expect_identical(r_asc$joint_assurance, r_desc$joint_assurance)
    expect_identical(r_asc$joint_assurance, r_shuf$joint_assurance)
  }

  # The specific, independently verified reproduction for decision = "point":
  # ascending selects the true minimum N in G that reaches target_assurance.
  r_point_asc <- suppressWarnings(do.call(ss_unified,
    c(list(N_range = G, decision = "point"), common)))
  expect_equal(r_point_asc$N_effective, 850)
})

test_that("full_grid = FALSE still finds the smallest N regardless of N_range's order, at lower cost (defect 1)", {
  # With full_grid = FALSE the search stops as soon as it finds an
  # accepting N; that N must still be the SMALLEST one in N_range that
  # reaches target_assurance, not merely the first one encountered in
  # N_range's own order (which, before the fix, could stop the search
  # early at an N far from the true minimum when N_range was not
  # ascending).
  G <- seq(400, 1200, by = 50)
  common <- list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
                 Se_ref = 1, Sp_ref = 1, delta_se = 0.07, delta_sp = 0.05,
                 delta_auc = 0, check_nb = FALSE, target_assurance = 0.80,
                 B = 1000, seed = 2026, decision = "point", full_grid = FALSE)

  r_asc  <- suppressWarnings(do.call(ss_unified, c(list(N_range = G), common)))
  r_desc <- suppressWarnings(do.call(ss_unified, c(list(N_range = rev(G)), common)))

  expect_identical(r_asc$N_effective, r_desc$N_effective)
  expect_equal(r_asc$N_effective, 850)
})

test_that("the \"did not converge\" fallback reports the assurance actually observed AT the returned N, regardless of N_range's order (defect 2)", {
  # Versions <= 0.6.1 reported joint_assurance/assurance_mcse/assurance_lower
  # left over from whichever N the search loop happened to evaluate LAST,
  # not from the N being returned (max(N_range)). Independently verified
  # against that code for this exact non-converging search (identical
  # priors, B, seed and SET of candidate N throughout, decision =
  # "isotonic"): the SAME returned N_effective = 400 was reported with
  # joint_assurance = 0.37850 given N_range ascending (correct), 0.06400
  # given descending (the assurance at N = 200, not N = 400), and 0.30000
  # given shuffled -- three different numbers attached to the identical N.
  low <- seq(200, 400, by = 20)
  run <- function(g) suppressWarnings(ss_unified(
    prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
    Se_ref = 1, Sp_ref = 1, delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
    check_nb = FALSE, target_assurance = 0.80, N_range = g, B = 2000,
    seed = 2026, decision = "isotonic"
  ))

  r_asc <- run(low)
  r_desc <- run(rev(low))
  set.seed(7)
  r_shuf <- run(sample(low))

  # N_effective is NA when the search does not converge (A2 fix), regardless
  # of N_range's order; the assurance actually observed at max(N_range) --
  # this test's own subject -- remains available via joint_assurance et al.
  # (checked below), unaffected by that change.
  expect_true(is.na(r_asc$N_effective))
  expect_true(is.na(r_desc$N_effective))
  expect_true(is.na(r_shuf$N_effective))
  expect_identical(r_asc$status, "grid_exhausted")

  expect_equal(r_desc$joint_assurance, r_asc$joint_assurance)
  expect_equal(r_shuf$joint_assurance, r_asc$joint_assurance)
  expect_equal(r_desc$assurance_mcse, r_asc$assurance_mcse)
  expect_equal(r_shuf$assurance_lower, r_asc$assurance_lower)
})

test_that("N_range is validated on entry, with a diagnostic message identifying the problem (defect 4)", {
  # Versions <= 0.6.1 did not validate N_range at all: a negative, NA,
  # infinite or non-integer element reached rbinom()/rmultinom() deep
  # inside the search loop and died with a generic, uninformative
  # "missing value where TRUE/FALSE needed" -- no indication of which
  # argument, or which element of it, was at fault.
  base_args <- list(B = 100, seed = 2026, delta_auc = 0, check_nb = FALSE)

  expect_error(
    do.call(ss_unified, c(list(N_range = c(-100, 200, 400)), base_args)),
    "below 1"
  )
  expect_error(
    do.call(ss_unified, c(list(N_range = c(200, NA, 400)), base_args)),
    "NA"
  )
  expect_error(
    do.call(ss_unified, c(list(N_range = c(200, Inf, 400)), base_args)),
    "non-finite"
  )
  expect_error(
    do.call(ss_unified, c(list(N_range = c(200.5, 300.2, 400)), base_args)),
    "non-integer"
  )
  expect_error(
    do.call(ss_unified, c(list(N_range = numeric(0)), base_args)),
    "length"
  )
  # None of these should reach the generic, undiagnostic internal R error
  # this defect used to produce.
  for (bad in list(c(-100, 200, 400), c(200, NA, 400), c(200, Inf, 400),
                    c(200.5, 300.2, 400))) {
    err <- tryCatch(
      do.call(ss_unified, c(list(N_range = bad), base_args)),
      error = function(e) conditionMessage(e)
    )
    expect_false(grepl("missing value where TRUE/FALSE needed", err, fixed = TRUE))
  }

  # A degenerate value such as 0 -- previously accepted, silently shifting
  # N_effective via the same block-pooling mechanism as defect 3 -- is now
  # rejected outright.
  expect_error(
    do.call(ss_unified, c(list(N_range = c(0, 200, 400)), base_args)),
    "below 1"
  )
})

test_that("decision = \"isotonic\" tolerates a duplicated N_range, without error or a stray warning, and stays order-invariant", {
  # A duplicated N must keep working (not crash, not warn) and must not
  # break the order-invariance above: a repeated N always looks up the
  # same value-keyed RNG stream, so it deterministically contributes the
  # same simulated result each time it appears, wherever in N_range it sits.
  base <- seq(200, 1200, by = 100)
  dup <- c(base, 700, 700, 200)
  dup_shuffled <- c(200, 700, 700, 700, 600, 400, 300, 1200, 900, 1100, 500,
                     1000, 800, 200)
  stopifnot(setequal(dup, dup_shuffled), length(dup) == length(dup_shuffled))

  common <- list(B = 2000, delta_auc = 0, seed = 2026)
  expect_no_warning(
    r1 <- do.call(ss_unified, c(list(N_range = dup), common))
  )
  expect_no_warning(
    r2 <- do.call(ss_unified, c(list(N_range = dup_shuffled), common))
  )

  expect_identical(r1$N_effective, r2$N_effective)
  expect_identical(r1$joint_assurance, r2$joint_assurance)
  expect_identical(r1$assurance_lower, r2$assurance_lower)
  expect_equal(nrow(r1$grid_results), length(dup))

  # --- defect 3: duplicated grid points must NOT fabricate precision -----
  # `dup`/`dup_shuffled` above only ever compare duplicated grids against
  # PERMUTATIONS OF THEMSELVES (same multiset of N, different order), which
  # cannot detect a defect that inflates the effective sample size behind
  # the margin -- that defect biases every permutation of `dup` identically,
  # since it depends on how many times each N repeats, not on the order.
  # The result for a duplicated grid must instead match the result for the
  # DEDUPLICATED base grid exactly: a repeated N carries no new information
  # (it always replays the identical stream), so it must not change
  # N_effective, joint_assurance, or assurance_lower at all. Versions
  # <= 0.6.2-pre (with the isoreg order fix already in place, but before
  # the pre-fit dedup) instead let each repeat inflate the pooled block's
  # effective sample size, so this comparison -- unlike the
  # permutation-only checks above -- fails against that code.
  r_base <- suppressWarnings(do.call(ss_unified, c(list(N_range = base), common)))
  expect_identical(r1$N_effective, r_base$N_effective)
  expect_identical(r1$joint_assurance, r_base$joint_assurance)
  expect_identical(r1$assurance_lower, r_base$assurance_lower)
})

test_that("decision = \"isotonic\" is unaffected by HOW MANY TIMES each N is duplicated (defect 3)", {
  # Versions <= 0.6.2-pre: N_effective fell from 857 (each N evaluated
  # once) to 835 (each N repeated 10x) as the replication count rose, with
  # assurance_lower wrongly INCREASING (0.800163 to 0.800275) -- fabricated
  # precision from re-counting identical, non-independent evidence.
  # Deduplicating before the isotonic fit removes the mechanism entirely:
  # N_effective and assurance_lower must be identical across every
  # replication count.
  g <- seq(700, 1000, by = 20)
  run_reps <- function(reps) {
    suppressWarnings(ss_unified(
      prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
      Se_ref = 1, Sp_ref = 1, delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
      check_nb = FALSE, target_assurance = 0.80,
      N_range = rep(g, each = reps), B = 4000, seed = 2026,
      decision = "isotonic"
    ))
  }
  results <- lapply(c(1, 2, 3, 5, 10), run_reps)
  n_eff <- vapply(results, function(r) r$N_effective, integer(1))
  al <- vapply(results, function(r) r$assurance_lower, numeric(1))

  expect_true(all(n_eff == n_eff[1]))
  expect_equal(al, rep(al[1], length(al)))
})
