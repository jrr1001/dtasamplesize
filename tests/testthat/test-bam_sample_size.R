
# NOTE ON `method`: bam_sample_size() defaults to method = "exact" (v0.5.0),
# a deterministic closed-form calculation for the headline joint search
# (N_total / joint_assurance / assurance_mcse). The tests below that were
# written against the older simulation-based joint search pass
# method = "monte_carlo" explicitly, so they keep exercising exactly the
# behavior (including the anti-noise acceptance rule and the degenerate-
# replication handling) they were designed and tuned for. Tests specific to
# method = "exact" are grouped further down.

test_that("BAM n_diseased is near Buderer with informative prior", {
  # n_range must reach n_sp (= 373 here), otherwise the Sp search falls back to
  # max(n_range) with a warning and this test would pass on that branch.
  # N_range is fixed to a small explicit grid so this test only exercises the
  # (cheap) per-arm searches, not the (more expensive) joint N search.
  expect_no_warning(
    result <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14,
                              B = 5000, seed = 2026, n_range = 20:600,
                              N_range = seq(300, 700, by = 20),
                              method = "monte_carlo")
  )
  # With prior Beta(17,3) and delta_se=0.14, BAM gives ~107 diseased
  expect_true(result$n_diseased > 90 & result$n_diseased < 130)
  # Converged strictly inside the range (not pinned to the ceiling)
  expect_lt(result$n_diseased, 600L)
  expect_lt(result$n_non_diseased, 600L)
})

test_that("BAM assurance meets target", {
  expect_no_warning(
    result <- bam_sample_size(B = 5000, seed = 2026, n_range = 20:600,
                              N_range = seq(300, 700, by = 20),
                              method = "monte_carlo")
  )
  expect_true(result$assurance_se >= 0.80)
  expect_true(result$assurance_sp >= 0.80)
})

test_that("BAM returns dtasamplesize class with expected fields", {
  expect_no_warning(
    result <- bam_sample_size(B = 1000, seed = 42, n_range = 20:600,
                              N_range = seq(300, 700, by = 20),
                              method = "monte_carlo")
  )
  expect_s3_class(result, "dtasamplesize")
  expect_true(!is.null(result$n_diseased))
  expect_true(!is.null(result$n_non_diseased))
  expect_true(!is.null(result$N_total_median))
  expect_true(!is.null(result$assurance_se))
  expect_true(!is.null(result$assurance_sp))
  # New in v0.5.0: the joint-assurance headline result and its diagnostics.
  expect_true(!is.null(result$N_total))
  expect_true(!is.null(result$joint_assurance))
  expect_true(!is.null(result$assurance_mcse))
  expect_identical(result$n_total, result$N_total)
  # New in the exact-mode wave: the method actually used, kept separate from
  # $method (which holds this object's human-readable method NAME, shared
  # across the package's print.dtasamplesize()).
  expect_identical(result$assurance_method, "monte_carlo")
  expect_identical(result$method, "Bayesian Assurance Method (BAM) for DTA Sample Size")
})

test_that("BAM vague prior requires more subjects than informative", {
  # The vague prior_se = c(2, 2) needs a much larger n_se (~188) than the
  # informative default (~107), so the joint search also needs a wider
  # N_range to converge without falling back (and warning).
  rng <- seq(300, 900, by = 10)
  expect_no_warning(
    res_informative <- bam_sample_size(prior_se = c(17, 3), delta_se = 0.14,
                                       B = 3000, seed = 2026, n_range = 20:600,
                                       N_range = rng, method = "monte_carlo")
  )
  expect_no_warning(
    res_vague <- bam_sample_size(prior_se = c(2, 2), delta_se = 0.14,
                                 B = 3000, seed = 2026, n_range = 20:600,
                                 N_range = rng, method = "monte_carlo")
  )
  expect_true(res_vague$n_diseased > res_informative$n_diseased)
})

test_that("BAM N_total covers BOTH arms (Se and Sp requirements)", {
  result <- bam_sample_size(B = 2000, seed = 2026, n_range = 20:600,
                            N_range = seq(300, 700, by = 20),
                            method = "monte_carlo")
  # The joint N_total must, on average, supply n_se diseased AND n_sp
  # non-diseased subjects; at the default prior_prev (E[prev] = 0.30) it can
  # never rationally fall below the naive geometric floor for the binding
  # (here, specificity) arm.
  expect_gte(result$N_total,
             ceiling(result$n_non_diseased / (1 - 0.30)))
  expect_true(is.finite(result$N_total_P90))
})

test_that("BAM still warns when n_range is genuinely too short", {
  # The fallback branch must keep announcing itself. n_range = 20:150 lets the
  # Se search converge (n_se ~ 107) but not the Sp one (n_sp ~ 373), so the
  # per-arm fallback warning fires. Because n_sp is then capped at 150 (far
  # below the true ~373 requirement), the joint total-N search also cannot
  # reach target_assurance within N_range = seq(150, 250, by = 10), so a
  # SECOND, independently legitimate warning fires too: a too-short n_range
  # genuinely cannot deliver a valid joint result. Both are expected.
  msgs <- character(0)
  res <- withCallingHandlers(
    bam_sample_size(B = 500, seed = 2026, n_range = 20:150,
                    N_range = seq(150, 250, by = 10), method = "monte_carlo"),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("achieved target assurance for Sp", msgs)))
  expect_true(any(grepl(
    "No N in N_range achieved the target JOINT assurance", msgs
  )))
  expect_equal(res$n_non_diseased, 150L)
  expect_lt(res$n_diseased, 150L)   # the Se arm did converge
})

test_that("BAM joint search accepts N_total only once the JOINT assurance genuinely meets target", {
  # A scenario deliberately unlike the package defaults, chosen to converge
  # inside a small, explicit N_range so the test stays fast.
  result <- bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(17, 3),
    delta_se = 0.14, delta_sp = 0.14,
    target_assurance = 0.80,
    prior_prev = c(10, 10),
    n_range = 20:200,
    N_range = seq(50, 400, by = 5),
    B = 3000, seed = 2026, method = "monte_carlo"
  )
  expect_true(is.finite(result$N_total))
  expect_true(result$N_total < max(seq(50, 400, by = 5)))  # genuine convergence, not a ceiling fallback
  # The point estimate must reach target_assurance ...
  expect_gte(result$joint_assurance, result$target_assurance)
  # ... AND survive the anti-noise rule: the lower bound of the one-sided 95%
  # Monte Carlo confidence interval must also reach target_assurance.
  expect_gte(
    result$joint_assurance - qnorm(0.95) * result$assurance_mcse,
    result$target_assurance
  )
})

test_that("BAM N_total grows when precision is tightened (smaller delta_se)", {
  common <- list(
    prior_se = c(17, 3), prior_sp = c(17, 3),
    target_assurance = 0.80, prior_prev = c(10, 10),
    B = 3000, seed = 2026, method = "monte_carlo"
  )
  res_loose <- do.call(bam_sample_size, c(common, list(
    delta_se = 0.14, delta_sp = 0.14,
    n_range = 20:200, N_range = seq(50, 400, by = 5)
  )))
  res_tight <- do.call(bam_sample_size, c(common, list(
    delta_se = 0.08, delta_sp = 0.14,
    n_range = 20:400, N_range = seq(50, 900, by = 5)
  )))
  expect_gt(res_tight$N_total, res_loose$N_total)
})

test_that("BAM N_total and joint_assurance are reproducible with a fixed seed", {
  args <- list(
    prior_se = c(17, 3), prior_sp = c(17, 3),
    delta_se = 0.14, delta_sp = 0.14, target_assurance = 0.80,
    prior_prev = c(10, 10), n_range = 20:200,
    N_range = seq(50, 400, by = 5), B = 1000, seed = 99,
    method = "monte_carlo"
  )
  res_a <- do.call(bam_sample_size, args)
  res_b <- do.call(bam_sample_size, args)
  expect_identical(res_a$N_total, res_b$N_total)
  expect_identical(res_a$joint_assurance, res_b$joint_assurance)
  expect_identical(res_a$assurance_mcse, res_b$assurance_mcse)
})

test_that("BAM treats a degenerate replication (n_d = 0 or n_nd = 0) as a failure, not a drop", {
  # At N = 1, EVERY replication is degenerate: n_d ~ Binomial(1, prev) is
  # either 0 (n_nd = 1) or 1 (n_nd = 0), so one arm always has zero subjects.
  # Priors are made deliberately very tight (far tighter than delta_se /
  # delta_sp), so that a posterior computed from the prior ALONE (i.e., with
  # no data, exactly what happens on the degenerate arm) would numerically
  # satisfy the width target on its own. If degenerate replications were
  # merely left to "fail naturally" rather than being explicitly forced to
  # fail, joint_assurance at N = 1 would come out near 1 here. The explicit
  # override must instead force it to exactly 0. method = "exact" (see
  # further down) needs an analogous override: it enumerates n_d = 0 and
  # n_d = N exhaustively, but without forcing the arm-size-0 term to 0 it
  # would credit those enumerated terms using the prior-only credible
  # interval width, which is not generally 0 (see the regression test
  # below).
  expect_warning(
    result <- bam_sample_size(
      prior_se = c(1000, 100), prior_sp = c(1000, 100),
      delta_se = 0.14, delta_sp = 0.14,
      target_assurance = 0.80,
      prior_prev = c(10, 10),
      n_range = 5:10,
      N_range = 1,
      B = 2000, seed = 2026, method = "monte_carlo"
    ),
    "No N in N_range achieved the target JOINT assurance"
  )
  expect_identical(result$joint_assurance, 0)
  expect_identical(result$N_total, 1L)
})

## ---------------------------------------------------------------------
## method = "exact": closed-form joint search (default since v0.5.0)
## ---------------------------------------------------------------------

test_that("BAM exact mode reproduces independently verified reference values", {
  # This is the headline scenario carried into the paper: with only
  # prior_prev and target_assurance changed from the package defaults, the
  # true (exact, not simulated) first N with joint assurance >= 0.80 is 678,
  # one step above N = 677 (0.7996849, just short of target). All three
  # values below were verified independently via an exact beta-binomial sum.
  expect_no_warning(
    res_678 <- bam_sample_size(
      prior_prev = c(4, 16), target_assurance = 0.80, method = "exact"
    )
  )
  expect_identical(res_678$N_total, 678L)
  expect_equal(res_678$joint_assurance, 0.800349, tolerance = 1e-6)
  expect_identical(res_678$assurance_mcse, 0)
  expect_identical(res_678$assurance_method, "exact")

  # N = 591 and N = 677 both fall short of target_assurance, so the search
  # (restricted to that single candidate) falls back to max(N_range) with a
  # warning; only the resulting joint_assurance VALUE is of interest here.
  # n_range is deliberately tiny: N_range is explicit, so n_se / n_sp never
  # feed into this exact search, and there is no need to pay for the (much
  # more expensive, and here irrelevant) default per-arm search.
  res_591 <- suppressWarnings(bam_sample_size(
    prior_prev = c(4, 16), target_assurance = 0.80, method = "exact",
    n_range = 20:20, N_range = 591
  ))
  expect_equal(res_591$joint_assurance, 0.7244942, tolerance = 1e-6)

  res_677 <- suppressWarnings(bam_sample_size(
    prior_prev = c(4, 16), target_assurance = 0.80, method = "exact",
    n_range = 20:20, N_range = 677
  ))
  expect_equal(res_677$joint_assurance, 0.7996849, tolerance = 1e-6)
})

test_that("BAM exact mode defaults to method = \"exact\" when method is not passed", {
  # N_range must actually reach the true exact crossing point (275, per the
  # "agree within Monte Carlo error" test below) so this converges cleanly
  # with no fallback warning; n_range is tiny on purpose (N_range is
  # explicit, so it cannot affect this exact search either way) but must
  # still be wrapped in suppressWarnings() since it deliberately never
  # converges on its own.
  result <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(17, 3),
    delta_se = 0.14, delta_sp = 0.14, target_assurance = 0.80,
    prior_prev = c(10, 10), n_range = 20:22,
    N_range = seq(50, 400, by = 5)
    # method intentionally omitted
  ))
  expect_identical(result$assurance_method, "exact")
  expect_identical(result$assurance_mcse, 0)
  expect_true(is.finite(result$N_total))
})

test_that("BAM exact mode is deterministic: unaffected by seed or B", {
  # n_range is tiny on purpose (N_range is explicit, so the legacy per-arm
  # search's outcome cannot feed back into the exact joint search); this
  # keeps the test fast even at the large B used below.
  common <- list(
    prior_se = c(17, 3), prior_sp = c(17, 3),
    delta_se = 0.14, delta_sp = 0.14, target_assurance = 0.80,
    prior_prev = c(10, 10), n_range = 20:22,
    N_range = seq(50, 400, by = 5), method = "exact"
  )
  res_a <- suppressWarnings(do.call(bam_sample_size, c(common, list(seed = 1, B = 10))))
  res_b <- suppressWarnings(do.call(bam_sample_size, c(common, list(seed = 424242, B = 50000))))
  expect_identical(res_a$N_total, res_b$N_total)
  expect_identical(res_a$joint_assurance, res_b$joint_assurance)
  expect_identical(res_a$assurance_mcse, res_b$assurance_mcse)

  # Two calls with identical arguments (including seed/B) must be identical
  # bit-for-bit -- there is no simulation to reseed.
  res_c <- suppressWarnings(do.call(bam_sample_size, c(common, list(seed = 1, B = 10))))
  expect_identical(res_a, res_c)
})

test_that("BAM exact and monte_carlo methods agree within Monte Carlo error", {
  common <- list(
    prior_se = c(17, 3), prior_sp = c(17, 3),
    delta_se = 0.14, delta_sp = 0.14, target_assurance = 0.80,
    prior_prev = c(10, 10), n_range = 20:22,
    N_range = seq(50, 400, by = 5)
  )
  res_exact <- suppressWarnings(do.call(bam_sample_size, c(common, list(method = "exact"))))
  res_mc <- suppressWarnings(do.call(bam_sample_size, c(
    common, list(method = "monte_carlo", B = 20000, seed = 2026)
  )))

  # The Monte Carlo search (with its anti-noise rule) can land a few grid
  # steps above the exact crossing point, never below it and never far above.
  expect_gte(res_mc$N_total, res_exact$N_total)
  expect_lte(res_mc$N_total - res_exact$N_total, 5 * 5L)  # at most 5 grid steps (grid step = 5)

  # The real test of agreement: what the exact calculation says the TRUE
  # probability is at the Monte Carlo method's own chosen N should match the
  # Monte Carlo estimate there, within a small number of Monte Carlo standard
  # errors (assurance_mcse) -- not within an arbitrary absolute tolerance.
  exact_at_mc_N <- suppressWarnings(do.call(bam_sample_size, c(
    common[setdiff(names(common), "N_range")],
    list(method = "exact", N_range = res_mc$N_total)
  )))
  expect_lt(
    abs(exact_at_mc_N$joint_assurance - res_mc$joint_assurance),
    4 * res_mc$assurance_mcse
  )
})

test_that("BAM exact mode treats a degenerate replication (n_d = 0 or n_nd = 0) as a failure, matching Monte Carlo", {
  # Companion to "BAM treats a degenerate replication ... as a failure, not a
  # drop" above, but for method = "exact". Priors are deliberately tight
  # enough that the credible interval width computed from the PRIOR ALONE
  # (i.e., an arm that received zero subjects) already satisfies delta_se /
  # delta_sp on its own: P_se(0) and P_sp(0), the per-arm probabilities at
  # arm size 0, would both be 1 if the degenerate-replication rule were not
  # enforced. Under prior_prev = c(1, 19), a small N such as 10 puts most of
  # the prevalence mass near 0, so P(n_d = 0) alone is about 0.66 there:
  # unless degenerate splits are excluded, the exact joint assurance comes
  # out at exactly 1 regardless of N, which the assertions below rule out.
  common <- list(
    prior_se = c(200, 20), prior_sp = c(400, 20), prior_prev = c(1, 19),
    delta_se = 0.10, delta_sp = 0.06, target_assurance = 0.999,
    n_range = 20:20
  )

  for (N in c(10L, 20L, 30L)) {
    res_exact <- suppressWarnings(do.call(bam_sample_size, c(common, list(
      method = "exact", N_range = N
    ))))
    res_mc <- suppressWarnings(do.call(bam_sample_size, c(common, list(
      method = "monte_carlo", N_range = N, B = 100000, seed = 2026
    ))))

    # A joint assurance of 1 at these small N, under informative priors this
    # tight, is only possible if degenerate splits are being credited as
    # successes -- the exact defect this test guards against.
    expect_lt(res_exact$joint_assurance, 0.9)
    # The exact and Monte Carlo calculations must agree on the same
    # generative model, within a small number of Monte Carlo standard
    # errors (not merely both be "less than 1").
    expect_lt(
      abs(res_exact$joint_assurance - res_mc$joint_assurance),
      5 * res_mc$assurance_mcse
    )
  }
})

## ---------------------------------------------------------------------
## $assurance_method (computation mode, kept distinct from the descriptive
## $method name) and the small-B advisory scoping fix
## ---------------------------------------------------------------------

test_that("BAM $assurance_method records the computation mode actually used, in both modes, while $method keeps the descriptive name", {
  res_exact <- bam_sample_size(prior_prev = c(4, 16), target_assurance = 0.80)
  expect_identical(res_exact$assurance_method, "exact")
  expect_identical(
    res_exact$method,
    "Bayesian Assurance Method (BAM) for DTA Sample Size"
  )

  res_mc <- bam_sample_size(
    prior_prev = c(4, 16), target_assurance = 0.80, seed = 99, B = 2000,
    N_range = seq(600, 700, by = 5), method = "monte_carlo"
  )
  expect_identical(res_mc$assurance_method, "monte_carlo")
  expect_identical(
    res_mc$method,
    "Bayesian Assurance Method (BAM) for DTA Sample Size"
  )

  # A user who selects a mode via `method` and later inspects the result
  # must recover their own choice from $assurance_method, not the shared
  # descriptive $method string.
  expect_false(identical(res_exact$assurance_method, res_exact$method))
})

test_that("print.dtasamplesize() runs without error on a bam_sample_size() result", {
  res <- bam_sample_size(prior_prev = c(4, 16), target_assurance = 0.80)
  expect_no_error(print(res))
  expect_no_error(capture.output(print(res)))
})

test_that("BAM method = \"exact\" (the default) never emits the spurious small-B Monte Carlo warning, regardless of B", {
  # method = "exact" ignores B entirely for the headline joint search (see
  # @details / @param B), so no warning about it should ever claim Monte
  # Carlo error in the reported (headline) assurance under this mode. This
  # must hold with the small-B advisory option turned back on (it is off for
  # the rest of the suite via setup.R).
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)

  msgs <- character(0)
  res <- withCallingHandlers(
    bam_sample_size(
      prior_prev = c(4, 16), target_assurance = 0.80, seed = 99, B = 100
    ),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # N_total / joint_assurance are exactly what the fully-converged default
  # call gives -- B = 100 must not have perturbed them under method = "exact".
  expect_identical(res$N_total, 678L)
  expect_equal(res$joint_assurance, 0.800349, tolerance = 1e-6)

  # None of the emitted warnings (if any) may claim Monte Carlo error in the
  # headline "assurance estimate" -- that phrasing is specific to the old,
  # spurious message and must not appear under method = "exact".
  expect_false(any(grepl(
    "Monte Carlo error in the assurance estimate", msgs, fixed = TRUE
  )))
})

test_that("BAM small-B warning under method = \"exact\" (when B is small enough to fire) names only the legacy per-arm fields it actually affects", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)

  expect_warning(
    bam_sample_size(
      prior_prev = c(4, 16), target_assurance = 0.80, seed = 99, B = 100,
      method = "exact"
    ),
    "legacy per-arm diagnostic fields"
  )
})

test_that("BAM small-B warning under method = \"monte_carlo\" keeps the package-wide generic message (B affects the headline result there)", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)

  # B = 100 makes the anti-noise Monte Carlo margin wide enough that this
  # small N_range may ALSO legitimately fail to converge (a second, unrelated
  # warning) -- that is expected behavior of the anti-noise rule at small B,
  # not something this test cares about. Capture every warning explicitly
  # (rather than relying on expect_warning()'s partial match) so an
  # incidental second warning cannot leak into the test reporter as noise;
  # this test only asserts that the small-B message itself is the generic,
  # package-wide one.
  msgs <- character(0)
  withCallingHandlers(
    bam_sample_size(
      prior_prev = c(4, 16), target_assurance = 0.80, seed = 99, B = 100,
      N_range = seq(670, 720, by = 2), method = "monte_carlo"
    ),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl(
    "Monte Carlo error in the assurance estimate", msgs, fixed = TRUE
  )))
})

test_that("documented default prevalence prior is Beta(6, 14) and gives N = 583 in exact mode", {
  # Guards against a repeat of the mismatch between the @param default doc
  # (c(6, 14)) and the worked examples in @details, which use prior_prev =
  # c(4, 16) (the accompanying article's worked example) and must not be
  # mistaken for the actual default. This checks the real formal default and
  # the N_total it produces under the fully-default call (method = "exact").
  # Runtime of the fully-default exact call was ~19s when this test was
  # written, well under 60s, so skip_on_cran() is not used here.
  expect_equal(eval(formals(bam_sample_size)$prior_prev), c(6, 14))

  result <- suppressWarnings(bam_sample_size(method = "exact"))
  expect_identical(result$N_total, 583L)
})

test_that("BAM issues no small-B warning at all when B >= 1000, in either mode", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)

  expect_no_warning(
    bam_sample_size(
      prior_prev = c(4, 16), target_assurance = 0.80, B = 5000,
      method = "exact"
    )
  )
  expect_no_warning(
    bam_sample_size(
      prior_prev = c(4, 16), target_assurance = 0.80, B = 5000, seed = 2026,
      N_range = seq(670, 720, by = 2), method = "monte_carlo"
    )
  )
})
