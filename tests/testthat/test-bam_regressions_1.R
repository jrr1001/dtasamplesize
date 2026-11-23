# Regression tests for the bam_* functions of version 0.7.0 (part 1).

collect_w <- function(expr) {
  w <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = w)
}

## ---------------------------------------------------------------------
## Monte Carlo with a marginal objective
## ---------------------------------------------------------------------

headline <- function(r) {
  r[c("N_total", "joint_assurance", "assurance_mcse", "N_range_used",
      "search_type", "max_assurance_evaluated", "N_at_max_assurance",
      "target_reached", "n_diseased", "assurance_se")]
}

test_that("Monte Carlo, assurance_for = 'se', N_range = NULL: the unused Sp arm does not matter", {
  common <- list(prior_se = c(25.9, 2.1), prior_prev = c(29, 98),
                 delta_se = 0.16, assurance_for = "se",
                 method = "monte_carlo", B = 2000, n_range = 20:200,
                 seed = 5)
  a <- suppressWarnings(do.call(bam_sample_size, c(common, list(
    prior_sp = c(2, 2), delta_sp = 0.10))))
  b <- suppressWarnings(do.call(bam_sample_size, c(common, list(
    prior_sp = c(40, 3), delta_sp = 0.30))))
  expect_identical(headline(a), headline(b))
  # the grid does not start at the unused arm's requirement
  expect_identical(a$search_type, "monte_carlo_auto")
  expect_lt(min(a$N_range_used), 120L)
})

test_that("Monte Carlo, assurance_for = 'sp', N_range = NULL: the unused Se arm does not matter", {
  # Mirror image of the VAP example: the Sp arm is the minority one
  common <- list(prior_sp = c(25.9, 2.1), prior_prev = c(98, 29),
                 delta_sp = 0.16, assurance_for = "sp",
                 method = "monte_carlo", B = 2000, n_range = 20:200,
                 seed = 5)
  a <- suppressWarnings(do.call(bam_sample_size, c(common, list(
    prior_se = c(2, 2), delta_se = 0.04))))
  b <- suppressWarnings(do.call(bam_sample_size, c(common, list(
    prior_se = c(40, 3), delta_se = 0.30))))
  keep <- setdiff(names(headline(a)), c("n_diseased", "assurance_se"))
  expect_identical(headline(a)[keep], headline(b)[keep])
  expect_true(a$target_reached)
  expect_lt(min(a$N_range_used), 120L)
})

test_that("VAP example, Monte Carlo marginal N is close to the exact N (104)", {
  skip_on_cran()
  args <- list(prior_se = c(25.9, 2.1), prior_prev = c(29, 98),
               delta_se = 0.16, assurance_for = "se", n_range = 20:200,
               B = 5000, seed = 2026)
  mc <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    method = "monte_carlo"))))
  expect_true(mc$target_reached)
  expect_lte(mc$N_total, 140L)
  ex <- suppressWarnings(do.call(bam_sample_size, c(
    args, list(method = "exact", N_range = mc$N_total))))
  expect_gte(ex$max_assurance_evaluated, 0.79)
  # the exact answer itself
  exa <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    method = "exact"))))
  expect_identical(exa$N_total, 104L)
})

test_that("analogous Sp-only example (mirror of VAP), Monte Carlo N close to exact", {
  skip_on_cran()
  args <- list(prior_sp = c(25.9, 2.1), prior_prev = c(98, 29), delta_sp = 0.16,
               prior_se = c(2, 2), delta_se = 0.04, assurance_for = "sp",
               n_range = 20:200, B = 5000, seed = 2026)
  exa <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    method = "exact"))))
  expect_identical(exa$N_total, 104L)   # same distribution of n_nd as n_d in VAP
  mc <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    method = "monte_carlo"))))
  expect_true(mc$target_reached)
  expect_lte(mc$N_total, 140L)
  chk <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    method = "exact", N_range = mc$N_total))))
  expect_gte(chk$max_assurance_evaluated, 0.79)
})

test_that("marginal objective does not warn about the unused arm's legacy search", {
  args <- list(prior_se = c(25.9, 2.1), prior_prev = c(29, 98),
               delta_se = 0.16, prior_sp = c(2, 2), delta_sp = 0.001,
               n_range = 20:200, B = 1000, seed = 1)
  for (m in c("exact", "monte_carlo")) {
    extra <- if (m == "exact") list(N_range = 90:110) else list()
    w <- collect_w(do.call(bam_sample_size, c(args, extra, list(
      assurance_for = "se", method = m))))$warnings
    expect_false(any(grepl("for Sp", w)), info = m)
  }
  # the joint objective still warns (the Sp arm is used)
  wj <- collect_w(do.call(bam_sample_size, c(args, list(
    assurance_for = "joint", method = "exact", N_range = 90:110))))$warnings
  expect_true(any(grepl("for Sp", wj)))
})

test_that("joint Monte Carlo keeps the legacy grid from max(20, n_se, n_sp)", {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.30, delta_sp = 0.30, target_assurance = 0.85,
    B = 300, n_range = 20:60, N_range = NULL, seed = 7,
    method = "monte_carlo"))
  expect_identical(r$search_type, "monte_carlo_auto")
  expect_identical(r$N_range_used[1], as.integer(max(20L, r$n_diseased,
                                                    r$n_non_diseased)))
})

## ---------------------------------------------------------------------
## validation of alpha_ci
## ---------------------------------------------------------------------

test_that("alpha_ci must be a single number in (0, 1)", {
  for (bad in list(0, 1, -0.1, 1.5, NA_real_, "0.95", c(0.9, 0.95), NULL)) {
    expect_error(
      bam_sample_size(alpha_ci = bad, N_range = 100L, n_range = 20:30,
                      B = 1000),
      "alpha_ci", info = paste(format(bad), collapse = ","))
  }
  expect_error(bam_assurance_curve(N_range = 10L, alpha_ci = 0), "alpha_ci")
  expect_error(bam_assurance_curve(N_range = 10L, alpha_ci = c(0.9, 0.95)),
               "alpha_ci")
})

test_that("lower_half_width requires alpha_ci > 0.5 (both functions)", {
  for (a in c(0.5, 0.3, 0.01)) {
    expect_error(
      bam_sample_size(alpha_ci = a, width_type = "lower_half_width",
                      N_range = 100L, n_range = 20:30, B = 1000),
      "median|0\\.5")
    expect_error(
      bam_assurance_curve(N_range = 100L, alpha_ci = a,
                          width_type = "lower_half_width"),
      "median|0\\.5")
  }
  # a two-sided criterion at level <= 0.5 is meaningful and still allowed
  expect_error(
    suppressWarnings(bam_sample_size(alpha_ci = 0.4, N_range = 100L,
                                     n_range = 20:30, B = 1000)),
    NA)
  # and 0.6 is fine for the one-sided criterion
  expect_error(
    suppressWarnings(bam_sample_size(alpha_ci = 0.6, N_range = 100L,
                                     width_type = "lower_half_width",
                                     n_range = 20:30, B = 1000)),
    NA)
})

## ---------------------------------------------------------------------
## warning about the EFFECTIVE analysis prior
## ---------------------------------------------------------------------

test_that("a design prior that already meets the target warns even when passed by default/design", {
  # Beta(170, 30): prior-only 95% full width ~ 0.099 <= 0.14
  r <- collect_w(bam_sample_size(
    prior_se = c(170, 30), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, method = "exact",
    n_range = 20:30, B = 1000, N_range = 300:400))
  expect_true(any(grepl("analysis prior", r$warnings) &
                    grepl("equal to the design prior", r$warnings) &
                    grepl("Se", r$warnings)))
  # same with assurance_for = 'se'; an unused arm does not warn
  r2 <- collect_w(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(180, 20), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.14, method = "exact", assurance_for = "se",
    n_range = 20:30, B = 1000, N_range = 300:400))
  expect_false(any(grepl("analysis prior for Sp", r2$warnings)))
  # explicit analysis prior: the old message (no "equal to the design prior")
  r3 <- collect_w(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, method = "exact",
    analysis_prior_se = c(500, 100),
    n_range = 20:30, B = 1000, N_range = 300:400))
  w3 <- r3$warnings[grepl("analysis prior for Se", r3$warnings)]
  expect_length(w3, 1L)
  expect_false(grepl("equal to the design prior", w3))
})

test_that("the analysis prior taken from a priors object also warns", {
  po <- bam_prior_from_pilot(tp = 90, fn = 10, tn = 80, fp = 20,
                             base_prior_se = c(1, 1), base_prior_sp = c(1, 1),
                             base_prior_prev = c(1, 1), a0 = 1)
  # With a very informative pilot the updated (analysis) prior of Se, Beta(91, 11),
  # has a prior-only width below a lenient delta.
  r <- collect_w(bam_sample_size(
    priors = po, delta_se = 0.30, delta_sp = 0.30, method = "exact",
    n_range = 20:30, B = 1000, N_range = 100:200))
  expect_true(any(grepl("analysis prior for Se", r$warnings)))
})

test_that("bam_assurance_curve() warns for the effective analysis prior too", {
  r <- collect_w(bam_assurance_curve(
    prior_se = c(170, 30), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, N_range = 300:310))
  expect_true(any(grepl("analysis prior for Se", r$warnings) &
                    grepl("equal to the design prior", r$warnings)))
  # no warning with the defaults
  r0 <- collect_w(bam_assurance_curve(N_range = 300:310))
  expect_false(any(grepl("analysis prior", r0$warnings)))
})

## ---------------------------------------------------------------------
## N_range is reserved in bam_sensitivity()
## ---------------------------------------------------------------------

test_that("bam_sensitivity() rejects N_range with a clear error", {
  scen <- list(a = list(prior_se = c(17, 3)))
  expect_error(
    bam_sensitivity(scen, prior_prev = c(4, 16), delta_se = 0.30,
                    delta_sp = 0.30, N_range = 100:200),
    "N_range")
  expect_error(
    bam_sensitivity_local(
      base_prior_se = c(9.9, 1.1), base_prior_sp = c(2, 2),
      base_prior_prev = c(12, 43), epsilon = 0.00354, n_phi = 3,
      N_range = 100:200),
    "N_range")
})

## ---------------------------------------------------------------------
## precision of the Hellinger contour
## ---------------------------------------------------------------------

test_that("the contour warns when the numerical precision is not enough", {
  expect_warning(
    ct <- bam_prior_contour(c(1e5, 1e4), epsilon = 1e-6, distance = "H",
                            n_phi = 12L),
    "precision")
  expect_true(any(attr(ct, "imprecise")))
  expect_length(attr(ct, "imprecise"), 12L)
})

test_that("the VAP contour does not warn and is precise", {
  expect_warning(
    ct <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354, n_phi = 36L),
    NA)
  expect_false(any(attr(ct, "imprecise")))
  expect_lt(max(abs(ct$distance - 0.00354)) / 0.00354, 1e-6)
})

test_that("bam_sensitivity_local() reports the discarded directions", {
  skip_on_cran()
  # Beta(1000, 100), H distance, epsilon 1e-4: the estimated rounding error of
  # the distance exceeds the tolerance in every direction (deterministic, see
  # test-bam_regressions_2.R), so no direction is left and the error says so and
  # counts them. The partial case (some directions discarded, some kept) is
  # tested in test-bam_regressions_2.R with a mocked contour.
  expect_error(suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(1000, 100), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 1e-4,
    distance = "H", n_phi = 24L, n_ref = 50L, delta_se = 0.30,
    assurance_for = "se")),
    "24 with insufficient numerical precision")
  # all directions imprecise: an error that says why
  expect_error(suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(1e5, 1e4), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 1e-6,
    distance = "H", n_phi = 6L, n_ref = 50L, delta_se = 0.30,
    assurance_for = "se")), "numerical precision")
})

test_that("a regular local sensitivity discards no direction", {
  skip_on_cran()
  res <- suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(9.9, 1.1), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
    n_phi = 4L, n_ref = 106L, delta_se = 0.16, assurance_for = "se"))
  expect_identical(res$summary$n_dir_no_root, 0L)
  expect_identical(res$summary$n_dir_imprecise, 0L)
  expect_identical(nrow(res$table), 4L)
})

## ---------------------------------------------------------------------
## Buderer comparator under the one-sided criterion
## ---------------------------------------------------------------------

test_that("lower_half_width uses z = qnorm(alpha_ci) in the Buderer comparator", {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(25, 7), delta_se = 0.10, width_type = "lower_half_width",
    N_range = 100:150, n_range = 20:60, B = 1000))
  expect_identical(r$buderer_n_se,
                   buderer_n(25 / 32, 0.10, alpha = 2 * (1 - 0.95)))
  expect_equal(r$buderer_n_se, buderer_n(25 / 32, 0.10, alpha = 0.10))
  # alpha_ci = 0.90: z = qnorm(0.90)
  r9 <- suppressWarnings(bam_sample_size(
    prior_se = c(25, 7), delta_se = 0.10, width_type = "lower_half_width",
    alpha_ci = 0.90, N_range = 100:150, n_range = 20:60, B = 1000))
  expect_identical(r9$buderer_n_se,
                   ceiling(stats::qnorm(0.90)^2 * (25 / 32) * (7 / 32) / 0.10^2))
})

test_that("full_width keeps the 95% z whatever alpha_ci is (inherited behaviour)", {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(25, 7), delta_se = 0.10, alpha_ci = 0.90,
    N_range = 100:150, n_range = 20:60, B = 1000))
  expect_identical(r$buderer_n_se, buderer_n(25 / 32, 0.10 / 2))
})

## ---------------------------------------------------------------------
## print method
## ---------------------------------------------------------------------

test_that("the print labels the assurance as marginal when assurance_for is se / sp", {
  se <- suppressWarnings(bam_sample_size(
    prior_se = c(25.9, 2.1), prior_prev = c(29, 98), delta_se = 0.16,
    assurance_for = "se", n_range = 20:60, B = 1000, N_range = 90:110))
  out <- paste(capture.output(print(se)), collapse = "\n")
  expect_match(out, "marginal Se assurance:")
  expect_false(grepl("joint_assurance:", out, fixed = TRUE))
  sp <- suppressWarnings(bam_sample_size(
    prior_sp = c(18, 2), prior_prev = c(60, 40), delta_sp = 0.20,
    assurance_for = "sp", n_range = 20:60, B = 1000, N_range = 40:80))
  out2 <- paste(capture.output(print(sp)), collapse = "\n")
  expect_match(out2, "marginal Sp assurance:")
  # the joint mode prints exactly what it printed before
  jo <- suppressWarnings(bam_sample_size(
    prior_prev = c(4, 16), n_range = 20:60, B = 1000, N_range = 600:700))
  out3 <- paste(capture.output(print(jo)), collapse = "\n")
  expect_match(out3, "joint_assurance:", fixed = TRUE)
})

test_that("print says 'same as design' when the given analysis prior equals the design prior", {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), analysis_prior_se = c(17, 3), prior_prev = c(4, 16),
    n_range = 20:60, B = 1000, N_range = 600:700))
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_match(out, "same as design", ignore.case = TRUE)
  expect_false(grepl("differ", out, ignore.case = TRUE))
  # a priors object with use_in = "both" has analysis = design
  po <- bam_prior_from_pilot(tp = 8, fn = 2, tn = 9, fp = 1,
                             base_prior_se = c(1, 1), base_prior_sp = c(1, 1),
                             base_prior_prev = c(1, 1))
  r2 <- suppressWarnings(bam_sample_size(
    priors = po, delta_se = 0.30, delta_sp = 0.30,
    n_range = 20:60, B = 1000, N_range = 100:200))
  out2 <- paste(capture.output(print(r2)), collapse = "\n")
  expect_match(out2, "same as design", ignore.case = TRUE)
  expect_false(grepl("differ", out2, ignore.case = TRUE))
  # a different analysis prior is still reported as different
  r3 <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), analysis_prior_se = c(1, 1), prior_prev = c(4, 16),
    n_range = 20:60, B = 1000, N_range = 600:700))
  out3 <- paste(capture.output(print(r3)), collapse = "\n")
  expect_match(out3, "differ", ignore.case = TRUE)
})
