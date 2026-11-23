# Tests for bam_prior_data_conflict()

# Study of Wilson et al. 2022 (VAP): 53 diseased of 150; TP 51 of 53;
# TN 42, FP 55 (Table 1C). The Sp prior below is arbitrary: nothing is
# asserted about the specificity row.
vap_conf <- function(...) {
  args <- list(tp = 51, fn = 2, tn = 42, fp = 55,
               prior_se = c(25.9, 2.1), prior_prev = c(29, 98),
               prior_sp = c(21, 36))
  suppressWarnings(do.call(bam_prior_data_conflict,
                           utils::modifyList(args, list(...))))
}

get_row <- function(res, quantity_pat, ref = NULL) {
  r <- res[grepl(quantity_pat, res$quantity), , drop = FALSE]
  if (!is.null(ref)) r <- r[r$reference == ref, , drop = FALSE]
  r
}

test_that("VAP example: tail probabilities match the reference values (4 decimals)", {
  res <- vap_conf()
  expect_s3_class(res, "dta_bam_conflict")
  expect_s3_class(res, "data.frame")
  prev <- get_row(res, "diseased", "prevalence")
  expect_equal(nrow(prev), 1)
  expect_equal(prev$x_obs, 53); expect_equal(prev$n, 150)
  expect_equal(prev$a, 29); expect_equal(prev$b, 98)
  expect_equal(round(prev$p_lower, 4), 0.9913)
  expect_equal(round(prev$p_upper, 4), 0.0119)
  se <- get_row(res, "true positives")
  expect_equal(se$x_obs, 51); expect_equal(se$n, 53)
  expect_equal(round(se$p_lower, 4), 0.7589)
  expect_equal(round(se$p_upper, 4), 0.3897)
  expect_true(prev$flag_upper)
  expect_false(prev$flag_lower)
  expect_equal(prev$p_two_sided, min(1, 2 * min(prev$p_lower, prev$p_upper)))
  expect_true(prev$flag_two_sided)
  expect_false(se$flag_two_sided)
})

test_that("the two tails sum to 1 + P(X = x_obs)", {
  res <- vap_conf()
  for (i in seq_len(nrow(res))) {
    pm <- dtasamplesize:::.bb_pmf(res$x_obs[i], res$n[i], res$a[i], res$b[i])
    expect_equal(res$p_lower[i] + res$p_upper[i], 1 + pm, tolerance = 1e-12)
  }
})

test_that("pmf sums to one in several cases", {
  cases <- list(c(10, 1, 1), c(53, 25.9, 2.1), c(150, 29, 98),
                c(200, 0.5, 0.5), c(1, 3, 7), c(97, 21, 36))
  for (cs in cases) {
    pm <- dtasamplesize:::.bb_pmf(0:cs[1], cs[1], cs[2], cs[3])
    expect_equal(sum(pm), 1, tolerance = 1e-12)
  }
})

test_that("distinct analysis priors give two rows per arm and differ", {
  res <- vap_conf(analysis_prior_se = c(9.9, 1.1),
                  analysis_prior_sp = c(1, 1))
  se <- get_row(res, "true positives")
  sp <- get_row(res, "true negatives")
  expect_equal(sort(se$reference), c("analysis", "design"))
  expect_equal(sort(sp$reference), c("analysis", "design"))
  expect_false(isTRUE(all.equal(se$p_lower[se$reference == "analysis"],
                                se$p_lower[se$reference == "design"])))
  expect_equal(se$a[se$reference == "analysis"], 9.9)
  expect_equal(se$a[se$reference == "design"], 25.9)
})

test_that("equal priors give a single row per arm", {
  res <- vap_conf()
  expect_equal(nrow(get_row(res, "true positives")), 1)
  expect_equal(nrow(get_row(res, "true negatives")), 1)
  expect_true(all(get_row(res, "true positives")$same_as_design))
})

test_that("flat prior with n = 10: warning and no flag", {
  expect_warning(
    res <- bam_prior_data_conflict(tp = 5, fn = 0, tn = 5, fp = 0,
                                   prior_se = c(1, 1), prior_sp = c(1, 1)),
    "cannot flag")
  expect_false(any(res$flag_lower | res$flag_upper | res$flag_two_sided))
  expect_equal(res$min_attainable_p, rep(2 / 6, nrow(res)),
               tolerance = 1e-12)
  res2 <- suppressWarnings(
    bam_prior_data_conflict(tp = 10, fn = 0, tn = 10, fp = 0,
                            prior_se = c(1, 1), prior_sp = c(1, 1)))
  expect_equal(res2$p_upper[1], 1 / 11, tolerance = 1e-12)
  expect_false(any(res2$flag_upper))
})

test_that("flat prior, n = 30: only extremes flagged (30/30 -> 1/31)", {
  res <- suppressWarnings(
    bam_prior_data_conflict(tp = 30, fn = 0, tn = 15, fp = 15,
                            prior_se = c(1, 1), prior_sp = c(1, 1)))
  se <- get_row(res, "true positives")
  expect_equal(se$p_upper, 1 / 31, tolerance = 1e-12)
  expect_true(se$flag_upper)
  expect_equal(se$min_attainable_p, 2 / 31, tolerance = 1e-12)
})

test_that("flat prior, n = 20: a single tail can flag but the two-sided rule cannot", {
  expect_warning(
    res <- bam_prior_data_conflict(tp = 20, fn = 0, tn = 10, fp = 10,
                                   prior_se = c(1, 1), prior_sp = c(1, 1)),
    "two-sided check cannot flag")
  se <- get_row(res, "true positives")
  expect_equal(se$min_attainable_tail, 1 / 21, tolerance = 1e-12)
  expect_true(se$flag_upper)
  expect_false(se$flag_two_sided)
})

test_that("a flat prior triggers a warning even when flags are possible", {
  expect_warning(
    bam_prior_data_conflict(tp = 30, fn = 0, tn = 15, fp = 15,
                            prior_se = c(1, 1), prior_sp = c(1, 1)),
    "flat")
})

test_that("case_control has no prevalence row and says why", {
  res <- vap_conf(design = "case_control")
  expect_false(any(res$reference == "prevalence"))
  expect_true(any(grepl("prevalence", attr(res, "notes"))))
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "case_control")
})

test_that("priors object is accepted and cannot be mixed", {
  pr <- bam_prior_from_pilot(16, 1, 20, 35, c(9.9, 1.1), c(1, 1), c(12, 43),
                             use_in = "design")
  res <- suppressWarnings(bam_prior_data_conflict(51, 2, 42, 55, priors = pr))
  se <- get_row(res, "true positives")
  expect_equal(sort(se$reference), c("analysis", "design"))
  expect_error(bam_prior_data_conflict(51, 2, 42, 55, priors = pr,
                                       prior_se = c(1, 1)), "priors")
  expect_error(bam_prior_data_conflict(51, 2, 42, 55,
                                       priors = list(prior_se = c(1, 1))),
               "dta_bam_priors")
})

test_that("input errors", {
  ok <- list(prior_se = c(2, 2), prior_sp = c(2, 2), prior_prev = c(2, 2))
  f <- function(...) do.call(bam_prior_data_conflict,
                             utils::modifyList(c(list(tp = 5, fn = 1, tn = 5,
                                                      fp = 1), ok), list(...)))
  expect_error(f(tp = -1), "tp")
  expect_error(f(fn = 0.5), "fn")
  expect_error(f(alpha = 0), "alpha")
  expect_error(f(alpha = 1), "alpha")
  expect_error(f(prior_se = c(0, 1)), "prior_se")
  expect_error(f(prior_sp = c(1, NA)), "prior_sp")
  expect_error(f(design = "zz"))
  expect_error(bam_prior_data_conflict(5, 1, 5, 1), "prior")
})

test_that("print method runs and states the flagging rule", {
  res <- vap_conf()
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "alpha")
  expect_match(out, "P\\(X <= x_obs\\)")
})
