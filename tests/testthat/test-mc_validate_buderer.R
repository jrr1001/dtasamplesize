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

# --- exact (binomial enumeration) mode ----------------------------------------
# Reference values: validation/width_attainment_exact.R (a separate script
# written by the same team); full enumeration of x = 0..n, no Monte Carlo.

exact_row <- function(res, m) res$results[res$results$ci_method == m, ]

test_that("method defaults to monte_carlo and is validated", {
  res <- mc_validate_buderer(B = 200, seed = 1)
  expect_identical(res$method_type, "monte_carlo")
  expect_identical(res$method, "Monte Carlo Validation of Buderer Formula")
  expect_error(mc_validate_buderer(method = "bogus"), "should be one of")
})

test_that("exact mode reproduces the reference of width_attainment_exact.R at Se = 0.85, d = 0.07, n = 100", {
  res <- mc_validate_buderer(Se = 0.85, d = 0.07, method = "exact")
  expect_s3_class(res, "dtasamplesize")
  expect_identical(res$method_type, "exact")
  expect_equal(res$n_diseased, 100)
  wald <- exact_row(res, "wald")
  wilson <- exact_row(res, "wilson")
  expect_equal(round(wald$P_width_target, 4), 0.5683)
  expect_equal(round(wilson$P_width_target, 4), 0.5683)
  expect_equal(round(wald$coverage, 7), 0.9330186)
  expect_equal(round(wilson$coverage, 7), 0.9332462)
  # full precision against the CSV values of width_attainment_exact.R
  expect_equal(wald$P_width_target, 0.568315165125715, tolerance = 1e-12)
  expect_equal(wald$coverage, 0.933018636072062, tolerance = 1e-12)
  expect_equal(wilson$coverage, 0.933246217165038, tolerance = 1e-12)
  # P(width <= target) is not the coverage
  expect_false(isTRUE(all.equal(wald$P_width_target, wald$coverage)))
})

test_that("exact mode reproduces selected points of the reference curve (n = 40..200)", {
  ref <- data.frame(
    n = c(40, 60, 80, 100, 113, 114, 115, 116, 117, 150, 200),
    wald = c(0.0485986656928547, 0.0967985087269186, 0.221119227532048,
             0.568315165125715, 0.754485617571245, 0.815631710633949,
             0.804508325619351, 0.793062947880646, 0.846801721453803,
             0.999672728266301, 1),
    wilson = c(0.0121067806608898, 0.0423720469869597, 0.221119227532048,
               0.568315165125715, 0.754485617571245, 0.815631710633949,
               0.804508325619351, 0.856285034437302, 0.846801721453803,
               0.999844174349534, 1)
  )
  for (i in seq_len(nrow(ref))) {
    r <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = ref$n[i],
                             method = "exact")
    expect_equal(exact_row(r, "wald")$P_width_target, ref$wald[i],
                 tolerance = 1e-12, info = paste("wald n =", ref$n[i]))
    expect_equal(exact_row(r, "wilson")$P_width_target, ref$wilson[i],
                 tolerance = 1e-12, info = paste("wilson n =", ref$n[i]))
  }
})

test_that("first n with P >= 0.80 and smallest n stable to 200", {
  ns <- 40:200
  P <- vapply(ns, function(n) {
    r <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = n,
                             method = "exact")
    c(exact_row(r, "wald")$P_width_target, exact_row(r, "wilson")$P_width_target)
  }, numeric(2))
  first_ge <- function(p) ns[which(p >= 0.80 - 1e-12)[1L]]
  stable_ge <- function(p) {
    bad <- which(p < 0.80 - 1e-12)
    if (length(bad) == 0L) ns[1L] else ns[max(bad) + 1L]
  }
  expect_equal(c(first_ge(P[1, ]), first_ge(P[2, ])), c(114, 114))
  expect_equal(c(stable_ge(P[1, ]), stable_ge(P[2, ])), c(117, 114))
})

test_that("exact mode reproduces Se = 0.98, d = 0.05, n = 31 (Wald truncated)", {
  res <- mc_validate_buderer(Se = 0.98, d = 0.05, method = "exact")
  expect_equal(res$n_diseased, 31)
  expect_equal(round(exact_row(res, "wald")$P_width_target, 4), 0.8728)
  expect_equal(exact_row(res, "wald")$P_width_target, 0.872774911011899,
               tolerance = 1e-12)
  # Wilson never reaches the target here (cell table: P_width = 0)
  expect_equal(exact_row(res, "wilson")$P_width_target, 0)
})

test_that("no enumerated width lies within 1e-12 of the target in the test cases", {
  min_gap <- function(Se, d, n, m) {
    z <- stats::qnorm(0.975); x <- 0:n; ph <- x / n
    if (m == "wald") {
      h <- z * sqrt(ph * (1 - ph) / n)
      w <- pmin(1, ph + h) - pmax(0, ph - h)
    } else {
      den <- 1 + z^2 / n; ctr <- (ph + z^2 / (2 * n)) / den
      mar <- z * sqrt((ph * (1 - ph) + z^2 / (4 * n)) / n) / den
      w <- pmin(1, ctr + mar) - pmax(0, ctr - mar)
    }
    pr <- stats::dbinom(x, n, Se)
    min(abs(w[pr > 0] - 2 * d))
  }
  for (m in c("wald", "wilson")) {
    expect_gt(min_gap(0.85, 0.07, 100, m), 1e-12)
    expect_gt(min_gap(0.98, 0.05, 31, m), 1e-12)
    for (n in 40:200) expect_gt(min_gap(0.85, 0.07, n, m), 1e-12)
  }
})

test_that("exact mode does not depend on B or seed and leaves the RNG alone", {
  r1 <- mc_validate_buderer(method = "exact", B = 10, seed = 1)
  r2 <- mc_validate_buderer(method = "exact", B = 99999, seed = 7)
  r3 <- mc_validate_buderer(method = "exact")
  expect_identical(r1$results, r2$results)
  expect_identical(r1$results, r3$results)
  set.seed(123)
  before <- .Random.seed
  invisible(mc_validate_buderer(method = "exact"))
  expect_identical(.Random.seed, before)
})

test_that("exact mode does not warn about small B", {
  op <- options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(op), add = TRUE)
  expect_no_warning(mc_validate_buderer(method = "exact", B = 5))
  expect_warning(mc_validate_buderer(method = "monte_carlo", B = 5, seed = 1),
                 "B = 5 is small")
})

test_that("exact quantiles are discrete-distribution quantiles", {
  res <- mc_validate_buderer(Se = 0.85, d = 0.07, method = "exact")
  n <- 100; z <- stats::qnorm(0.975); x <- 0:n; ph <- x / n
  h <- z * sqrt(ph * (1 - ph) / n)
  w <- pmin(1, ph + h) - pmax(0, ph - h)
  pr <- stats::dbinom(x, n, 0.85)
  ord <- order(w)
  cw <- cumsum(pr[ord])
  q <- function(p) w[ord][which(cw >= p - 1e-12)[1L]]
  wald <- exact_row(res, "wald")
  expect_equal(wald$median_width, q(0.5), tolerance = 1e-12)
  expect_equal(wald$P75_width, q(0.75), tolerance = 1e-12)
  expect_equal(wald$P90_width, q(0.90), tolerance = 1e-12)
  expect_true(wald$median_width <= wald$P75_width &&
                wald$P75_width <= wald$P90_width)
  # the quantile is a realisable width
  expect_true(any(abs(w - wald$P90_width) < 1e-14))
})

test_that("exact mode print states that figures are exact", {
  res <- mc_validate_buderer(method = "exact", ci_method = "wald")
  out <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(out, "Exact")
  expect_false(grepl("Monte Carlo", out))
  expect_false(grepl("B:", out, fixed = TRUE))
  res_mc <- mc_validate_buderer(B = 200, seed = 1, ci_method = "wald")
  out_mc <- paste(utils::capture.output(print(res_mc)), collapse = "\n")
  expect_match(out_mc, "Monte Carlo")
})

test_that("exact mode honours alpha, ci_method and explicit n_diseased", {
  r <- mc_validate_buderer(Se = 0.9, d = 0.05, n_diseased = 50, alpha = 0.10,
                           ci_method = "wilson", method = "exact")
  expect_equal(nrow(r$results), 1)
  expect_equal(r$results$ci_method, "wilson")
  expect_equal(r$n_diseased, 50)
  expect_equal(r$alpha, 0.10)
  expect_true(all(c("ci_method", "P_width_target", "coverage", "median_width",
                    "P75_width", "P90_width") %in% names(r$results)))
})

test_that("Monte Carlo path keeps the contract of the earlier tests (v0.6.9)", {
  # Same calls and bounds as the pre-existing tests; the bit-for-bit
  # comparison against v0.6.9 is done in a separate process with the v0.6.9
  # library (see validation/compat_v069_check.R).
  res <- mc_validate_buderer(Se = 0.85, d = 0.07, B = 5000, seed = 2026)
  wald <- exact_row(res, "wald")
  expect_true(wald$coverage < 0.95 && wald$coverage > 0.90)
  expect_true(wald$P_width_target > 0.50 && wald$P_width_target < 0.65)
  res2 <- mc_validate_buderer(Se = 0.85, d = 0.07, B = 5000, seed = 2026)
  expect_identical(res$results, res2$results)
  expect_identical(res$B, 5000)
})

test_that("input errors in exact mode", {
  expect_error(mc_validate_buderer(Se = 0, method = "exact"))
  expect_error(mc_validate_buderer(Se = 1, method = "exact"))
  expect_error(mc_validate_buderer(d = 0.5, method = "exact"))
  expect_error(mc_validate_buderer(d = 0, method = "exact"))
  expect_error(mc_validate_buderer(alpha = 1, method = "exact"))
  expect_error(mc_validate_buderer(ci_method = "x", method = "exact"))
  expect_error(mc_validate_buderer(method = c("exact", "monte_carlo", "x")))
  expect_error(mc_validate_buderer(n_diseased = 0, method = "exact"))
  expect_error(mc_validate_buderer(n_diseased = c(50, 60), method = "exact"))
  expect_error(mc_validate_buderer(n_diseased = 10.5, method = "exact"))
  expect_error(mc_validate_buderer(n_diseased = NA_real_, method = "exact"))
})
