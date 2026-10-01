# M-01 (v0.6.6) contract tests for bam_sample_size()'s exact-mode joint
# search: the deterministic integer scan (N_range = NULL), the explicit-
# N_range contract, the shared "no crossing found" contract (both methods),
# RNG preservation, minimality of the returned N, degeneracy at N = 1 vs.
# N = 2, and print.dtasamplesize() on a non-crossing result. See
# R/bam_sample_size.R and NEWS.md ("dtasamplesize 0.6.6") for the defect
# this corrects (versions <= 0.6.5 built the automatic N_range from n_se,
# n_sp and N_total_P90 -- all Monte Carlo quantities -- so the supposedly
# B/seed-independent exact search could silently return different N_total
# depending on B and seed; and non-crossing searches returned max(N_range)
# as if it were a validated solution).

published_args <- list(
  prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
  delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
  method = "exact"
)

## ---------------------------------------------------------------------
## Invariance under B / seed (N_range = NULL, the auto integer scan)
## ---------------------------------------------------------------------

test_that("exact auto scan (N_range = NULL): N_total and joint_assurance are identical across B and seed, reduced battery", {
  results <- lapply(c(5L, 5000L), function(B) {
    lapply(c(1L, 2L), function(seed) {
      do.call(bam_sample_size, c(published_args, list(B = B, seed = seed)))
    })
  })
  flat <- unlist(results, recursive = FALSE)
  N_totals <- vapply(flat, function(r) r$N_total, integer(1))
  assurances <- vapply(flat, function(r) r$joint_assurance, double(1))
  expect_true(all(N_totals == 678L))
  expect_true(all(vapply(assurances, identical, logical(1), assurances[[1]])))
})

test_that("exact auto scan (N_range = NULL): N_total and joint_assurance are identical across B and seed, full battery", {
  skip_on_cran()
  B_grid <- c(5L, 10L, 50L, 5000L)
  seed_grid <- c(1L, 2L, 3L, 4L, 2026L)
  results <- list()
  for (B in B_grid) {
    for (seed in seed_grid) {
      results[[length(results) + 1L]] <- do.call(
        bam_sample_size, c(published_args, list(B = B, seed = seed))
      )
    }
  }
  N_totals <- vapply(results, function(r) r$N_total, integer(1))
  assurances <- vapply(results, function(r) r$joint_assurance, double(1))
  expect_true(all(N_totals == 678L))
  expect_true(all(assurances == assurances[1]))
  expect_equal(assurances[1], 0.800349, tolerance = 1e-6)
})

test_that("exact auto scan reproduces the harmonized-priors published case (N = 672), invariant to B/seed", {
  skip_on_cran()
  harmonized_args <- list(
    prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    method = "exact"
  )
  res_a <- do.call(bam_sample_size, c(harmonized_args, list(B = 5L, seed = 2L)))
  res_b <- do.call(bam_sample_size, c(harmonized_args, list(B = 5000L, seed = 2026L)))
  expect_identical(res_a$N_total, 672L)
  expect_identical(res_a$N_total, res_b$N_total)
  expect_identical(res_a$joint_assurance, res_b$joint_assurance)
  expect_equal(res_a$joint_assurance, 0.8002692084, tolerance = 1e-8)
})

## ---------------------------------------------------------------------
## Explicit N_range: B/seed irrelevant (already true pre-v0.6.6; kept as a
## regression guard), ascending/unique traversal
## ---------------------------------------------------------------------

test_that("explicit N_range (exact): B/seed do not change the result", {
  args <- c(published_args, list(N_range = 600:700))
  res_a <- do.call(bam_sample_size, c(args, list(B = 5L, seed = 1L)))
  res_b <- do.call(bam_sample_size, c(args, list(B = 5000L, seed = 999L)))
  expect_identical(res_a$N_total, res_b$N_total)
  expect_identical(res_a$joint_assurance, res_b$joint_assurance)
  expect_identical(res_a$N_total, 678L)
})

test_that("explicit N_range (exact) is searched as sort(unique(.)), ascending", {
  # A shuffled, duplicated N_range must give the exact same result as its
  # sorted-unique form: search_type = "user_N_range" always evaluates
  # sort(unique(as.integer(N_range))) regardless of the order/duplication
  # the caller supplied.
  res_sorted <- do.call(bam_sample_size, c(published_args, list(N_range = 600:700)))
  res_shuffled <- do.call(bam_sample_size, c(published_args, list(
    N_range = c(700:600, 650, 650, 678, 601)
  )))
  expect_identical(res_sorted$N_total, res_shuffled$N_total)
  expect_identical(res_sorted$joint_assurance, res_shuffled$joint_assurance)
  expect_identical(res_shuffled$N_range_used, 600:700)
  expect_identical(res_sorted$search_type, "user_N_range")
})

## ---------------------------------------------------------------------
## No-crossing contract: shared by both methods, never max(N_range)
## ---------------------------------------------------------------------

test_that("insufficient explicit N_range (exact): target_reached = FALSE, N_total/joint_assurance are NA, never max(N_range)", {
  expect_warning(
    res <- bam_sample_size(
      prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
      delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
      method = "exact", N_range = 100:200
    ),
    "No N in N_range achieved the target JOINT assurance"
  )
  expect_false(res$target_reached)
  expect_identical(res$N_total, NA_integer_)
  expect_identical(res$n_total, NA_integer_)
  expect_identical(res$joint_assurance, NA_real_)
  # The old defect: N_total == max(N_range) returned as if it were a
  # solution. Guard explicitly against a regression to that behavior.
  expect_false(isTRUE(res$N_total == 200L))
  expect_identical(res$N_at_max_assurance, 200L)
  expect_true(is.finite(res$max_assurance_evaluated))
})

test_that("insufficient N_max (exact, auto scan): target_reached = FALSE, N_total/joint_assurance are NA, never N_max", {
  expect_warning(
    res <- bam_sample_size(
      prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
      delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
      method = "exact", N_max = 200L
    ),
    "No integer N in 2:200"
  )
  expect_false(res$target_reached)
  expect_identical(res$N_total, NA_integer_)
  expect_identical(res$joint_assurance, NA_real_)
  expect_identical(res$search_type, "integer_scan_auto")
  expect_identical(res$N_max, 200L)
  expect_identical(res$N_range_used, 2L:200L)
  expect_false(isTRUE(res$N_total == 200L))
})

test_that("insufficient N_range (monte_carlo): target_reached = FALSE, N_total/joint_assurance are NA, never max(N_range)", {
  expect_warning(
    res <- bam_sample_size(
      prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
      delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
      method = "monte_carlo", N_range = 100:200, B = 2000, seed = 2026
    ),
    "No N in N_range achieved the target JOINT assurance"
  )
  expect_false(res$target_reached)
  expect_identical(res$N_total, NA_integer_)
  expect_identical(res$joint_assurance, NA_real_)
  expect_false(isTRUE(res$N_total == 200L))
})

## ---------------------------------------------------------------------
## RNG preservation under the exact auto scan (no RNG use at all)
## ---------------------------------------------------------------------

test_that("exact auto scan (N_range = NULL) never touches the caller's RNG stream", {
  set.seed(20260714)
  invisible(runif(3))
  expected <- runif(2)

  set.seed(20260714)
  invisible(runif(3))
  invisible(do.call(bam_sample_size, published_args))
  actual <- runif(2)

  expect_identical(actual, expected)
})

test_that(".Random.seed is restored byte for byte by the exact auto scan, including on non-crossing", {
  set.seed(4242)
  invisible(rnorm(1))
  before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  invisible(suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    method = "exact", N_max = 10L
  )))
  after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  expect_identical(after, before)
})

## ---------------------------------------------------------------------
## Crossing and minimality: A(N* - 1) < target <= A(N*), nothing smaller
## in 2:(N* - 1) also crosses
## ---------------------------------------------------------------------

test_that("the returned N* is minimal: A(N* - 1) < target <= A(N*), and no smaller integer crosses", {
  res <- do.call(bam_sample_size, published_args)
  N_star <- res$N_total
  expect_identical(N_star, 678L)

  a_se <- published_args$prior_se[1]; b_se <- published_args$prior_se[2]
  a_sp <- published_args$prior_sp[1]; b_sp <- published_args$prior_sp[2]
  ci_lower_q <- (1 - 0.95) / 2
  ci_upper_q <- 1 - ci_lower_q

  P_se <- dtasamplesize:::.bam_exact_width_prob(
    N_star, a_se, b_se, published_args$delta_se, ci_lower_q, ci_upper_q
  )
  P_sp <- dtasamplesize:::.bam_exact_width_prob(
    N_star, a_sp, b_sp, published_args$delta_sp, ci_lower_q, ci_upper_q
  )
  P_se[1] <- 0
  P_sp[1] <- 0

  assurance_at <- function(N) {
    dtasamplesize:::.bam_exact_joint_assurance(
      N, published_args$prior_prev[1], published_args$prior_prev[2], P_se, P_sp
    )
  }

  A_star <- assurance_at(N_star)
  A_below <- assurance_at(N_star - 1L)
  expect_gte(A_star, published_args$target_assurance)
  expect_lt(A_below, published_args$target_assurance)
  expect_identical(A_star, res$joint_assurance)

  all_assurances <- vapply(2:(N_star - 1L), assurance_at, double(1))
  expect_true(all(all_assurances < published_args$target_assurance))
})

## ---------------------------------------------------------------------
## Degeneracy: N = 1 always fails (A(1) = 0); n = 1 is NOT degenerate
## ---------------------------------------------------------------------

test_that("N = 1 is always a failure (A(1) = 0) under the exact calculation, regardless of priors", {
  res <- suppressWarnings(bam_sample_size(
    prior_se = c(1, 1), prior_sp = c(1, 1), prior_prev = c(1, 1),
    delta_se = 0.99, delta_sp = 0.99, target_assurance = 0.01,
    method = "exact", N_range = 1
  ))
  expect_identical(res$max_assurance_evaluated, 0)
})

test_that("the n = 1 cache entry is not forced to 0 and matches the direct calculation", {
  a <- 17; b <- 3; delta <- 0.14
  ci_lower_q <- 0.025; ci_upper_q <- 0.975
  P <- dtasamplesize:::.bam_exact_width_prob(5L, a, b, delta, ci_lower_q, ci_upper_q)
  # Direct calculation of P(width <= delta | n = 1): x in {0, 1}.
  post_a <- a + c(0, 1); post_b <- b + 1 - c(0, 1)
  width <- stats::qbeta(ci_upper_q, post_a, post_b) - stats::qbeta(ci_lower_q, post_a, post_b)
  log_pmf <- lchoose(1, c(0, 1)) + lbeta(post_a, post_b) - lbeta(a, b)
  direct_P1 <- sum(exp(log_pmf)[width <= delta])
  expect_equal(P[2], direct_P1) # P[1] is n = 0 (always 0 by convention, but
  # NOT forced inside .bam_exact_width_prob itself -- only by the caller);
  # P[2] is n = 1 and must come out of the raw enumeration unmodified, not
  # be hard-coded to any particular value (0 or otherwise) by
  # .bam_exact_width_prob() itself. With this prior/delta, a one-subject
  # posterior is in fact too wide to meet delta = 0.14 (P[2] happens to be
  # 0 here too) -- the point of this test is that P[2] is computed, not
  # that it is nonzero; see the N = 2 arm-of-size-1 test elsewhere (looser
  # priors/deltas) for a case where an n = 1 posterior DOES meet its
  # target width.
  expect_identical(P[1], 0) # n = 0 is NOT forced to 0 inside this raw
  # enumeration (only the caller zeroes it for the degenerate-arm
  # override) -- it is 0 here only because the PRIOR-only width also
  # exceeds delta, confirmed by computing it directly.
  post_a0 <- a; post_b0 <- b
  width0 <- stats::qbeta(ci_upper_q, post_a0, post_b0) - stats::qbeta(ci_lower_q, post_a0, post_b0)
  expect_identical(P[1], as.numeric(width0 <= delta))
})

## ---------------------------------------------------------------------
## .bam_exact_extend_cache(): block extension matches a from-scratch build
## ---------------------------------------------------------------------

test_that(".bam_exact_extend_cache() extending a cache matches building it from scratch", {
  a <- 17; b <- 3; delta <- 0.14
  ci_lower_q <- 0.025; ci_upper_q <- 0.975

  from_scratch <- dtasamplesize:::.bam_exact_width_prob(50L, a, b, delta, ci_lower_q, ci_upper_q)

  first_block <- dtasamplesize:::.bam_exact_width_prob(20L, a, b, delta, ci_lower_q, ci_upper_q)
  extended <- dtasamplesize:::.bam_exact_extend_cache(
    first_block, 20L, 50L, a, b, delta, ci_lower_q, ci_upper_q
  )
  expect_identical(extended, from_scratch)
})

## ---------------------------------------------------------------------
## print.dtasamplesize() on a non-crossing bam_sample_size() result
## ---------------------------------------------------------------------

test_that("print.dtasamplesize() reports a non-crossing BAM result as no solution, not as N = NA", {
  res <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    method = "exact", N_range = 100:200
  ))
  out <- capture.output(print(res))
  expect_true(any(grepl("Target assurance NOT reached", out)))
  expect_true(any(grepl("max_assurance_evaluated", out)))
  # Must not print a bare "N_total: NA" line as if NA were a candidate value
  # under the normal "N_total:" label used by the converged branch.
  expect_false(any(grepl("^\\s*N_total: NA", out)))
})
