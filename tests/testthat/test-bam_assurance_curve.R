# bam_assurance_curve() (v0.7.0): exact curve A(N).
#
# The curve must match bam_sample_size(method = "exact",
# N_range = N)$max_assurance_evaluated BIT FOR BIT at every N (with a single
# N evaluated, max_assurance_evaluated is A(N)).

p1 <- list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
           delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80)

A_ref <- function(N, args) {
  vapply(N, function(n) {
    suppressWarnings(do.call(bam_sample_size, c(args, list(
      method = "exact", N_range = n, n_range = 20:30, B = 1000
    ))))$max_assurance_evaluated
  }, numeric(1))
}

test_that("A(N) identical to bam_sample_size() for P1 around the crossing", {
  skip_on_cran()
  Ns <- c(600L, 625L, 650L, 677L, 678L, 679L, 700L)
  cv <- do.call(bam_assurance_curve, c(p1, list(N_range = 600:700)))
  expect_s3_class(cv, "dta_bam_curve")
  got <- cv$curve$assurance[match(Ns, cv$curve$N)]
  expect_identical(got, A_ref(Ns, p1))
})

test_that("identical for width_type = 'lower_half_width'", {
  skip_on_cran()
  a <- list(prior_se = c(25, 7), prior_sp = c(30, 2),
            prior_prev = c(13.56, 122.06), delta_se = 0.10, delta_sp = 0.05,
            target_assurance = 0.80, width_type = "lower_half_width")
  Ns <- c(250L, 300L, 320L, 321L, 330L)
  cv <- do.call(bam_assurance_curve, c(a, list(N_range = 250:330)))
  expect_identical(cv$curve$assurance[match(Ns, cv$curve$N)], A_ref(Ns, a))
  expect_identical(cv$width_type, "lower_half_width")
})

test_that("identical for assurance_for = 'se' and 'sp'", {
  skip_on_cran()
  for (af in c("se", "sp")) {
    a <- c(p1, list(assurance_for = af))
    Ns <- c(2L, 5L, 60L, 100L, 150L)
    cv <- do.call(bam_assurance_curve, c(a, list(N_range = 2:150)))
    expect_identical(cv$curve$assurance[match(Ns, cv$curve$N)], A_ref(Ns, a))
    expect_identical(cv$assurance_for, af)
  }
})

test_that("identical with distinct analysis priors", {
  skip_on_cran()
  a <- c(p1, list(analysis_prior_se = c(1, 1), analysis_prior_sp = c(1, 1)))
  Ns <- c(2L, 10L, 100L, 300L, 400L)
  cv <- do.call(bam_assurance_curve, c(a, list(N_range = 2:400)))
  expect_identical(cv$curve$assurance[match(Ns, cv$curve$N)], A_ref(Ns, a))
  expect_identical(cv$analysis_prior_se, c(1, 1))
  expect_true(all(cv$analysis_prior_given))
})

test_that("the priors object is handled like bam_sample_size()", {
  skip_on_cran()
  pr <- bam_prior_from_pilot(
    tp = 16, fn = 1, tn = 20, fp = 35,
    base_prior_se = c(9.9, 1.1), base_prior_sp = c(1, 1),
    base_prior_prev = c(12, 43), a0 = 0.5, use_in = "design")
  cv <- suppressWarnings(bam_assurance_curve(
    priors = pr, delta_se = 0.30, delta_sp = 0.30, N_range = c(40L, 80L, 120L)))
  ref <- suppressWarnings(bam_sample_size(
    priors = pr, delta_se = 0.30, delta_sp = 0.30, method = "exact",
    N_range = 80L, n_range = 20:30, B = 1000))$max_assurance_evaluated
  expect_identical(cv$curve$assurance[2], ref)
  expect_identical(cv$prior_se, pr$prior_se)
  expect_identical(cv$analysis_prior_se, pr$analysis_prior_se)
  expect_error(bam_assurance_curve(priors = pr, prior_se = c(2, 2),
                                   N_range = 50L), "do not also pass")
  expect_error(bam_assurance_curve(priors = list(), N_range = 50L),
               "dta_bam_priors")
})

test_that("P1 first crossing N = 678 with the documented values", {
  skip_on_cran()
  cv <- do.call(bam_assurance_curve, c(p1, list(N_range = 670:690)))
  expect_identical(cv$N_first, 678L)
  # The reference values are rounded to 10 decimals: absolute tolerance
  # 5e-11 (half a unit of the last digit given).
  a678 <- cv$curve$assurance[cv$curve$N == 678]
  a677 <- cv$curve$assurance[cv$curve$N == 677]
  expect_lt(abs(a678 - 0.8003489948), 5e-11)
  expect_lt(abs(a677 - 0.7996848824), 5e-11)
  expect_identical(cv$target_assurance, 0.80)
})

test_that("Williamson et al. values", {
  skip_on_cran()
  cv <- bam_assurance_curve(
    prior_se = c(25, 7), prior_sp = c(30, 2), prior_prev = c(13.56, 122.06),
    width_type = "lower_half_width", delta_se = 0.10, delta_sp = 0.05,
    N_range = 320:322)
  expect_equal(cv$curve$assurance, c(0.7985633, 0.8007552, 0.8029264),
               tolerance = 5e-8)
  expect_identical(cv$N_first, 321L)
})

test_that("structure, empty-arm convention, N_first NA and monotonicity", {
  skip_on_cran()
  cv <- bam_assurance_curve(prior_se = c(17, 3), prior_sp = c(2, 2),
                            prior_prev = c(4, 16),
                            N_range = c(3L, 2L, 3L, 20L, 10L))
  expect_identical(cv$curve$N, c(2L, 3L, 10L, 20L))   # sorted, unique
  expect_named(cv$curve, c("N", "assurance"))
  expect_true(all(cv$curve$assurance >= 0 & cv$curve$assurance <= 1))
  expect_identical(cv$N_first, NA_integer_)           # nothing reaches 0.80
  expect_false(cv$target_reached)
  expect_true(is.numeric(cv$n_decreases) && is.numeric(cv$max_decrease))
  # N = 1 is degenerate (one empty arm) and not accepted: N_range >= 2
  expect_error(bam_assurance_curve(N_range = c(1L, 2L)), "N_range")
  # monotonicity bookkeeping is consistent with the curve
  d <- diff(cv$curve$assurance)
  expect_identical(cv$n_decreases, sum(d < 0))
  expect_identical(cv$max_decrease, max(0, -d))
  expect_identical(cv$monotone, cv$n_decreases == 0L)
})

test_that("m9: the empty-arm convention is applied by the curve (A(N) excludes the empty-arm term)", {
  skip_on_cran()
  # Beta(170, 30) already meets the width with no data; at small N the arm is
  # empty with appreciable probability, so crediting it would change A(N).
  prior <- c(170, 30); delta <- 0.14; pp <- c(2, 18)
  expect_lte(stats::qbeta(0.975, 170, 30) - stats::qbeta(0.025, 170, 30), delta)
  cv <- suppressWarnings(bam_assurance_curve(
    prior_se = prior, prior_prev = pp, delta_se = delta, assurance_for = "se",
    N_range = c(3L, 5L, 8L)))
  direct <- function(N, credit) {
    Parm <- function(n) {
      if (n == 0L) return(if (credit) 1 else 0)
      x <- 0:n
      pa <- prior[1] + x; pb <- prior[2] + n - x
      w <- stats::qbeta(0.975, pa, pb) - stats::qbeta(0.025, pa, pb)
      sum(exp(lchoose(n, x) + lbeta(prior[1] + x, prior[2] + n - x) -
                lbeta(prior[1], prior[2]))[w <= delta])
    }
    k <- 0:N
    pk <- exp(lchoose(N, k) + lbeta(pp[1] + k, pp[2] + N - k) - lbeta(pp[1], pp[2]))
    sum(pk * vapply(k, Parm, numeric(1)))
  }
  for (i in seq_along(cv$curve$N)) {
    N <- cv$curve$N[i]
    expect_lte(abs(cv$curve$assurance[i] - direct(N, FALSE)), 1e-12)
    expect_gt(direct(N, TRUE) - cv$curve$assurance[i], 0.1)
  }
})

test_that("print method names the conventions", {
  skip_on_cran()
  cv <- bam_assurance_curve(N_range = c(50L, 100L))
  out <- paste(capture.output(print(cv)), collapse = "\n")
  expect_match(out, "full width", ignore.case = TRUE)
  expect_match(out, "joint", ignore.case = TRUE)
  expect_match(out, "design", ignore.case = TRUE)
  expect_match(out, "analysis", ignore.case = TRUE)
  invisible(capture.output(r <- withVisible(print(cv))))
  expect_false(r$visible)
  expect_identical(r$value, cv)
})

test_that("input errors", {
  expect_error(bam_assurance_curve(), "N_range")
  expect_error(bam_assurance_curve(N_range = 1:3, prior_se = c(1, 2, 3)))
  expect_error(bam_assurance_curve(N_range = 10L, delta_se = -1))
  expect_error(bam_assurance_curve(N_range = 10L, target_assurance = 1))
  expect_error(bam_assurance_curve(N_range = 10L, width_type = "half"),
               "width_type")
  expect_error(bam_assurance_curve(N_range = 10L, assurance_for = "x"),
               "assurance_for")
  expect_error(bam_assurance_curve(N_range = 10L, analysis_prior_se = c(0, 1)),
               "analysis_prior_se")
  expect_error(bam_assurance_curve(N_range = c(0L, 5L)), "N_range")
  expect_error(bam_assurance_curve(N_range = c(5.5, 6)), "N_range")
  expect_error(bam_assurance_curve(N_range = c(NA, 6L)), "N_range")
  expect_error(bam_assurance_curve(N_range = integer(0)), "N_range")
  expect_error(bam_assurance_curve(N_range = 10L, priors = 3), "dta_bam_priors")
  expect_error(bam_assurance_curve(N_range = 10L, alpha_ci = 1), "alpha_ci")
})
