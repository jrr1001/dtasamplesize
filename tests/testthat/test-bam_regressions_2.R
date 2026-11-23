# Regression tests for the bam_* functions of version 0.7.0 (part 2).
# See NEWS.md (0.7.0) and the help pages of bam_prior_contour(),
# bam_sensitivity(), bam_sensitivity_local() and bam_sample_size().

## ---------------------------------------------------------------------
## deterministic precision diagnostic of bam_prior_contour()
## ---------------------------------------------------------------------

test_that("Beta(1000, 100), H2, epsilon 1e-8: every direction is flagged, with a warning", {
  expect_warning(
    ct <- bam_prior_contour(c(1000, 100), epsilon = 1e-8, distance = "H2",
                            n_phi = 24L),
    "precision")
  expect_true(all(attr(ct, "imprecise")))
  pb <- attr(ct, "precision_bound")
  expect_true(is.numeric(pb))
  expect_length(pb, 24L)
  expect_true(all(pb > dtasamplesize:::.bam_contour_tol))
})

test_that("Beta(1000, 100), H, epsilon 1e-4: the flags do not depend on n_phi", {
  res <- lapply(c(23L, 24L, 25L), function(n) {
    ct <- suppressWarnings(bam_prior_contour(c(1000, 100), epsilon = 1e-4,
                                             distance = "H", n_phi = n))
    list(imprecise = attr(ct, "imprecise"), pb = attr(ct, "precision_bound"))
  })
  for (r in res) expect_true(all(r$imprecise))
  # the (relative) bound itself is the same order of magnitude in all three
  rng <- vapply(res, function(r) range(r$pb), numeric(2))
  expect_lt(max(rng) / min(rng), 2)
})

test_that("the VAP contour is precise: nothing flagged, no warning, bound < 1e-9", {
  expect_warning(
    ct <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354, n_phi = 36L),
    NA)
  expect_false(any(attr(ct, "imprecise")))
  pb <- attr(ct, "precision_bound")
  expect_length(pb, 36L)
  expect_true(all(pb < 1e-9))
})

test_that("Beta(1e5, 1e4), H, epsilon 1e-6 (the example of the Rd) is flagged", {
  expect_warning(
    ct <- bam_prior_contour(c(1e5, 1e4), epsilon = 1e-6, distance = "H",
                            n_phi = 12L),
    "precision")
  expect_true(all(attr(ct, "imprecise")))
})

test_that("`imprecise` is the stated rule applied to precision_bound and to the recomputed distance", {
  cases <- list(
    list(c(9.9, 1.1), 0.00354, "H2"), list(c(9.9, 1.1), 0.00354, "H"),
    list(c(1000, 100), 1e-4, "H"), list(c(1000, 100), 1e-8, "H2"),
    list(c(300, 40), 1e-5, "H2"), list(c(300, 40), 1e-3, "H"),
    list(c(0.5, 0.5), 0.5, "H"), list(c(2, 2), 0.95, "H2"))
  tol <- dtasamplesize:::.bam_contour_tol
  for (cs in cases) {
    ct <- suppressWarnings(bam_prior_contour(cs[[1]], epsilon = cs[[2]],
                                             distance = cs[[3]], n_phi = 16L))
    imp <- attr(ct, "imprecise")
    pb <- attr(ct, "precision_bound")
    has <- !is.na(ct$a)
    expect_identical(is.na(pb), !has)
    rule <- !is.finite(ct$distance) | ct$distance == 0 |
      abs(ct$distance - cs[[2]]) > tol * cs[[2]] | is.na(pb) | pb > tol
    expect_identical(imp[has], rule[has])
    expect_false(any(imp[!has]))
  }
})

test_that("the bound is the stated rounding-error estimate (independent computation)", {
  a0 <- 1000; b0 <- 100; eps <- 1e-4
  ct <- suppressWarnings(bam_prior_contour(c(a0, b0), epsilon = eps,
                                           distance = "H", n_phi = 8L))
  a <- ct$a; b <- ct$b
  err <- 4 * .Machine$double.eps *
    ((4 + abs(lbeta((a + a0) / 2, (b + b0) / 2))) + (4 + abs(lbeta(a, b))) / 2 +
       (4 + abs(lbeta(a0, b0))) / 2)
  d2 <- ct$distance^2
  expect_equal(attr(ct, "precision_bound"), err / (2 * d2), tolerance = 1e-12)
  ct2 <- suppressWarnings(bam_prior_contour(c(a0, b0), epsilon = 1e-8,
                                            distance = "H2", n_phi = 8L))
  err2 <- 4 * .Machine$double.eps *
    ((4 + abs(lbeta((ct2$a + a0) / 2, (ct2$b + b0) / 2))) +
       (4 + abs(lbeta(ct2$a, ct2$b))) / 2 + (4 + abs(lbeta(a0, b0))) / 2)
  expect_equal(attr(ct2, "precision_bound"), err2 / ct2$distance,
               tolerance = 1e-12)
})

test_that("the contour points themselves are unchanged (phi, r, a, b, distance)", {
  ct <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354, n_phi = 36L)
  # the distance at every point is epsilon to within the root tolerance
  expect_lt(max(abs(ct$distance - 0.00354)) / 0.00354, 1e-9)
  expect_identical(names(ct), c("phi", "r", "a", "b", "distance"))
})

test_that("bam_sensitivity_local() counts the directions discarded for precision (partial case, mocked contour)", {
  skip_on_cran()
  real <- bam_prior_contour
  testthat::local_mocked_bindings(
    bam_prior_contour = function(...) {
      ct <- real(...)
      imp <- attr(ct, "imprecise")
      imp[c(2L, 5L)] <- TRUE
      attr(ct, "imprecise") <- imp
      ct
    })
  res <- suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(9.9, 1.1), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
    n_phi = 6L, n_ref = 106L, delta_se = 0.16, assurance_for = "se"))
  expect_identical(res$summary$n_dir_imprecise, 2L)
  expect_identical(res$summary$n_dir_no_root, 0L)
  expect_identical(nrow(res$table), 4L)
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "(0 without a root, 2 with insufficient numerical precision)",
               fixed = TRUE)
})

test_that("all directions imprecise is an error that says why", {
  skip_on_cran()
  expect_error(suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(1000, 100), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 1e-8,
    distance = "H2", n_phi = 6L, n_ref = 50L, delta_se = 0.30,
    assurance_for = "se")), "numerical precision")
})

## ---------------------------------------------------------------------
## one summary warning in bam_sensitivity() and bam_sensitivity_local()
## ---------------------------------------------------------------------

collect_warnings <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  })
  list(value = v, warnings = w)
}

test_that("bam_sensitivity_local() with a prior that already meets the width gives exactly one warning", {
  skip_on_cran()
  r <- collect_warnings(bam_sensitivity_local(
    base_prior_se = c(60, 6), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
    n_phi = 12, delta_se = 0.16, assurance_for = "se", n_ref = 50))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "analysis prior", ignore.case = TRUE)
  expect_match(r$warnings, "contour point")
  expect_match(r$warnings, "bam_sample_size")
  expect_match(r$warnings, "base prior")
  # same without n_ref (the base N is found first)
  r2 <- collect_warnings(bam_sensitivity_local(
    base_prior_se = c(60, 6), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
    n_phi = 6, delta_se = 0.16, assurance_for = "se"))
  expect_length(r2$warnings, 1L)
})

test_that("bam_sensitivity_local() with no affected point gives no warning", {
  skip_on_cran()
  r <- collect_warnings(bam_sensitivity_local(
    base_prior_se = c(9.9, 1.1), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
    n_phi = 4, delta_se = 0.16, assurance_for = "se", n_ref = 106))
  expect_length(r$warnings, 0L)
})

scen_f2 <- list(
  narrow_prior = list(prior_se = c(60, 6)),
  weak_prior = list(prior_se = c(8.5, 1.5))
)
args_f2 <- list(prior_sp = c(2, 2), prior_prev = c(4, 16), delta_se = 0.16,
                assurance_for = "se")

test_that("bam_sensitivity() gives one warning that names only the affected scenario", {
  skip_on_cran()
  r <- collect_warnings(do.call(bam_sensitivity, c(list(scen_f2), args_f2)))
  expect_length(r$warnings, 1L)
  expect_match(r$warnings, "narrow_prior", fixed = TRUE)
  expect_false(grepl("weak_prior", r$warnings, fixed = TRUE))
  expect_match(r$warnings, "bam_sample_size")
  expect_identical(r$value$scenario, c("narrow_prior", "weak_prior"))
})

test_that("warnings and results are the same with parallel = TRUE and FALSE", {
  skip_on_cran()
  rs <- collect_warnings(do.call(bam_sensitivity, c(list(scen_f2), args_f2)))
  rp <- collect_warnings(do.call(bam_sensitivity,
                                 c(list(scen_f2), args_f2,
                                   list(parallel = TRUE, cores = 2L))))
  expect_length(rp$warnings, 1L)
  expect_identical(rp$warnings, rs$warnings)
  expect_identical(rp$value, rs$value)

  la <- list(base_prior_se = c(60, 6), base_prior_sp = c(2, 2),
             base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
             n_phi = 4, delta_se = 0.16, assurance_for = "se", n_ref = 50)
  ls_ <- collect_warnings(do.call(bam_sensitivity_local, la))
  lp_ <- collect_warnings(do.call(bam_sensitivity_local,
                                  c(la, list(parallel = TRUE, cores = 2L))))
  expect_length(lp_$warnings, 1L)
  expect_identical(lp_$warnings, ls_$warnings)
  expect_identical(lp_$value$table, ls_$value$table)
  expect_identical(lp_$value$summary, ls_$value$summary)
})

test_that("the example of ?bam_sensitivity gives no warning", {
  skip_on_cran()
  scen <- list(
    informative = list(prior_se = c(17, 3), prior_sp = c(18, 2)),
    weaker      = list(prior_se = c(8.5, 1.5), prior_sp = c(9, 1)),
    se_only     = list(prior_se = c(17, 3), prior_sp = c(18, 2),
                       assurance_for = "se")
  )
  r <- collect_warnings(bam_sensitivity(scen, prior_prev = c(4, 16),
                                        delta_se = 0.22, delta_sp = 0.22))
  expect_length(r$warnings, 0L)
  df <- data.frame(
    scenario = c("informative", "weaker"),
    prior_se_a = c(17, 8.5), prior_se_b = c(3, 1.5),
    prior_sp_a = c(18, 9), prior_sp_b = c(2, 1)
  )
  r2 <- collect_warnings(bam_sensitivity(df, prior_prev = c(4, 16),
                                         delta_se = 0.22, delta_sp = 0.22))
  expect_length(r2$warnings, 0L)
})

## ---------------------------------------------------------------------
## empty-arm convention in the default exact scan (N_range = NULL)
## ---------------------------------------------------------------------

# Probability that an arm of size n meets the width, under the prior used as
# design and analysis prior; credit_empty = TRUE would credit the n = 0 term.
arm_ok <- function(n, prior, delta, credit_empty) {
  if (n == 0L) {
    w0 <- stats::qbeta(0.975, prior[1], prior[2]) -
      stats::qbeta(0.025, prior[1], prior[2])
    return(if (credit_empty) as.numeric(w0 <= delta) else 0)
  }
  x <- 0:n
  w <- stats::qbeta(0.975, prior[1] + x, prior[2] + n - x) -
    stats::qbeta(0.025, prior[1] + x, prior[2] + n - x)
  sum(exp(lchoose(n, x) + lbeta(prior[1] + x, prior[2] + n - x) -
            lbeta(prior[1], prior[2]))[w <= delta])
}
enum_A <- function(N, which, prior_se, prior_sp, delta_se, delta_sp, prev,
                   credit_empty) {
  k <- 0:N
  pk <- exp(lchoose(N, k) + lbeta(prev[1] + k, prev[2] + N - k) -
              lbeta(prev[1], prev[2]))
  se <- if (which == "sp") 1 else
    vapply(k, arm_ok, numeric(1), prior = prior_se, delta = delta_se,
           credit_empty = credit_empty)
  sp <- if (which == "se") 1 else
    vapply(N - k, arm_ok, numeric(1), prior = prior_sp, delta = delta_sp,
           credit_empty = credit_empty)
  sum(pk * se * sp)
}
enum_first_N <- function(..., target = 0.80, N_max = 60L, credit_empty) {
  for (N in 2:N_max) {
    A <- enum_A(N, ..., credit_empty = credit_empty)
    if (A >= target) return(c(N = N, A = A))
  }
  c(N = NA, A = NA)
}

test_that("default exact scan (N_range = NULL), marginal Se: the empty-arm term is excluded", {
  p <- c(170, 30)
  r <- suppressWarnings(bam_sample_size(
    prior_se = p, prior_prev = c(2, 18), delta_se = 0.14,
    assurance_for = "se", method = "exact", B = 1000, n_range = 20L))
  ref <- enum_first_N(which = "se", prior_se = p, prior_sp = p, delta_se = 0.14,
                      delta_sp = 0.14, prev = c(2, 18), credit_empty = FALSE)
  bad <- enum_first_N(which = "se", prior_se = p, prior_sp = p, delta_se = 0.14,
                      delta_sp = 0.14, prev = c(2, 18), credit_empty = TRUE)
  expect_identical(r$N_total, as.integer(ref[["N"]]))
  expect_lte(abs(r$joint_assurance - ref[["A"]]), 1e-12)
  # the convention matters here: including the empty arm gives another answer
  expect_false(identical(as.integer(bad[["N"]]), r$N_total))
  expect_gt(r$N_total, 2L)
  expect_lt(bad[["N"]], r$N_total)
})

test_that("default exact scan (N_range = NULL), marginal Sp (mirror)", {
  p <- c(170, 30)
  r <- suppressWarnings(bam_sample_size(
    prior_sp = p, prior_prev = c(18, 2), delta_sp = 0.14,
    assurance_for = "sp", method = "exact", B = 1000, n_range = 20L))
  ref <- enum_first_N(which = "sp", prior_se = p, prior_sp = p, delta_se = 0.14,
                      delta_sp = 0.14, prev = c(18, 2), credit_empty = FALSE)
  bad <- enum_first_N(which = "sp", prior_se = p, prior_sp = p, delta_se = 0.14,
                      delta_sp = 0.14, prev = c(18, 2), credit_empty = TRUE)
  expect_identical(r$N_total, as.integer(ref[["N"]]))
  expect_lte(abs(r$joint_assurance - ref[["A"]]), 1e-12)
  expect_false(identical(as.integer(bad[["N"]]), r$N_total))
  expect_gt(r$N_total, 2L)
})

test_that("default exact scan (N_range = NULL), joint objective", {
  p <- c(170, 30)
  r <- suppressWarnings(bam_sample_size(
    prior_se = p, prior_sp = p, prior_prev = c(2, 18), delta_se = 0.14,
    delta_sp = 0.14, assurance_for = "joint", method = "exact", B = 1000,
    n_range = 20L))
  ref <- enum_first_N(which = "joint", prior_se = p, prior_sp = p,
                      delta_se = 0.14, delta_sp = 0.14, prev = c(2, 18),
                      credit_empty = FALSE)
  bad <- enum_first_N(which = "joint", prior_se = p, prior_sp = p,
                      delta_se = 0.14, delta_sp = 0.14, prev = c(2, 18),
                      credit_empty = TRUE)
  expect_identical(r$N_total, as.integer(ref[["N"]]))
  expect_lte(abs(r$joint_assurance - ref[["A"]]), 1e-12)
  expect_false(identical(as.integer(bad[["N"]]), r$N_total))
  expect_gt(r$N_total, 2L)
})

## ---------------------------------------------------------------------
## wording of the prior-only warning is decided by the values
## ---------------------------------------------------------------------

test_that("an explicit analysis prior equal to the design prior says 'equal to the design prior'", {
  p <- c(170, 30)
  w <- character()
  withCallingHandlers(
    bam_sample_size(prior_se = p, analysis_prior_se = p, prior_prev = c(2, 18),
                    delta_se = 0.14, assurance_for = "se", method = "exact",
                    N_range = 30L, B = 1000, n_range = 20L),
    warning = function(c) { w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning") })
  w <- w[grepl("^The analysis prior for", w)]
  expect_length(w, 1L)
  expect_match(w, "equal to the design prior")
})

test_that("an analysis prior that differs from the design prior does not say so", {
  w <- character()
  withCallingHandlers(
    bam_sample_size(prior_se = c(170, 30), analysis_prior_se = c(160, 28),
                    prior_prev = c(2, 18), delta_se = 0.14,
                    assurance_for = "se", method = "exact", N_range = 30L,
                    B = 1000, n_range = 20L),
    warning = function(c) { w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning") })
  w <- w[grepl("^The analysis prior for", w)]
  expect_length(w, 1L)
  expect_false(grepl("equal to the design prior", w))
})

test_that("same wording by value in bam_assurance_curve()", {
  p <- c(170, 30)
  w <- character()
  withCallingHandlers(
    bam_assurance_curve(prior_se = p, analysis_prior_se = p,
                        prior_prev = c(2, 18), delta_se = 0.14,
                        assurance_for = "se", N_range = c(20L, 30L)),
    warning = function(c) { w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning") })
  w <- w[grepl("^The analysis prior for", w)]
  expect_length(w, 1L)
  expect_match(w, "equal to the design prior")
})

## ---------------------------------------------------------------------
## print of the marginal objective and of the discarded directions
## ---------------------------------------------------------------------

test_that("target not reached with a marginal objective names the marginal assurance", {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_prev = c(4, 16), delta_se = 0.001,
    assurance_for = "se", method = "exact", N_max = 40L, B = 1000,
    n_range = 20L))
  expect_false(r$target_reached)
  out <- paste(capture.output(print(r)), collapse = "\n")
  expect_false(grepl("joint_assurance are NA", out, fixed = TRUE))
  expect_match(out, "marginal Se assurance are NA", fixed = TRUE)
  expect_match(out, "max_assurance_evaluated (marginal Se assurance):",
               fixed = TRUE)
  r2 <- suppressWarnings(bam_sample_size(
    prior_sp = c(17, 3), prior_prev = c(4, 16), delta_sp = 0.001,
    assurance_for = "sp", method = "exact", N_max = 40L, B = 1000,
    n_range = 20L))
  out2 <- paste(capture.output(print(r2)), collapse = "\n")
  expect_match(out2, "marginal Sp assurance are NA", fixed = TRUE)
  expect_match(out2, "max_assurance_evaluated (marginal Sp assurance):",
               fixed = TRUE)
})

test_that("target not reached with the joint objective is unchanged", {
  r <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
    delta_se = 0.001, delta_sp = 0.001, method = "exact", N_max = 40L,
    B = 1000, n_range = 20L))
  expect_false(r$target_reached)
  out <- capture.output(print(r))
  expect_true(any(grepl("(N_total / n_total / joint_assurance are NA; see the warning",
                        out, fixed = TRUE)))
  expect_true(any(grepl("^  max_assurance_evaluated: ", out)))
})

test_that("'Directions discarded' reads (0 without a root, 1 with ...)", {
  skip_on_cran()
  real <- bam_prior_contour
  testthat::local_mocked_bindings(
    bam_prior_contour = function(...) {
      ct <- real(...)
      imp <- attr(ct, "imprecise")
      imp[3L] <- TRUE
      attr(ct, "imprecise") <- imp
      ct
    })
  res <- suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(9.9, 1.1), base_prior_sp = c(2, 2),
    base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
    n_phi = 6L, n_ref = 106L, delta_se = 0.16, assurance_for = "se"))
  out <- capture.output(print(res))
  expect_true(any(grepl(
    "Directions discarded from the minima and maxima: 1 (0 without a root, 1 with insufficient numerical precision)",
    out, fixed = TRUE)))
})

