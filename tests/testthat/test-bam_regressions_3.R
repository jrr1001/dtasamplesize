# Regression tests for the bam_* functions of version 0.7.0 (part 3).
# See NEWS.md (0.7.0) and the help pages of bam_prior_contour() and
# bam_sensitivity().

# Reference for the distance H2 = 1 - BC at a contour point WITHOUT
# cancellation: Taylor series of log BC in the differences of the parameters
# (terms up to psigamma(., 5)).
.ref_h2_taylor <- function(a, b, a0, b0) {
  da <- (a - a0) / 2; db <- (b - b0) / 2
  am <- a0 + da; bm <- b0 + db
  q <- function(k) psigamma(am, k) * da^(k + 1) + psigamma(bm, k) * db^(k + 1) -
    psigamma(am + bm, k) * (da + db)^(k + 1)
  lbc <- -(q(1) / 2 + q(3) / 24 + q(5) / 720)
  -expm1(lbc)
}

test_that("Beta(1, 1), phi = 0, H2: the bound covers the real relative error", {
  # On phi = 0, b stays at 1 and a = 1 + r; H2 has a closed form without
  # cancellation: r^2 / ((sqrt(1 + r) + 1)^2 (2 + r)).
  h2_ref <- function(r) r^2 / ((sqrt(1 + r) + 1)^2 * (2 + r))
  for (eps in c(1e-6, 1e-8, 1e-9, 1e-10, 1e-11)) {
    ct <- suppressWarnings(bam_prior_contour(c(1, 1), epsilon = eps,
                                             distance = "H2", phi = 0))
    pb <- attr(ct, "precision_bound")
    imp <- attr(ct, "imprecise")
    err <- abs(h2_ref(ct$r) - eps) / eps
    expect_lte(err, pb)
    if (!imp) expect_lte(err, 1e-6)
    if (eps == 1e-11) expect_true(imp)
  }
})

test_that("Taylor reference over priors, epsilons and 24 directions: bound >= error", {
  for (p in list(c(1, 1), c(0.5, 3.6), c(9.9, 1.1), c(60, 6))) {
    for (eps in c(1e-8, 1e-9)) {
      ct <- suppressWarnings(bam_prior_contour(p, epsilon = eps,
                                               distance = "H2", n_phi = 24L))
      ok <- !is.na(ct$a)
      td <- .ref_h2_taylor(ct$a[ok], ct$b[ok], p[1], p[2])
      err <- abs(ct$distance[ok] - td) / td
      pb <- attr(ct, "precision_bound")[ok]
      imp <- attr(ct, "imprecise")[ok]
      info <- paste0("prior ", paste(p, collapse = ","), ", epsilon ", eps)
      expect_lt(max(err / pb), 1, label = paste("max error/bound,", info))
      expect_true(all(err[!imp] <= 1e-6), info = info)
    }
  }
})

test_that("the bound uses the absolute floor of 4 per lbeta() term", {
  p <- c(1, 1); eps <- 1e-10
  ct <- suppressWarnings(bam_prior_contour(p, epsilon = eps, distance = "H2",
                                           n_phi = 12L))
  l1 <- lbeta((ct$a + p[1]) / 2, (ct$b + p[2]) / 2)
  l2 <- lbeta(ct$a, ct$b); l3 <- lbeta(p[1], p[2])
  err <- 4 * .Machine$double.eps *
    ((4 + abs(l1)) + (4 + abs(l2)) / 2 + (4 + abs(l3)) / 2)
  expect_equal(attr(ct, "precision_bound"), err / ct$distance, tolerance = 1e-12)
  # distance "H": half the relative bound, on the distance H^2 = d^2
  cth <- suppressWarnings(bam_prior_contour(p, epsilon = 1e-5, distance = "H",
                                            n_phi = 12L))
  l1 <- lbeta((cth$a + p[1]) / 2, (cth$b + p[2]) / 2)
  l2 <- lbeta(cth$a, cth$b)
  err <- 4 * .Machine$double.eps *
    ((4 + abs(l1)) + (4 + abs(l2)) / 2 + (4 + abs(l3)) / 2)
  expect_equal(attr(cth, "precision_bound"), err / (2 * cth$distance^2),
               tolerance = 1e-12)
})

## ---------------------------------------------------------------------
## the summary warning is emitted before the error of a scenario
## ---------------------------------------------------------------------

test_that("bam_sensitivity() with a failing scenario keeps the prior-only warning", {
  w <- character()
  e <- tryCatch(
    withCallingHandlers(
      bam_sensitivity(
        list(ok = list(prior_se = c(60, 6)),
             bad = list(prior_se = c(60, 6), N_max = 1)),
        prior_sp = c(2, 2), prior_prev = c(4, 16), delta_se = 0.16,
        assurance_for = "se"),
      warning = function(c) {
        w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
      }),
    error = function(e) conditionMessage(e))
  expect_match(e, "bad")
  expect_true(any(grepl("'ok'", w, fixed = TRUE)))
})

test_that("bam_sensitivity_local() with a failing contour point keeps the prior-only warning (mocked point)", {
  skip_on_cran()
  real <- dtasamplesize:::.bam_local_point
  k <- 0L
  testthat::local_mocked_bindings(
    .bam_local_point = function(...) {
      r <- real(...)
      k <<- k + 1L
      if (k == 2L) r$error <- "mocked failure"
      r
    })
  w <- character()
  e <- tryCatch(
    withCallingHandlers(
      bam_sensitivity_local(
        base_prior_se = c(60, 6), base_prior_sp = c(2, 2),
        base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
        n_phi = 6L, n_ref = 50L, delta_se = 0.16, assurance_for = "se"),
      warning = function(c) {
        w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
      }),
    error = function(e) conditionMessage(e))
  expect_match(e, "mocked failure")
  expect_true(any(grepl("analysis prior", w)))
})
