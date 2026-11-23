# width_type and assurance_for of bam_sample_size().
# Reference values: Williamson et al. 2023 (one-sided), Wilson et al. 2022
# (VAP, marginal objective for Se) and Figure 3 of Williamson et al.

quick_w <- function(...) {
  suppressWarnings(bam_sample_size(n_range = 20:60, B = 1000, seed = 3, ...))
}

## ---------------------------------------------------------------------
## One-sided criterion: example of Williamson et al. 2023
## ---------------------------------------------------------------------

will_args <- list(
  prior_se = c(25, 7), prior_sp = c(30, 2), prior_prev = c(13.56, 122.06),
  width_type = "lower_half_width", delta_se = 0.10, delta_sp = 0.05,
  target_assurance = 0.80, method = "exact"
)

test_that("Williamson et al. 2023 one-sided example gives N = 321, A = 0.8007552", {
  res <- do.call(quick_w, will_args)
  message("Williamson: N = ", res$N_total, ", A = ",
          sprintf("%.10f", res$joint_assurance))
  expect_identical(res$N_total, 321L)
  expect_lt(abs(res$joint_assurance - 0.8007552), 5e-8)
  expect_identical(res$width_type, "lower_half_width")
})

test_that("Williamson A(320) = 0.7985633 and A(322) = 0.8029264", {
  r320 <- do.call(quick_w, c(will_args, list(N_range = 320)))
  r322 <- do.call(quick_w, c(will_args, list(N_range = 322)))
  expect_lt(abs(r320$max_assurance_evaluated - 0.7985633), 5e-8)
  expect_lt(abs(r322$max_assurance_evaluated - 0.8029264), 5e-8)
  expect_false(r320$target_reached)
  expect_true(r322$target_reached)
})

test_that("default width_type is full_width and does not change the default result", {
  res <- quick_w(prior_prev = c(4, 16), N_range = 600:700)
  expect_identical(res$width_type, "full_width")
  expect_identical(res$N_total, 678L)
  res2 <- quick_w(prior_prev = c(4, 16), N_range = 600:700, width_type = "full_width")
  expect_identical(res2$joint_assurance, res$joint_assurance)
})

test_that("the one-sided criterion is median minus the one-sided lower bound (direct check, Monte Carlo)", {
  # With fixed N and a single n, the exact mode coincides with the direct sum
  a <- 25; b <- 7; delta <- 0.10; n <- 60
  x <- 0:n
  w <- stats::qbeta(0.5, a + x, b + n - x) - stats::qbeta(0.05, a + x, b + n - x)
  P <- sum(exp(lchoose(n, x) + lbeta(a + x, b + n - x) - lbeta(a, b))[w <= delta])
  Pfun <- dtasamplesize:::.bam_exact_width_prob(
    n, a, b, delta, 0.025, 0.975, width_type = "lower_half_width"
  )
  expect_identical(Pfun[n + 1L], P)
})

test_that("Buderer comparison uses the half-width implied by width_type", {
  full <- quick_w(prior_se = c(25, 7), delta_se = 0.10, N_range = 100:150)
  half <- quick_w(prior_se = c(25, 7), delta_se = 0.10, N_range = 100:150,
                  width_type = "lower_half_width")
  p <- 25 / 32
  expect_identical(full$buderer_n_se, buderer_n(p, 0.10 / 2))
  # one-sided criterion: half-width delta_se at the one-sided z = qnorm(alpha_ci)
  expect_identical(half$buderer_n_se, buderer_n(p, 0.10, alpha = 2 * (1 - 0.95)))
  expect_equal(half$buderer_n_se, buderer_n(p, 0.10, alpha = 0.10))
})

test_that("lower_half_width also drives the Monte Carlo search and the legacy fields", {
  skip_on_cran()
  args <- c(will_args, list(n_range = 20:300, B = 20000, seed = 2026))
  ex <- suppressWarnings(do.call(bam_sample_size, c(args, list(N_range = 321))))
  mc <- suppressWarnings(do.call(bam_sample_size, c(
    args[setdiff(names(args), "method")],
    list(method = "monte_carlo", N_range = 321)
  )))
  A <- ex$max_assurance_evaluated
  mcse <- sqrt(A * (1 - A) / 20000)
  expect_lte(abs(mc$max_assurance_evaluated - A), 3 * mcse)
  # the per-arm legacy field also uses the one-sided criterion: with the
  # same prior, the one-sided n_diseased is smaller than the two-sided one
  full <- suppressWarnings(bam_sample_size(
    prior_se = c(25, 7), delta_se = 0.10, n_range = 20:300, B = 3000, N_range = 321
  ))
  half <- suppressWarnings(bam_sample_size(
    prior_se = c(25, 7), delta_se = 0.10, n_range = 20:300, B = 3000, N_range = 321,
    width_type = "lower_half_width"
  ))
  expect_lt(half$n_diseased, full$n_diseased)
})

## ---------------------------------------------------------------------
## Marginal objective
## ---------------------------------------------------------------------

# Direct marginal sum, written without the package code
marginal_direct <- function(N, which, prior_d, prior_a, delta, prior_prev,
                            width = c("full_width", "lower_half_width")) {
  width <- match.arg(width)
  Parm <- function(n) {
    if (n == 0L) return(0)
    x <- 0:n
    pa <- prior_a[1] + x; pb <- prior_a[2] + n - x
    w <- if (width == "full_width") stats::qbeta(0.975, pa, pb) - stats::qbeta(0.025, pa, pb)
         else stats::qbeta(0.5, pa, pb) - stats::qbeta(0.05, pa, pb)
    sum(exp(lchoose(n, x) + lbeta(prior_d[1] + x, prior_d[2] + n - x) -
              lbeta(prior_d[1], prior_d[2]))[w <= delta])
  }
  k <- 0:N
  pk <- exp(lchoose(N, k) + lbeta(prior_prev[1] + k, prior_prev[2] + N - k) -
              lbeta(prior_prev[1], prior_prev[2]))
  arm <- if (which == "se") vapply(k, Parm, numeric(1)) else
    vapply(N - k, Parm, numeric(1))
  sum(pk * arm)
}

test_that("assurance_for = 'se' equals the direct marginal sum (<= 1e-12)", {
  for (N in c(60L, 150L)) {
    res <- quick_w(
      prior_se = c(17, 3), prior_prev = c(4, 16), delta_se = 0.14,
      assurance_for = "se", N_range = N
    )
    direct <- marginal_direct(N, "se", c(17, 3), c(17, 3), 0.14, c(4, 16))
    expect_lte(abs(res$max_assurance_evaluated - direct), 1e-12)
  }
})

test_that("assurance_for = 'sp' equals the direct marginal sum (<= 1e-12)", {
  for (N in c(60L, 150L)) {
    res <- quick_w(
      prior_sp = c(2, 2), prior_prev = c(4, 16), delta_sp = 0.10,
      assurance_for = "sp", N_range = N
    )
    direct <- marginal_direct(N, "sp", c(2, 2), c(2, 2), 0.10, c(4, 16))
    expect_lte(abs(res$max_assurance_evaluated - direct), 1e-12)
  }
})

test_that("marginal objective with separate design/analysis priors and one-sided width", {
  res <- quick_w(
    prior_se = c(17, 3), analysis_prior_se = c(1, 1), prior_prev = c(4, 16),
    delta_se = 0.10, width_type = "lower_half_width",
    assurance_for = "se", N_range = 120
  )
  direct <- marginal_direct(120L, "se", c(17, 3), c(1, 1), 0.10, c(4, 16),
                            width = "lower_half_width")
  expect_lte(abs(res$max_assurance_evaluated - direct), 1e-12)
})

test_that("marginal A(N) >= joint A(N) and 'joint' remains the default", {
  N <- 700L
  joint <- quick_w(prior_prev = c(4, 16), N_range = N)
  se <- quick_w(prior_prev = c(4, 16), N_range = N, assurance_for = "se")
  sp <- quick_w(prior_prev = c(4, 16), N_range = N, assurance_for = "sp")
  expect_identical(joint$assurance_for, "joint")
  expect_identical(se$assurance_for, "se")
  expect_identical(sp$assurance_for, "sp")
  expect_gte(se$max_assurance_evaluated, joint$max_assurance_evaluated)
  expect_gte(sp$max_assurance_evaluated, joint$max_assurance_evaluated)
})

test_that("the unused arm's prior and delta do not matter for the headline result", {
  a <- quick_w(prior_se = c(17, 3), prior_prev = c(4, 16), delta_se = 0.14,
               assurance_for = "se", N_range = 100:300)
  b <- quick_w(prior_se = c(17, 3), prior_sp = c(40, 3), delta_sp = 0.001,
               prior_prev = c(4, 16), delta_se = 0.14,
               assurance_for = "se", N_range = 100:300)
  expect_identical(a$N_total, b$N_total)
  expect_identical(a$joint_assurance, b$joint_assurance)
})

test_that("Monte Carlo marginal agrees with exact (<= 3 MCSE)", {
  skip_on_cran()
  for (af in c("se", "sp")) {
    args <- list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
                 delta_se = 0.14, delta_sp = 0.10, assurance_for = af,
                 N_range = 100L, n_range = 20:20)
    ex <- suppressWarnings(do.call(bam_sample_size, c(args, list(method = "exact"))))
    mc <- suppressWarnings(do.call(bam_sample_size, c(args, list(
      method = "monte_carlo", B = 20000, seed = 2026
    ))))
    A <- ex$max_assurance_evaluated
    mcse <- sqrt(A * (1 - A) / 20000)
    expect_lte(abs(mc$max_assurance_evaluated - A), 3 * mcse)
  }
})

# Marginal assurance with the empty-arm term either credited (what would
# happen WITHOUT the convention: an empty arm keeps the prior-only interval)
# or set to 0 (the convention). `prior` is both design and analysis prior.
marginal_empty <- function(N, which, prior, delta, prior_prev, credit_empty) {
  w0 <- stats::qbeta(0.975, prior[1], prior[2]) -
    stats::qbeta(0.025, prior[1], prior[2])
  Parm <- function(n) {
    if (n == 0L) return(if (credit_empty) as.numeric(w0 <= delta) else 0)
    x <- 0:n
    pa <- prior[1] + x; pb <- prior[2] + n - x
    w <- stats::qbeta(0.975, pa, pb) - stats::qbeta(0.025, pa, pb)
    sum(exp(lchoose(n, x) + lbeta(prior[1] + x, prior[2] + n - x) -
              lbeta(prior[1], prior[2]))[w <= delta])
  }
  k <- 0:N
  pk <- exp(lchoose(N, k) + lbeta(prior_prev[1] + k, prior_prev[2] + N - k) -
              lbeta(prior_prev[1], prior_prev[2]))
  arm <- if (which == "se") vapply(k, Parm, numeric(1)) else
    vapply(N - k, Parm, numeric(1))
  sum(pk * arm)
}

test_that("m9: the empty-arm convention matters here, and the package applies it (exact)", {
  # Beta(170, 30): the prior-only 95% full width (~0.099) already meets 0.14,
  # and with N small the arm is empty with appreciable probability.
  prior <- c(170, 30)
  w0 <- stats::qbeta(0.975, 170, 30) - stats::qbeta(0.025, 170, 30)
  expect_lte(w0, 0.14)
  for (N in c(3L, 5L, 8L)) {
    # Se arm empty when n_d = 0 (prevalence prior mean 0.1)
    se <- suppressWarnings(bam_sample_size(
      prior_se = prior, prior_prev = c(2, 18), delta_se = 0.14,
      assurance_for = "se", N_range = N, n_range = 20L, B = 1000, seed = 3))
    with0 <- marginal_empty(N, "se", prior, 0.14, c(2, 18), credit_empty = TRUE)
    wout0 <- marginal_empty(N, "se", prior, 0.14, c(2, 18), credit_empty = FALSE)
    expect_gt(with0 - wout0, 0.1)                       # the convention matters
    expect_lte(abs(se$max_assurance_evaluated - wout0), 1e-12)
    expect_gt(abs(se$max_assurance_evaluated - with0), 0.1)
    # Sp arm empty when n_nd = 0 (prevalence prior mean 0.9)
    sp <- suppressWarnings(bam_sample_size(
      prior_sp = prior, prior_prev = c(18, 2), delta_sp = 0.14,
      assurance_for = "sp", N_range = N, n_range = 20L, B = 1000, seed = 3))
    with0 <- marginal_empty(N, "sp", prior, 0.14, c(18, 2), credit_empty = TRUE)
    wout0 <- marginal_empty(N, "sp", prior, 0.14, c(18, 2), credit_empty = FALSE)
    expect_gt(with0 - wout0, 0.1)
    expect_lte(abs(sp$max_assurance_evaluated - wout0), 1e-12)
    expect_gt(abs(sp$max_assurance_evaluated - with0), 0.1)
  }
})

test_that("m9: Monte Carlo marginal counts an empty arm as a failure (within 3 MCSE of exact, far from the credited value)", {
  skip_on_cran()
  B <- 20000
  prior <- c(170, 30)
  for (af in c("se", "sp")) {
    pp <- if (af == "se") c(2, 18) else c(18, 2)
    args <- list(prior_prev = pp, assurance_for = af, N_range = 5L,
                 n_range = 20L, B = B, seed = 2026)
    args[[paste0("prior_", af)]] <- prior
    args[[paste0("delta_", af)]] <- 0.14
    ex <- suppressWarnings(do.call(bam_sample_size, c(args, list(method = "exact"))))
    mc <- suppressWarnings(do.call(bam_sample_size, c(args, list(method = "monte_carlo"))))
    A <- ex$max_assurance_evaluated
    mcse <- sqrt(A * (1 - A) / B)
    credited <- marginal_empty(5L, af, prior, 0.14, pp, credit_empty = TRUE)
    expect_lte(abs(mc$max_assurance_evaluated - A), 3 * mcse)
    expect_gt(abs(mc$max_assurance_evaluated - credited), 0.1)
  }
})

test_that("VAP example of Wilson et al. 2022 (Se only, full width 0.16) gives N = 104, A = 0.8005", {
  res <- quick_w(
    prior_se = c(25.9, 2.1), prior_prev = c(29, 98), delta_se = 0.16,
    assurance_for = "se"
  )
  message("VAP: N = ", res$N_total, ", A = ", sprintf("%.10f", res$joint_assurance))
  expect_identical(res$N_total, 104L)
  expect_identical(round(res$joint_assurance, 4), 0.8005)
})

# Helper that fits a Beta prior to a given mean and interval width
fitb <- function(mu, wid = 0.1) {
  f <- function(s) {
    a <- mu * s; b <- (1 - mu) * s
    diff(qbeta(c(.025, .975), a, b)) - wid
  }
  s <- uniroot(f, c(2, 1e4))$root
  c(mu * s, (1 - mu) * s)
}

test_that("Fig 3 of Williamson et al. (rho_T = 0.30) gives N = 182, 154, 101", {
  skip_on_cran()
  pp <- fitb(0.30)
  priors_se <- list(c(9, 3), c(10, 2), c(11, 1))
  expected <- c(182L, 154L, 101L)
  got <- integer(3)
  for (i in 1:3) {
    res <- quick_w(
      prior_se = priors_se[[i]], prior_prev = pp, width_type = "lower_half_width",
      delta_se = 0.10, assurance_for = "se"
    )
    got[i] <- res$N_total
  }
  message("Fig 3, rho_T = 0.30: N = ", paste(got, collapse = ", "))
  expect_identical(got, expected)
})

test_that("Fig 3 of Williamson et al. (rho_T = 0.05) reports N (informative only; no assertion on the tail extrapolation)", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("DTA_RUN_SLOW_REF")), "set DTA_RUN_SLOW_REF to run (N up to ~1600, O(N^2) qbeta)")
  pp <- fitb(0.05)
  got <- integer(3)
  for (i in 1:3) {
    pr <- list(c(9, 3), c(10, 2), c(11, 1))[[i]]
    res <- quick_w(prior_se = pr, prior_prev = pp, width_type = "lower_half_width",
                   delta_se = 0.10, assurance_for = "se")
    got[i] <- res$N_total
  }
  message("Fig 3, rho_T = 0.05: N = ", paste(got, collapse = ", "),
          " (reference with tail extrapolation: 1628, 1294, 782)")
  succeed()
})

## ---------------------------------------------------------------------
## Errores de entrada
## ---------------------------------------------------------------------

test_that("invalid width_type and assurance_for are rejected", {
  expect_error(bam_sample_size(width_type = "half_width"), "width_type")
  expect_error(bam_sample_size(width_type = c("full_width", "x")), "width_type")
  expect_error(bam_sample_size(assurance_for = "both"), "assurance_for")
  expect_error(bam_sample_size(assurance_for = 1), "assurance_for")
})

## ---------------------------------------------------------------------
## Impresion
## ---------------------------------------------------------------------

test_that("print states width_type, assurance_for, level and the analysis priors in words", {
  res <- quick_w(prior_prev = c(4, 16), N_range = 600:700,
                 analysis_prior_se = c(1, 1))
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "full width of the equal-tailed 95% interval")
  expect_match(out, "joint")
  expect_match(out, "analysis prior", ignore.case = TRUE)
  expect_match(out, "a \\+ b = 2")

  res2 <- quick_w(prior_se = c(25, 7), delta_se = 0.10, N_range = 100:400,
                  width_type = "lower_half_width", assurance_for = "se")
  out2 <- paste(capture.output(print(res2)), collapse = "\n")
  expect_match(out2, "lower half-width: median minus the one-sided 95% lower bound")
  expect_match(out2, "marginal Se", ignore.case = TRUE)
})

test_that("print says analysis priors equal the design priors when they do", {
  res <- quick_w(prior_prev = c(4, 16), N_range = 600:700)
  out <- paste(capture.output(print(res)), collapse = "\n")
  expect_match(out, "same as design", ignore.case = TRUE)
})
