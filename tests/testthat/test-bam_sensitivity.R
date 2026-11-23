# Tests for bam_sensitivity(), bam_prior_contour() and bam_sensitivity_local()
# (v0.7.0).

# Fast scenarios (lenient widths): N is a few tens, so each exact call is
# well under a second.
fast_common <- list(prior_prev = c(4, 16), delta_se = 0.30, delta_sp = 0.30)

# ---------------------------------------------------------------------------
# (a) bam_sensitivity()
# ---------------------------------------------------------------------------

test_that("each row equals the direct bam_sample_size() call (identical)", {
  scen <- list(
    A = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
             delta_se = 0.30, delta_sp = 0.25),
    B = list(prior_se = c(9, 1), prior_sp = c(8, 2), prior_prev = c(6, 14),
             delta_se = 0.30, delta_sp = 0.30, assurance_for = "se",
             alpha_ci = 0.90),
    C = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
             delta_se = 0.15, delta_sp = 0.15, width_type = "lower_half_width",
             target_assurance = 0.70)
  )
  out <- suppressWarnings(bam_sensitivity(scen))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 3L)
  expect_identical(out$scenario, c("A", "B", "C"))
  for (nm in names(scen)) {
    d <- suppressWarnings(do.call(
      bam_sample_size, c(scen[[nm]], list(method = "exact"))
    ))
    dm1 <- suppressWarnings(do.call(
      bam_sample_size,
      c(scen[[nm]], list(method = "exact", N_range = d$N_total - 1L))
    ))
    row <- out[out$scenario == nm, ]
    expect_identical(row$N_total, d$N_total)
    expect_identical(row$A_N, d$joint_assurance)
    expect_identical(row$A_N_minus_1, dm1$max_assurance_evaluated)
    expect_true(row$target_reached)
    expect_true(row$A_N >= row$target_assurance)
    expect_true(row$A_N_minus_1 < row$target_assurance)
  }
  # conventions are reported explicitly
  expect_identical(out$assurance_for, c("joint", "se", "joint"))
  expect_identical(out$width_type,
                   c("full_width", "full_width", "lower_half_width"))
  expect_equal(out$alpha_ci, c(0.95, 0.90, 0.95))
})

test_that("P1 reproduces N = 678, A(N) and A(N - 1)", {
  skip_on_cran()
  out <- suppressWarnings(bam_sensitivity(list(
    P1 = list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
              delta_se = 0.14, delta_sp = 0.10)
  )))
  expect_identical(out$N_total, 678L)
  expect_equal(out$A_N, 0.8003489948, tolerance = 1e-10)
  expect_equal(out$A_N_minus_1, 0.7996848824, tolerance = 1e-10)
  expect_true(out$target_reached)
})

test_that("common arguments in ... fill what the scenario does not fix", {
  scen <- list(s1 = list(prior_se = c(17, 3)),
               s2 = list(prior_se = c(9, 1), delta_se = 0.35))
  out <- suppressWarnings(bam_sensitivity(
    scen, prior_sp = c(18, 2), prior_prev = c(4, 16),
    delta_se = 0.30, delta_sp = 0.30
  ))
  d1 <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
    delta_se = 0.30, delta_sp = 0.30, method = "exact"))
  d2 <- suppressWarnings(bam_sample_size(
    prior_se = c(9, 1), prior_sp = c(18, 2), prior_prev = c(4, 16),
    delta_se = 0.35, delta_sp = 0.30, method = "exact"))
  expect_identical(out$N_total, c(d1$N_total, d2$N_total))
  expect_identical(out$A_N, c(d1$joint_assurance, d2$joint_assurance))
  # the scenario value wins over the common one
  expect_equal(out$delta_se, c(0.30, 0.35))
})

test_that("a list and an equivalent data frame give the same result", {
  lst <- list(
    one = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
               delta_se = 0.30, delta_sp = 0.25),
    two = list(prior_se = c(9, 1), prior_sp = c(8, 2), prior_prev = c(6, 14),
               delta_se = 0.30, delta_sp = 0.30, assurance_for = "se",
               analysis_prior_se = c(1, 1)),
    three = list(prior_se = c(17, 3), prior_sp = c(18, 2),
                 prior_prev = c(4, 16), delta_se = 0.15, delta_sp = 0.15,
                 width_type = "lower_half_width")
  )
  df <- data.frame(
    scenario = c("one", "two", "three"),
    prior_se_a = c(17, 9, 17), prior_se_b = c(3, 1, 3),
    prior_sp_a = c(18, 8, 18), prior_sp_b = c(2, 2, 2),
    prior_prev_a = c(4, 6, 4), prior_prev_b = c(16, 14, 16),
    analysis_prior_se_a = c(NA, 1, NA), analysis_prior_se_b = c(NA, 1, NA),
    delta_se = c(0.30, 0.30, 0.15), delta_sp = c(0.25, 0.30, 0.15),
    width_type = c(NA, NA, "lower_half_width"),
    assurance_for = c(NA, "se", NA),
    stringsAsFactors = FALSE
  )
  o1 <- suppressWarnings(bam_sensitivity(lst))
  o2 <- suppressWarnings(bam_sensitivity(df))
  expect_identical(o1, o2)
})

test_that("an infeasible scenario does not abort the batch", {
  scen <- list(
    ok = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
              delta_se = 0.30, delta_sp = 0.30),
    impossible = list(prior_se = c(17, 3), prior_sp = c(18, 2),
                      prior_prev = c(4, 16), delta_se = 0.001,
                      delta_sp = 0.001, N_max = 40L)
  )
  out <- suppressWarnings(bam_sensitivity(scen))
  expect_equal(nrow(out), 2L)
  expect_true(out$target_reached[1])
  expect_false(out$target_reached[2])
  expect_true(is.na(out$N_total[2]))
  expect_true(is.na(out$A_N[2]))
  expect_true(is.na(out$A_N_minus_1[2]))
  expect_true(is.finite(out$max_assurance_evaluated[2]))
  expect_true(is.finite(out$N_at_max_assurance[2]))
  expect_true(out$max_assurance_evaluated[2] < out$target_assurance[2])
})

test_that("a dta_bam_priors object can be used as a scenario", {
  pri <- bam_prior_from_pilot(tp = 16, fn = 1, tn = 20, fp = 35,
                              base_prior_se = c(9.9, 1.1),
                              base_prior_sp = c(2, 2),
                              base_prior_prev = c(12, 43))
  out <- suppressWarnings(bam_sensitivity(
    list(pilot = list(priors = pri, delta_se = 0.30, delta_sp = 0.30))
  ))
  d <- suppressWarnings(bam_sample_size(
    priors = pri, delta_se = 0.30, delta_sp = 0.30, method = "exact"))
  expect_identical(out$N_total, d$N_total)
  expect_identical(out$A_N, d$joint_assurance)
})

test_that("bam_sensitivity() validates its input", {
  ok <- list(prior_se = c(17, 3))
  expect_error(bam_sensitivity(list(a = list(prior_xx = c(1, 1)))),
               "prior_xx")
  expect_error(bam_sensitivity(list(a = list(method = "monte_carlo"))),
               "method")
  expect_error(bam_sensitivity(list(a = list(seed = 1))), "seed")
  expect_error(bam_sensitivity("a"), "list|data frame")
  expect_error(bam_sensitivity(list()), "at least one")
  expect_error(bam_sensitivity(list(ok, 3)), "list")
  # a mix of named and unnamed entries is not allowed
  expect_error(bam_sensitivity(list(a = ok, ok)), "names")
  # duplicated names
  expect_error(bam_sensitivity(list(a = ok, a = ok)), "unique|duplicated")
  # unknown argument in ...
  expect_error(bam_sensitivity(list(a = ok), foo = 1), "foo")
  expect_error(bam_sensitivity(list(a = ok), method = "monte_carlo"),
               "method")
  # a priors object together with explicit priors in the same scenario
  pri <- bam_prior_from_pilot(1, 1, 1, 1, c(1, 1), c(1, 1), c(1, 1))
  expect_error(bam_sensitivity(list(a = list(priors = pri,
                                             prior_se = c(2, 2)))),
               "priors")
  # bad parallel / cores
  expect_error(bam_sensitivity(list(a = ok), parallel = NA), "parallel")
  expect_error(bam_sensitivity(list(a = ok), parallel = TRUE, cores = 0),
               "cores")
  # data frame: a half-specified pair
  df <- data.frame(prior_se_a = 17, delta_se = 0.3)
  expect_error(bam_sensitivity(df), "prior_se_b|both")
  # data frame: unknown column
  expect_error(bam_sensitivity(data.frame(foo = 1)), "foo")
  # an error raised inside bam_sample_size() names the scenario
  expect_error(
    suppressWarnings(bam_sensitivity(list(bad = list(prior_se = c(-1, 3))))),
    "bad"
  )
})

test_that("parallel = TRUE gives results identical to the sequential run", {
  skip_on_cran()
  scen <- list(
    A = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
             delta_se = 0.30, delta_sp = 0.25),
    B = list(prior_se = c(9, 1), prior_sp = c(8, 2), prior_prev = c(6, 14),
             delta_se = 0.30, delta_sp = 0.30, assurance_for = "se"),
    C = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
             delta_se = 0.001, delta_sp = 0.001, N_max = 30L)
  )
  seq_out <- suppressWarnings(bam_sensitivity(scen))
  par_out <- tryCatch(
    suppressWarnings(bam_sensitivity(scen, parallel = TRUE, cores = 2L)),
    error = function(e) e
  )
  if (inherits(par_out, "error")) {
    skip(paste("PSOCK cluster could not be started here:",
               conditionMessage(par_out)))
  }
  expect_identical(par_out, seq_out)
})

# ---------------------------------------------------------------------------
# (b) bam_prior_contour()
# ---------------------------------------------------------------------------

test_that("contour points are at distance epsilon (H2 and H), incl. a or b < 1", {
  priors <- list(c(9.9, 1.1), c(12, 43), c(17, 3), c(0.5, 0.5), c(0.4, 5),
                 c(7, 0.8), c(1, 1), c(2, 2))
  for (pr in priors) {
    for (dist in c("H2", "H")) {
      for (eps in c(0.00354, 0.02, 0.1)) {
        ct <- suppressWarnings(
          bam_prior_contour(pr, epsilon = eps, distance = dist))
        expect_equal(nrow(ct), 72L)
        expect_true(all(ct$a > 0) && all(ct$b > 0))
        d2 <- dtasamplesize:::.bam_hellinger_d2(ct$a, ct$b, pr[1], pr[2])
        dd <- if (dist == "H2") d2 else sqrt(d2)
        expect_lt(max(abs(dd - eps)), 1e-10)
        # reported distance agrees
        expect_lt(max(abs(ct$distance - eps)), 1e-10)
        # polar parametrisation
        expect_equal(ct$a, pr[1] + ct$r * cos(ct$phi), tolerance = 1e-12)
        expect_equal(ct$b, pr[2] + ct$r * sin(ct$phi), tolerance = 1e-12)
      }
    }
  }
})

test_that("the default angular grid is documented and reproducible", {
  c1 <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354)
  c2 <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354)
  expect_identical(c1, c2)
  expect_equal(c1$phi, -pi + 2 * pi * (0:71) / 72)
  # number of directions and explicit vector of angles
  c3 <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354, n_phi = 12L)
  expect_equal(nrow(c3), 12L)
  c4 <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354,
                          phi = c(-1, 0, 1))
  expect_equal(c4$phi, c(-1, 0, 1))
  expect_equal(nrow(c4), 3L)
})

test_that("H = sqrt(H2): the H contour at epsilon equals the H2 contour at epsilon^2", {
  ch <- bam_prior_contour(c(17, 3), epsilon = 0.06, distance = "H")
  c2 <- bam_prior_contour(c(17, 3), epsilon = 0.06^2, distance = "H2")
  expect_equal(ch$a, c2$a, tolerance = 1e-8)
  expect_equal(ch$b, c2$b, tolerance = 1e-8)
})

test_that("Wilson et al. 2022 contour extremes (H2 = 0.00354)", {
  skip_on_cran()
  se <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354, n_phi = 2000L)
  pv <- bam_prior_contour(c(12, 43), epsilon = 0.00354, n_phi = 2000L)
  expect_lt(max(abs(range(se$a) - c(7.63, 12.80))), 0.01)
  expect_lt(max(abs(range(se$b) - c(0.892, 1.362))), 0.01)
  expect_lt(max(abs(range(pv$a) - c(9.49, 15.19))), 0.01)
  expect_lt(max(abs(range(pv$b) - c(33.86, 54.60))), 0.01)
  message(sprintf(
    paste0("contour extremes (2000 directions): Se a [%.4f, %.4f] b ",
           "[%.4f, %.4f]; prev a [%.4f, %.4f] b [%.4f, %.4f]"),
    min(se$a), max(se$a), min(se$b), max(se$b),
    min(pv$a), max(pv$a), min(pv$b), max(pv$b)))
})

test_that("bam_prior_contour() validates its input", {
  expect_error(bam_prior_contour(c(1, 1)), "epsilon")
  expect_error(bam_prior_contour(c(1, 1), epsilon = 0), "epsilon")
  expect_error(bam_prior_contour(c(1, 1), epsilon = 1), "epsilon")
  expect_error(bam_prior_contour(c(1, 1), epsilon = -0.1), "epsilon")
  expect_error(bam_prior_contour(c(1, -1), epsilon = 0.01), "prior")
  expect_error(bam_prior_contour(1, epsilon = 0.01), "prior")
  expect_error(bam_prior_contour(c(1, 1), epsilon = 0.01, n_phi = 2),
               "n_phi")
  expect_error(bam_prior_contour(c(1, 1), epsilon = 0.01, phi = NA_real_),
               "phi")
  expect_error(bam_prior_contour(c(1, 1), epsilon = 0.01, distance = "TV"))
})

test_that("an unreachable epsilon gives a warning and NA points", {
  # H2 < 1 always, so epsilon very close to 1 may be unreachable within
  # numerical limits in the bounded directions
  expect_warning(
    ct <- bam_prior_contour(c(1, 1), epsilon = 1 - 1e-15, n_phi = 8L),
    "epsilon|reach|root"
  )
  expect_true(is.data.frame(ct))
})

# ---------------------------------------------------------------------------
# (b) bam_sensitivity_local()
# ---------------------------------------------------------------------------

local_args <- function(...) {
  utils::modifyList(list(
    base_prior_se = c(9.9, 1.1), base_prior_sp = c(18, 2),
    base_prior_prev = c(12, 43), delta_se = 0.30, delta_sp = 0.30,
    pilot = list(tp = 16, fn = 1, tn = 20, fp = 35)
  ), list(...))
}

  # epsilon 1e-7 (not 1e-8): with the absolute floor of the precision bound,
  # Beta(9.9, 1.1) at 1e-8 has a bound of 1.2e-6 > 1e-6 and
  # every direction is flagged; at 1e-7 the bound is 1.2e-7.
test_that("epsilon ~ 0 reproduces the base prior (A and n*)", {
  res <- suppressWarnings(do.call(bam_sensitivity_local, c(
    local_args(), list(which = "se", epsilon = 1e-7, n_phi = 6L)
  )))
  expect_s3_class(res, "dta_bam_sens_local")
  base <- res$base
  expect_equal(res$table$A_nref, rep(base$A_nref, 6L), tolerance = 1e-6)
  expect_true(all(res$table$n_star == base$n_star))
  expect_equal(res$summary$A_min, base$A_nref, tolerance = 1e-6)
  expect_equal(res$summary$A_max, base$A_nref, tolerance = 1e-6)
  expect_identical(res$summary$n_star_min, base$n_star)
  expect_identical(res$summary$n_star_max, base$n_star)
})

test_that("end-to-end example (few directions) is fast and consistent", {
  t0 <- proc.time()[["elapsed"]]
  res <- suppressWarnings(do.call(bam_sensitivity_local, c(
    local_args(assurance_for = "se"),
    list(which = "se", epsilon = 0.00354, n_phi = 8L)
  )))
  el <- proc.time()[["elapsed"]] - t0
  # Loose wall-clock bound: it only catches a gross slowdown. The run takes
  # about 5 to 11 s on the reference laptop depending on load, so a tight
  # bound would test the machine, not the code.
  expect_lt(el, 60)
  tb <- res$table
  expect_equal(nrow(tb), 8L)
  expect_true(all(c("phi", "a_init", "b_init", "distance", "A_nref",
                    "n_star", "target_reached") %in% names(tb)))
  # contour distance is epsilon at every point
  expect_lt(max(abs(tb$distance - 0.00354)), 1e-10)
  # summary agrees with the table
  expect_equal(res$summary$A_min, min(tb$A_nref))
  expect_equal(res$summary$A_max, max(tb$A_nref))
  expect_equal(res$summary$n_star_min, min(tb$n_star, na.rm = TRUE))
  expect_equal(res$summary$n_star_max, max(tb$n_star, na.rm = TRUE))
  i <- which.min(tb$A_nref)
  expect_equal(res$summary$A_min_a_init, tb$a_init[i])
  expect_equal(res$summary$A_min_b_init, tb$b_init[i])
  # perturbed prior is re-updated with the pilot: check one point by hand
  k <- 3L
  pri <- bam_prior_from_pilot(16, 1, 20, 35, c(tb$a_init[k], tb$b_init[k]),
                              c(18, 2), c(12, 43))
  d <- suppressWarnings(bam_sample_size(
    priors = pri, delta_se = 0.30, delta_sp = 0.30, assurance_for = "se",
    method = "exact"))
  expect_identical(tb$n_star[k], d$N_total)
  dn <- suppressWarnings(bam_sample_size(
    priors = pri, delta_se = 0.30, delta_sp = 0.30, assurance_for = "se",
    method = "exact", N_range = res$n_ref))
  expect_identical(tb$A_nref[k], dn$max_assurance_evaluated)
  # default n_ref is the base n*
  expect_identical(res$n_ref, res$base$n_star)
  # printing works
  expect_output(print(res), "Hellinger")
})

test_that("n_ref, no pilot, sp and prev parameters work", {
  r_sp <- suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(17, 3), base_prior_sp = c(18, 2),
    base_prior_prev = c(4, 16), which = "sp", epsilon = 0.01, n_phi = 6L,
    n_ref = 60, delta_se = 0.30, delta_sp = 0.30))
  expect_identical(r_sp$n_ref, 60L)
  expect_equal(nrow(r_sp$table), 6L)
  # no pilot: the perturbed prior is used as is, for design and analysis
  k <- 2L
  d <- suppressWarnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(r_sp$table$a_init[k], r_sp$table$b_init[k]),
    prior_prev = c(4, 16), delta_se = 0.30, delta_sp = 0.30,
    method = "exact"))
  expect_identical(r_sp$table$n_star[k], d$N_total)

  r_pv <- suppressWarnings(bam_sensitivity_local(
    base_prior_se = c(17, 3), base_prior_sp = c(18, 2),
    base_prior_prev = c(4, 16), which = "prev", epsilon = 0.01, n_phi = 6L,
    n_ref = 60, delta_se = 0.30, delta_sp = 0.30))
  expect_equal(nrow(r_pv$table), 6L)
  expect_true(all(is.finite(r_pv$table$A_nref)))
})

test_that("bam_sensitivity_local() parallel equals sequential", {
  skip_on_cran()
  args <- c(local_args(assurance_for = "se"),
            list(which = "se", epsilon = 0.00354, n_phi = 6L))
  s <- suppressWarnings(do.call(bam_sensitivity_local, args))
  p <- tryCatch(
    suppressWarnings(do.call(bam_sensitivity_local,
                             c(args, list(parallel = TRUE, cores = 2L)))),
    error = function(e) e)
  if (inherits(p, "error")) {
    skip(paste("PSOCK cluster could not be started here:",
               conditionMessage(p)))
  }
  expect_identical(p$table, s$table)
  expect_identical(p$summary, s$summary)
})

test_that("bam_sensitivity_local() validates its input", {
  a <- local_args()
  expect_error(do.call(bam_sensitivity_local, c(a, list(which = "se"))),
               "epsilon")
  expect_error(do.call(bam_sensitivity_local,
                       c(a, list(which = "xx", epsilon = 0.01))))
  expect_error(do.call(bam_sensitivity_local,
                       c(a, list(epsilon = 0.01, prior_se = c(1, 1)))),
               "prior_se")
  expect_error(do.call(bam_sensitivity_local,
                       c(a, list(epsilon = 0.01, method = "monte_carlo"))),
               "method")
  expect_error(do.call(bam_sensitivity_local,
                       c(a, list(epsilon = 0.01, n_ref = 1))), "n_ref")
  expect_error(do.call(bam_sensitivity_local,
                       c(a, list(epsilon = 0.01, n_ref = c(10, 20)))),
               "n_ref")
  bad <- a; bad$pilot <- list(tp = 1)
  expect_error(do.call(bam_sensitivity_local, c(bad, list(epsilon = 0.01))),
               "pilot")
  expect_error(bam_sensitivity_local(base_prior_se = c(1, 1)), "base_prior")
})
