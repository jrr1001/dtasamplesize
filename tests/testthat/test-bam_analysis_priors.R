# Separate design and analysis priors in bam_sample_size() (v0.7.0) and
# acceptance of a "dta_bam_priors" object.
#
# Convention: the DESIGN prior generates the true values (and, in exact mode,
# the Beta-Binomial predictive of the count); the ANALYSIS prior defines the
# posterior Beta(a_A + x, b_A + n - x) whose width is compared with delta.
# Prevalence acts in the design only.

p1_args <- list(
  prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
  delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
  method = "exact"
)
p2_args <- list(
  prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
  delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
  method = "exact"
)

# Numeric fields that must match bit for bit
num_fields <- c(
  "N_total", "n_total", "joint_assurance", "max_assurance_evaluated",
  "N_at_max_assurance", "assurance_mcse", "n_diseased", "n_non_diseased",
  "assurance_se", "assurance_sp", "N_total_median", "N_total_P75",
  "N_total_P90", "buderer_n_se", "N_range_used", "target_reached"
)

quick <- function(args, ...) {
  # Small n_range and B: they only affect the per-arm legacy fields
  suppressWarnings(do.call(bam_sample_size, c(args, list(
    n_range = 20:150, B = 1500, seed = 11
  ), list(...))))
}

## ---------------------------------------------------------------------
## Identity: analysis_prior_* equal to the design prior (explicit) == NULL
## ---------------------------------------------------------------------

test_that("analysis priors equal to the design priors reproduce NULL bit by bit (P1)", {
  skip_on_cran()
  ref <- quick(p1_args)
  alt <- quick(p1_args, analysis_prior_se = c(17, 3), analysis_prior_sp = c(2, 2))
  expect_identical(ref$N_total, 678L)
  for (f in num_fields) expect_identical(alt[[f]], ref[[f]], info = f)
})

test_that("analysis priors equal to the design priors reproduce NULL bit by bit (P2)", {
  skip_on_cran()
  ref <- quick(p2_args)
  alt <- quick(p2_args, analysis_prior_se = c(17, 3), analysis_prior_sp = c(18, 2))
  expect_identical(ref$N_total, 672L)
  for (f in num_fields) expect_identical(alt[[f]], ref[[f]], info = f)
})

test_that("identity also holds in Monte Carlo mode (legacy fields and headline)", {
  args <- c(p1_args[setdiff(names(p1_args), "method")], list(method = "monte_carlo"))
  ref <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    n_range = 20:120, N_range = c(400, 600), B = 800, seed = 5
  ))))
  alt <- suppressWarnings(do.call(bam_sample_size, c(args, list(
    n_range = 20:120, N_range = c(400, 600), B = 800, seed = 5,
    analysis_prior_se = c(17, 3), analysis_prior_sp = c(2, 2)
  ))))
  for (f in num_fields) expect_identical(alt[[f]], ref[[f]], info = f)
})

test_that("the effective analysis priors and the 'given' indicator are stored", {
  r0 <- quick(p1_args, N_range = 600:700)
  expect_identical(r0$analysis_prior_se, c(17, 3))
  expect_identical(r0$analysis_prior_sp, c(2, 2))
  expect_identical(unname(r0$analysis_prior_given), c(FALSE, FALSE))
  expect_identical(names(r0$analysis_prior_given), c("se", "sp"))
  expect_identical(unname(r0$analysis_prior_size), c(20, 4))

  r1 <- quick(p1_args, N_range = 600:700, analysis_prior_se = c(1, 1))
  expect_identical(r1$analysis_prior_se, c(1, 1))
  expect_identical(r1$analysis_prior_sp, c(2, 2))
  expect_identical(unname(r1$analysis_prior_given), c(TRUE, FALSE))
  expect_identical(unname(r1$analysis_prior_size), c(2, 4))
  # the design prior is kept
  expect_identical(r1$prior_se, c(17, 3))
})

## ---------------------------------------------------------------------
## Bridge: near-point design prior + flat analysis prior
## ---------------------------------------------------------------------

bridge_args <- function(analysis_se, analysis_sp) {
  s <- 1e6
  list(
    prior_se = c(0.85 * s, 0.15 * s), prior_sp = c(0.90 * s, 0.10 * s),
    prior_prev = c(0.20 * s, 0.80 * s),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    method = "exact",
    analysis_prior_se = analysis_se, analysis_prior_sp = analysis_sp
  )
}

test_that("point design with Jeffreys analysis prior gives N close to 576", {
  skip_on_cran()
  res <- quick(bridge_args(c(0.5, 0.5), c(0.5, 0.5)))
  message("BRIDGE Jeffreys: N = ", res$N_total, ", A = ",
          sprintf("%.7f", res$joint_assurance))
  expect_lte(abs(res$N_total - 576L), 1L)
})

test_that("point design with uniform analysis prior gives N close to 579", {
  skip_on_cran()
  res <- quick(bridge_args(c(1, 1), c(1, 1)))
  message("BRIDGE uniform: N = ", res$N_total, ", A = ",
          sprintf("%.7f", res$joint_assurance))
  expect_lte(abs(res$N_total - 579L), 1L)
})

test_that("point design with informative analysis priors (17,3)/(18,2) gives N close to 466", {
  skip_on_cran()
  res <- quick(bridge_args(c(17, 3), c(18, 2)))
  message("BRIDGE informative: N = ", res$N_total, ", A = ",
          sprintf("%.7f", res$joint_assurance))
  expect_lte(abs(res$N_total - 466L), 1L)
})

## ---------------------------------------------------------------------
## Warning when the analysis prior already meets the target
## ---------------------------------------------------------------------

collect_warnings <- function(expr) {
  w <- character()
  withCallingHandlers(
    expr,
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    }
  )
  w
}

test_that("a narrow analysis prior whose prior-only width is <= delta warns", {
  w <- collect_warnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, method = "exact",
    analysis_prior_se = c(500, 100), n_range = 20:30, B = 1000, N_range = 50:60
  ))
  expect_true(any(grepl("analysis prior", w) & grepl("Se", w)))
  expect_warning(
    withCallingHandlers(
      bam_sample_size(
        prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
        delta_se = 0.14, delta_sp = 0.10, method = "exact",
        analysis_prior_se = c(500, 100), n_range = 20:30, B = 1000,
        N_range = 50:60
      ),
      # silences only the inherited warnings (other causes)
      warning = function(cond) {
        if (!grepl("analysis prior", conditionMessage(cond)))
          invokeRestart("muffleWarning")
      }
    ),
    "analysis prior"
  )
  # Without an explicit analysis prior, and with deltas that the design prior
  # does not meet without data, this warning is not issued
  w0 <- collect_warnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, method = "exact",
    n_range = 20:30, B = 1000, N_range = 50:80
  ))
  expect_false(any(grepl("analysis prior", w0)))
  # The warning follows the EFFECTIVE analysis prior: if none was given
  # explicitly (equal to the design prior) and it already meets the target
  # without data, it also warns, and the message says so
  w1 <- collect_warnings(bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.30, delta_sp = 0.30, method = "exact",
    n_range = 20:30, B = 1000, N_range = 50:80
  ))
  expect_true(any(grepl("analysis prior for Se", w1) &
                    grepl("equal to the design prior", w1)))
})

test_that("the warning names the arm; an unused arm does not warn", {
  args <- list(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, method = "exact",
    analysis_prior_sp = c(2000, 2000), n_range = 20:30, B = 1000, N_range = 50:60
  )
  w <- collect_warnings(do.call(bam_sample_size, args))
  expect_true(any(grepl("analysis prior", w) & grepl("Sp", w)))
  # with assurance_for = "se" the Sp arm is not used: no warning
  w2 <- collect_warnings(do.call(bam_sample_size, c(args, list(assurance_for = "se"))))
  expect_false(any(grepl("analysis prior", w2)))
})

## ---------------------------------------------------------------------
## Input errors for the analysis priors
## ---------------------------------------------------------------------

test_that("malformed analysis priors are rejected", {
  f <- function(...) bam_sample_size(
    prior_prev = c(4, 16), n_range = 20:30, B = 1000, N_range = 50:60, ...
  )
  expect_error(f(analysis_prior_se = c(1, 2, 3)), "analysis_prior_se")
  expect_error(f(analysis_prior_se = c(-1, 2)), "analysis_prior_se")
  expect_error(f(analysis_prior_se = c(0, 2)), "analysis_prior_se")
  expect_error(f(analysis_prior_sp = "a"), "analysis_prior_sp")
  expect_error(f(analysis_prior_sp = c(NA, 2)), "analysis_prior_sp")
  expect_error(f(analysis_prior_sp = 5), "analysis_prior_sp")
})

## ---------------------------------------------------------------------
## Monte Carlo with separate priors vs. exact, at fixed N (<= 3 MCSE)
## ---------------------------------------------------------------------

test_that("Monte Carlo with separate priors agrees with exact at fixed N (<= 3 MCSE)", {
  skip_on_cran()
  base <- list(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    analysis_prior_se = c(1, 1), analysis_prior_sp = c(1, 1)
  )
  for (N in c(300L, 678L)) {
    ex <- suppressWarnings(do.call(bam_sample_size, c(base, list(
      method = "exact", N_range = N, n_range = 20:20
    ))))
    mc <- suppressWarnings(do.call(bam_sample_size, c(base, list(
      method = "monte_carlo", N_range = N, B = 20000, seed = 2026, n_range = 20:20
    ))))
    A <- ex$max_assurance_evaluated
    mcse <- sqrt(A * (1 - A) / 20000)
    message("MC vs exact (separate priors), N = ", N, ": exact ",
            sprintf("%.6f", A), ", MC ", sprintf("%.6f", mc$max_assurance_evaluated),
            ", diff/MCSE = ",
            sprintf("%.2f", (mc$max_assurance_evaluated - A) / mcse))
    expect_lte(abs(mc$max_assurance_evaluated - A), 3 * mcse)
  }
})

## ---------------------------------------------------------------------
## dta_bam_priors object (fixed contract; built by hand)
## ---------------------------------------------------------------------

make_priors_obj <- function(use_in) {
  base <- list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16))
  design <- list(prior_se = c(21, 4), prior_sp = c(6, 3), prior_prev = c(7, 22))
  analysis <- if (use_in == "both") {
    list(analysis_prior_se = design$prior_se, analysis_prior_sp = design$prior_sp)
  } else {
    list(analysis_prior_se = base$prior_se, analysis_prior_sp = base$prior_sp)
  }
  structure(c(
    design, analysis,
    list(base = base,
         pilot = list(tp = 4, fn = 1, tn = 4, fp = 1),
         a0 = 0.5, design = "cohort", use_in = use_in)
  ), class = "dta_bam_priors")
}

cmp_args <- list(
  delta_se = 0.25, delta_sp = 0.20, target_assurance = 0.80,
  method = "exact", n_range = 20:60, B = 1000, seed = 3
)

test_that("priors object, use_in = 'both' (analysis = design) equals passing the priors by hand", {
  po <- make_priors_obj("both")
  via_obj <- suppressWarnings(do.call(bam_sample_size, c(list(priors = po), cmp_args)))
  by_hand <- suppressWarnings(do.call(bam_sample_size, c(list(
    prior_se = c(21, 4), prior_sp = c(6, 3), prior_prev = c(7, 22)
  ), cmp_args)))
  for (f in num_fields) expect_identical(via_obj[[f]], by_hand[[f]], info = f)
  expect_identical(via_obj$prior_se, c(21, 4))
  expect_identical(via_obj$prior_prev, c(7, 22))
  expect_identical(via_obj$analysis_prior_se, c(21, 4))
})

test_that("priors object, use_in = 'design' (analysis = base) equals passing analysis_prior_* by hand", {
  po <- make_priors_obj("design")
  via_obj <- suppressWarnings(do.call(bam_sample_size, c(list(priors = po), cmp_args)))
  by_hand <- suppressWarnings(do.call(bam_sample_size, c(list(
    prior_se = c(21, 4), prior_sp = c(6, 3), prior_prev = c(7, 22),
    analysis_prior_se = c(17, 3), analysis_prior_sp = c(2, 2)
  ), cmp_args)))
  for (f in num_fields) expect_identical(via_obj[[f]], by_hand[[f]], info = f)
  expect_identical(via_obj$analysis_prior_se, c(17, 3))
  expect_identical(via_obj$analysis_prior_sp, c(2, 2))
  # and it differs from the use_in = "both" case
  both <- suppressWarnings(do.call(bam_sample_size, c(list(priors = make_priors_obj("both")), cmp_args)))
  expect_false(identical(both$max_assurance_evaluated, via_obj$max_assurance_evaluated))
})

test_that("priors object conflicts with explicitly passed prior arguments", {
  po <- make_priors_obj("both")
  call_with <- function(...) do.call(bam_sample_size, c(list(priors = po), cmp_args, list(...)))
  expect_error(call_with(prior_se = c(17, 3)), "priors")
  expect_error(call_with(prior_sp = c(2, 2)), "priors")
  expect_error(call_with(prior_prev = c(4, 16)), "priors")
  expect_error(call_with(analysis_prior_se = c(1, 1)), "priors")
  expect_error(call_with(analysis_prior_sp = c(1, 1)), "priors")
})

test_that("priors that is neither NULL nor a dta_bam_priors object is an error", {
  expect_error(bam_sample_size(priors = list(prior_se = c(1, 1))), "dta_bam_priors")
  expect_error(bam_sample_size(priors = c(1, 1)), "dta_bam_priors")
})

## ---------------------------------------------------------------------
## Argument-name compatibility with v0.6.9
## ---------------------------------------------------------------------

test_that("target = 0.85 still partial-matches target_assurance (no new argument starts with 'target')", {
  expect_false(any(startsWith(setdiff(names(formals(bam_sample_size)),
                                      "target_assurance"), "target")))
  res <- suppressWarnings(bam_sample_size(
    target = 0.85, prior_prev = c(4, 16), delta_se = 0.30, delta_sp = 0.30,
    n_range = 20:60, B = 1000, N_range = 40:200
  ))
  expect_identical(res$target_assurance, 0.85)
})
