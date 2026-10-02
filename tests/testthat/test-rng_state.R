
# Every exported function that calls set.seed() internally.
seeded_calls <- list(
  mc_validate_buderer = function() {
    mc_validate_buderer(B = 100, seed = 7)
  },
  bam_sample_size = function() {
    bam_sample_size(B = 100, seed = 7, n_range = 20:600)
  },
  joint_sample_size = function() {
    joint_sample_size(B = 100, seed = 7, N_range = seq(100, 800, by = 50))
  },
  ss_net_benefit = function() {
    ss_net_benefit(pt_range = 0.30, B = 100, seed = 7)
  },
  ss_adaptive_prevalence = function() {
    ss_adaptive_prevalence(B = 20, seed = 7, prev_true_range = 0.30)
  },
  ss_imperfect_ref = function() {
    ss_imperfect_ref(B = 100, seed = 7, sensitivity_table = FALSE)
  },
  ss_unified = function() {
    ss_unified(B = 20, seed = 7, N_range = seq(200, 1000, by = 100),
               delta_auc = 0, check_nb = FALSE)
  }
)

test_that("no exported function disturbs the caller's RNG stream", {
  for (nm in names(seeded_calls)) {
    # Reference: what the user's next draw should be.
    set.seed(20260714)
    invisible(runif(3))
    expected <- runif(2)

    # Same stream, but with a package call interposed.
    set.seed(20260714)
    invisible(runif(3))
    invisible(suppressWarnings(seeded_calls[[nm]]()))
    actual <- runif(2)

    expect_identical(actual, expected,
                     info = paste0(nm, "() perturbed the caller's RNG stream"))
  }
})

test_that(".Random.seed itself is restored byte for byte", {
  for (nm in names(seeded_calls)) {
    set.seed(4242)
    invisible(rnorm(1))
    before <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)

    invisible(suppressWarnings(seeded_calls[[nm]]()))

    after <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    expect_identical(after, before,
                     info = paste0(nm, "() left .Random.seed modified"))
  }
})

test_that("a function called twice in a row still reproduces itself", {
  # Restoring the caller's seed must not break the function's OWN determinism:
  # it seeds itself internally, so two identical calls must agree exactly.
  set.seed(1)
  a <- mc_validate_buderer(B = 200, seed = 11)$results$P_width_target
  invisible(runif(5))
  b <- mc_validate_buderer(B = 200, seed = 11)$results$P_width_target
  expect_identical(a, b)
})

test_that("ss_time_dependent_roc also restores the caller's RNG stream", {
  # The eighth seeded function. Kept separate because it requires timeROC;
  # without this block NEWS's claim of "all eight exported functions" is untrue.
  skip_if_not_installed("timeROC")
  skip_if_not_installed("survival")

  set.seed(20260714)
  invisible(runif(3))
  expected <- runif(2)

  set.seed(20260714)
  invisible(runif(3))
  invisible(suppressWarnings(
    ss_time_dependent_roc(N_range = c(100, 150), B = 20,
                          censoring_rates = 0.2, seed = 7)
  ))
  actual <- runif(2)

  expect_identical(actual, expected)
})

test_that("the functions do not CREATE .Random.seed in a virgin session", {
  # If the caller has never drawn a random number, .Random.seed does not exist
  # in globalenv. A well-behaved package must not leave one behind.
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
    rm(".Random.seed", envir = .GlobalEnv)
  }
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  invisible(mc_validate_buderer(B = 50, seed = 7))

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("ss_unified() restores RNGkind() from a virgin session, not just .Random.seed's value", {
  # ss_unified() switches to the L'Ecuyer-CMRG generator internally (see its
  # @details) to draw an independent RNG stream per N via
  # parallel::nextRNGStream(). Restoring only the .Random.seed VALUE on exit
  # is not enough to undo that: when the calling session has never drawn a
  # random number before (.Random.seed does not exist), there is nothing
  # for a later RNGkind()/random-draw query to resync from, so the RNG kind
  # stays stuck on L'Ecuyer-CMRG even though .Random.seed itself is
  # correctly left absent again. (When .Random.seed already exists at
  # entry, as in every other test_that() above, a later query happens to
  # self-heal the kind from it, which is why this defect only shows up
  # from a virgin session -- exactly the condition under which
  # validation/make_manuscript_assets.R calls ss_unified() for Figure 4.)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
    rm(".Random.seed", envir = .GlobalEnv)
  }
  old_kind <- RNGkind()

  invisible(suppressWarnings(
    ss_unified(B = 20, seed = 7, N_range = seq(200, 1000, by = 100),
               delta_auc = 0, check_nb = FALSE)
  ))

  expect_identical(RNGkind(), old_kind)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("a later seeded call is unaffected by an earlier ss_unified() call from a virgin session", {
  # The concrete, measured consequence of the RNGkind() leak above: every
  # seeded function in this package reseeds itself with set.seed(seed) and
  # no explicit `kind` argument, which reuses whatever kind is CURRENTLY
  # active rather than whatever kind the caller's (restored) .Random.seed
  # happens to encode. So a downstream call with its own fixed seed gave a
  # different result purely because ss_unified() had run earlier in the
  # same session -- exactly the corruption that reached the manuscript's
  # "Joint Se/Sp + AUC" comparison-table row (0.8136 instead of 0.8041).
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
    rm(".Random.seed", envir = .GlobalEnv)
  }

  isolated <- joint_sample_size(Se = 0.85, Sp = 0.90, delta_se = 0.07,
                                 delta_sp = 0.05, prev = 0.20,
                                 design = "cohort", target_prob = 0.80,
                                 B = 2000, seed = 2026)

  invisible(suppressWarnings(
    ss_unified(prior_se = c(17, 3), prior_sp = c(18, 2),
               Se_ref = 0.95, Sp_ref = 0.98,
               delta_se = 0.07, delta_sp = 0.06, delta_auc = 0,
               N_range = seq(300, 900, by = 100), B = 100)
  ))

  after <- joint_sample_size(Se = 0.85, Sp = 0.90, delta_se = 0.07,
                              delta_sp = 0.05, prev = 0.20,
                              design = "cohort", target_prob = 0.80,
                              B = 2000, seed = 2026)

  expect_identical(after$n_total, isolated$n_total)
  expect_identical(after$joint_prob_se_sp, isolated$joint_prob_se_sp)
})

# --- defect 2: generator-inheritance (results depended on the caller's ---
# --- active RNGkind, not just on `seed`) ----------------------------------

test_that("the six generator-inheritance-fixed functions give identical results under Mersenne-Twister vs Wichmann-Hill", {
  # Pre-fix, each of these called set.seed(seed) with no `kind`, so it
  # reused whichever generator the CALLER had active -- for a package whose
  # whole point is reproducibility from `seed` alone, the result silently
  # depended on something the caller's `seed` argument said nothing about.
  # Each now fixes kind = "Mersenne-Twister" explicitly, so the result must
  # be identical regardless of what the caller had active beforehand.
  generator_invariant_calls <- list(
    mc_validate_buderer = function() mc_validate_buderer(B = 300, seed = 11),
    bam_sample_size = function() {
      bam_sample_size(B = 100, seed = 11, n_range = 20:300,
                       N_range = seq(100, 400, by = 50))
    },
    joint_sample_size = function() {
      joint_sample_size(B = 300, seed = 11, N_range = seq(100, 700, by = 50))
    },
    ss_net_benefit = function() {
      ss_net_benefit(pt_range = 0.30, B = 300, seed = 11,
                      N_range = seq(50, 600, by = 25))
    },
    ss_adaptive_prevalence = function() {
      ss_adaptive_prevalence(B = 50, seed = 11, prev_true_range = 0.30)
    }
  )

  old_kind <- RNGkind()
  on.exit(suppressWarnings(
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
  ), add = TRUE)

  for (nm in names(generator_invariant_calls)) {
    RNGkind("Mersenne-Twister")
    ref <- suppressWarnings(generator_invariant_calls[[nm]]())

    RNGkind("Wichmann-Hill")
    alt <- suppressWarnings(generator_invariant_calls[[nm]]())

    expect_equal(alt, ref,
                 info = paste0(nm, "() changed when the caller's RNGkind was ",
                                "Wichmann-Hill instead of Mersenne-Twister"))
  }
})

test_that("ss_time_dependent_roc also gives identical results under Mersenne-Twister vs Wichmann-Hill", {
  # Kept separate because it requires timeROC (as in the block above).
  skip_if_not_installed("timeROC")
  skip_if_not_installed("survival")

  old_kind <- RNGkind()
  on.exit(suppressWarnings(
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
  ), add = TRUE)

  call_it <- function() {
    suppressWarnings(ss_time_dependent_roc(
      N_range = c(100, 150), B = 20, censoring_rates = 0.2, seed = 11
    ))
  }

  RNGkind("Mersenne-Twister")
  ref <- call_it()
  RNGkind("Wichmann-Hill")
  alt <- call_it()

  expect_equal(alt, ref)
})

test_that("ss_unified()'s check_nb ceiling (nb_assurance_ceiling) is invariant to the caller's active RNGkind (defect 1)", {
  # nb_assurance_ceiling() (internal to R/ss_unified.R) saved/restored
  # .Random.seed but called set.seed(seed) with no `kind`, the same class
  # of defect as the six functions above, applied to the check_nb ceiling
  # helper rather than to ss_unified()'s own engine (which already fixed
  # its generator explicitly). nb_B_ceiling is kept small here for speed;
  # this test is about generator-independence, not ceiling precision (see
  # test-ss_unified.R for the precision fix).
  old_kind <- RNGkind()
  on.exit(suppressWarnings(
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
  ), add = TRUE)

  call_ceiling <- function() {
    suppressWarnings(ss_unified(
      B = 10, seed = 11, N_range = 300, delta_auc = 0,
      check_nb = TRUE, nb_B_ceiling = 20000
    ))$nb_ceiling
  }

  RNGkind("Mersenne-Twister")
  ref <- call_ceiling()
  RNGkind("Wichmann-Hill")
  alt <- call_ceiling()

  expect_identical(alt, ref)
})

test_that("published values reproduce identically under all five RNG generators (defect 2, mandatory check)", {
  # The two figures explicitly checked against adversarial verification:
  # ss_net_benefit(B = 4000, seed = 2026)'s N_conservative = 240 (the
  # manuscript's Figure 3/Table headline number), and
  # mc_validate_buderer()'s P_width_target at the Buderer n for
  # Se = 0.85, d = 0.07 = 0.573 (the 57% figure in Figure 1). Both must be
  # bit-identical across generators now that every set.seed() call in
  # these functions fixes kind = "Mersenne-Twister" explicitly.
  kinds <- c("Mersenne-Twister", "L'Ecuyer-CMRG", "Wichmann-Hill",
             "Marsaglia-Multicarry", "Knuth-TAOCP-2002")
  old_kind <- RNGkind()
  on.exit(suppressWarnings(
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
  ), add = TRUE)

  nb_by_kind <- vapply(kinds, function(k) {
    suppressWarnings(RNGkind(k))
    suppressWarnings(ss_net_benefit(B = 4000, seed = 2026)$N_conservative)
  }, integer(1))

  n_bud <- buderer_n(0.85, 0.07)
  mc_by_kind <- vapply(kinds, function(k) {
    suppressWarnings(RNGkind(k))
    mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = n_bud, B = 4000,
                         ci_method = "wald", seed = 2026)$results$P_width_target[1]
  }, numeric(1))

  # The test's real purpose is RNGkind-invariance *within this R version*
  # (every kind must agree with every other), not reproducing one fixed
  # literal captured under R <= 4.5's sampler -- R-devel/4.6.x changed the
  # RNG stream used by some sampling primitives, which shifts the absolute
  # figure without breaking cross-kind agreement. Compare every kind to the
  # first instead of to a hardcoded literal; keep the historical literals
  # as a documented, version-guarded sanity check.
  expect_true(all(nb_by_kind == nb_by_kind[1]))
  expect_true(all(mc_by_kind == mc_by_kind[1]))
  if (getRversion() < "4.6.0") {
    expect_true(all(nb_by_kind == 240L))
    expect_true(all(mc_by_kind == 0.573))
  }
})

test_that("calling ss_unified() twice in a row still reproduces itself", {
  # Mirrors the "a function called twice in a row" check above, specifically
  # for ss_unified() -- restoring RNGkind() on exit must not, itself, make
  # the function's OWN internal seeding path non-reproducible.
  a <- suppressWarnings(ss_unified(
    B = 20, seed = 7, N_range = seq(200, 1000, by = 100),
    delta_auc = 0, check_nb = FALSE
  ))
  b <- suppressWarnings(ss_unified(
    B = 20, seed = 7, N_range = seq(200, 1000, by = 100),
    delta_auc = 0, check_nb = FALSE
  ))
  expect_identical(a$grid_results, b$grid_results)
  expect_identical(a$N_effective, b$N_effective)
})
