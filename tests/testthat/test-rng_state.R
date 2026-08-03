
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
