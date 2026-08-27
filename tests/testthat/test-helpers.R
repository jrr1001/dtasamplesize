test_that("buderer_n reproduces known values", {
  # Se=0.85, d=0.07 -> n=100
  expect_equal(buderer_n(0.85, 0.07), 100)
  # Se=0.90, d=0.05 -> n=139
  expect_equal(buderer_n(0.90, 0.05), 139)
})

test_that("wilson_ci returns valid interval", {
  ci <- wilson_ci(85, 100)
  expect_true(ci["lower"] >= 0)
  expect_true(ci["upper"] <= 1)
  expect_true(ci["lower"] < ci["upper"])
  expect_true(ci["width"] > 0)
})

test_that("wald_ci returns valid interval", {
  ci <- wald_ci(85, 100)
  expect_true(ci["lower"] >= 0)
  expect_true(ci["upper"] <= 1)
  expect_true(ci["lower"] < ci["upper"])
  expect_true(ci["width"] > 0)
})

test_that("hanley_mcneil_var is positive", {
  v <- hanley_mcneil_var(0.75, 100, 100)
  expect_true(v > 0)
})

# --- wilson_lower() (internal): the B = 1 vacuous-bound fix -------------

wilson_lower <- get("wilson_lower", envir = asNamespace("dtasamplesize"))

test_that("wilson_lower stays in [0, 1] and does not collapse at the extremes", {
  z <- stats::qnorm(0.95)
  expect_true(wilson_lower(1, 1, z) > 0 && wilson_lower(1, 1, z) < 1)
  expect_true(wilson_lower(0, 1, z) >= 0 && wilson_lower(0, 1, z) < 1)
  expect_equal(wilson_lower(0, 1, z), 0)
})

test_that("wilson_lower is always <= phat (a valid interval contains its own point estimate)", {
  z <- stats::qnorm(0.95)
  for (phat in c(0, 0.1, 0.5, 0.8, 0.95, 1)) {
    for (n in c(1, 5, 50, 1000)) {
      expect_lte(wilson_lower(phat, n, z), phat + 1e-12)
    }
  }
})

test_that("wilson_lower converges to the Wald bound as n -> Inf", {
  z <- stats::qnorm(0.95)
  phat <- 0.8
  wald <- function(n) phat - z * sqrt(phat * (1 - phat) / n)
  gap <- function(n) abs(wilson_lower(phat, n, z) - wald(n))
  expect_true(gap(1e5) < gap(1e3))
  expect_true(gap(1e3) < gap(10))
})

test_that("wilson_lower is vectorized over phat and n", {
  z <- stats::qnorm(0.95)
  out <- wilson_lower(c(0, 0.5, 1), c(10, 100, 1000), z)
  expect_length(out, 3)
  expect_true(all(out >= 0 & out <= 1))
})

# --- save_rng_state()/restore_rng_state() (internal) ---------------------

save_rng_state <- get("save_rng_state", envir = asNamespace("dtasamplesize"))
restore_rng_state <- get("restore_rng_state", envir = asNamespace("dtasamplesize"))

test_that("save_rng_state()/restore_rng_state() round-trip both the kind and the seed", {
  old_kind <- RNGkind()
  on.exit(suppressWarnings(
    RNGkind(kind = old_kind[1], normal.kind = old_kind[2], sample.kind = old_kind[3])
  ), add = TRUE)

  RNGkind("Mersenne-Twister")
  set.seed(123)
  invisible(runif(1))
  state <- save_rng_state()

  RNGkind("Wichmann-Hill")
  set.seed(999)
  invisible(runif(5))

  restore_rng_state(state)
  expect_identical(RNGkind()[1], "Mersenne-Twister")
  expect_identical(get(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                    state$seed)
})

test_that("save_rng_state()/restore_rng_state() handle a virgin session (.Random.seed absent)", {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    saved <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
    rm(".Random.seed", envir = .GlobalEnv)
  }
  state <- save_rng_state()
  expect_false(state$seed_present)

  invisible(runif(1))  # creates .Random.seed
  expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))

  restore_rng_state(state)
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})
