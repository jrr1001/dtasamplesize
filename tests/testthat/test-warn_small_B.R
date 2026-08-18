
# The small-B advisory warning (warn_small_B() in R/helpers.R) is controlled
# by the 'dtasamplesize.warn_small_B' option. The suite-wide setup.R sets it
# to FALSE so the rest of the suite (which uses small B deliberately, for
# speed) stays quiet; here the option is flipped back to TRUE to verify the
# warning actually fires, and restored (via on.exit) after every test.

test_that("warn_small_B() itself warns below 1000 and not at/above it, when the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)

  expect_warning(warn_small_B(999), "B >= 1000")
  expect_warning(warn_small_B(1), "small")
  expect_no_warning(warn_small_B(1000))
  expect_no_warning(warn_small_B(5000))
  expect_no_warning(warn_small_B(0))  # B = 0 means "skip", not "small MC"
})

test_that("warn_small_B() is silent when the option is FALSE", {
  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(warn_small_B(100))
})

test_that("warn_small_B() is silent when the option is left at its default (TRUE) is opt-in only via the option", {
  # getOption(..., TRUE) means the warning is ON by default if the caller
  # never touches the option at all (i.e. outside of this package's own
  # test suite, which forces it off in setup.R).
  old <- getOption("dtasamplesize.warn_small_B")
  options(dtasamplesize.warn_small_B = NULL)
  on.exit(options(dtasamplesize.warn_small_B = old), add = TRUE)
  expect_warning(warn_small_B(100), "small")
})

# --- exported functions: warning present/absent per the option ---------

test_that("mc_validate_buderer warns on small B iff the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(mc_validate_buderer(B = 100, seed = 1), "small")

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(mc_validate_buderer(B = 100, seed = 1))
})

test_that("bam_sample_size warns on small B iff the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(bam_sample_size(B = 100, seed = 1, n_range = 20:600), "small")

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(bam_sample_size(B = 100, seed = 1, n_range = 20:600))
})

test_that("joint_sample_size warns on small B iff the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(
    joint_sample_size(B = 100, seed = 1, N_range = seq(100, 800, by = 50)),
    "small"
  )

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(
    joint_sample_size(B = 100, seed = 1, N_range = seq(100, 800, by = 50))
  )
})

test_that("ss_net_benefit warns on small B iff the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(ss_net_benefit(pt_range = 0.20, B = 100, seed = 1), "small")

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(ss_net_benefit(pt_range = 0.20, B = 100, seed = 1))
})

test_that("ss_adaptive_prevalence warns on small B iff the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(
    ss_adaptive_prevalence(B = 100, seed = 1, prev_true_range = 0.30),
    "small"
  )

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(
    ss_adaptive_prevalence(B = 100, seed = 1, prev_true_range = 0.30)
  )
})

test_that("ss_imperfect_ref warns on small B iff the option is TRUE, and never when B = 0", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(ss_imperfect_ref(B = 100, seed = 1), "small")
  expect_no_warning(ss_imperfect_ref(B = 0))  # B = 0 skips MC validation entirely

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(ss_imperfect_ref(B = 100, seed = 1))
})

test_that("ss_unified warns on small B iff the option is TRUE", {
  options(dtasamplesize.warn_small_B = TRUE)
  on.exit(options(dtasamplesize.warn_small_B = FALSE), add = TRUE)
  expect_warning(
    ss_unified(B = 100, seed = 1, N_range = seq(200, 1000, by = 100),
              delta_auc = 0, check_nb = FALSE),
    "small"
  )

  options(dtasamplesize.warn_small_B = FALSE)
  expect_no_warning(
    ss_unified(B = 100, seed = 1, N_range = seq(200, 1000, by = 100),
              delta_auc = 0, check_nb = FALSE)
  )
})
