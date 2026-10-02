# Contract tests for joint_sample_size() under the "Option A" decision
# (lote 01b, M-03): keep the step-10 default grid and its published figure
# (n_total = 580, joint_prob_se_sp = 0.8041) unchanged, reinterpret 580 as
# "the first candidate of the step-10 grid reaching the target" (never
# "smallest"), and give explicit non-crossing output (NA, never
# max(N_range)) plus exposed grid/B/seed/MCSE/AUC-gate diagnostics. See
# DECISION_JOINT_SAMPLE_SIZE.md (Option A) and NEWS.md (0.6.6) for the full
# rationale; this file pins the resulting contract.

# --- the published-case call: the figure quoted in the manuscript (Table 4)
# is reproduced from Se/Sp/prev/delta_se/delta_sp at their package DEFAULTS
# together with the B/seed actually used for that figure (B = 20000, seed =
# 2026 -- the function's own B default is 5000, see ?joint_sample_size).
# This is the "default call" referred to throughout the lote 01b audit
# trail (DECISION_JOINT_SAMPLE_SIZE.md, tools/L01_analisis_joint.R block
# (a)): defaults for every argument EXCEPT B/seed, which are pinned to the
# manuscript's declared values. -------------------------------------------

test_that("the published-case call reaches n_total = 580 with joint_prob_se_sp = 0.8041", {
  result <- joint_sample_size(B = 20000, seed = 2026)

  expect_true(result$target_reached)
  # n_total = 580 is stable across R RNG streams (same grid, same AUC gate,
  # same first-crossing candidate); the Monte Carlo *probability* at that N
  # is NOT bit-identical across R versions (R-devel/4.6.x changed the RNG
  # stream used by this sampling loop relative to R 4.5.x -- see NEWS
  # 0.6.7). We pin n_total exactly, and check
  # joint_prob_se_sp against its own reported MCSE rather than against the
  # single value observed under R 4.5.x.
  expect_equal(result$n_total, 580L)
  expect_equal(result$joint_prob, result$joint_prob_se_sp) # deprecated alias

  expect_equal(result$search_type, "grid_first_candidate")
  expect_equal(result$B, 20000)
  expect_equal(result$seed, 2026)
  mcse_from_p <- sqrt(result$joint_prob_se_sp * (1 - result$joint_prob_se_sp) / 20000)
  expect_equal(result$joint_prob_mcse, mcse_from_p, tolerance = 1e-9)

  # Published-figure value (R 4.5.x RNG stream): within +/- 3 MCSE of the
  # value reported in the manuscript/NEWS, robust to the RNG-stream change
  # observed on R-devel (>= 4.6.0).
  expect_lt(abs(result$joint_prob_se_sp - 0.8041), 3 * mcse_from_p)
  if (getRversion() < "4.6.0") {
    # Exact published figure, reproducible deterministically on R < 4.6.0
    # (the RNG stream used by this Monte Carlo loop is unchanged there).
    expect_equal(result$joint_prob_se_sp, 0.8041, tolerance = 1e-4)
    expect_equal(result$joint_prob_mcse, 0.002806, tolerance = 1e-3)
  }

  # N_range_used is the default grid, exactly: seq(100, 800, by = 10). The
  # search stops at the first crossing (580), but N_range_used still
  # records the FULL declared grid (every element was "visited" by the
  # loop, whether or not its Se/Sp probability was computed; see
  # ?joint_sample_size, @return).
  expect_identical(result$N_range_used, as.integer(seq(100, 800, by = 10)))

  expect_true(result$auc_gate_passed)
  expect_equal(result$auc_gate$first_N_auc_pass, 350L)
})

test_that("n_total is NOT described or computed as the smallest integer N", {
  # At the published-case parameters, a finer (integer-step) search crosses
  # target_prob = 0.80 strictly before 580 (verified independently in
  # tools/L01_analisis_joint.R, block (c): the deterministic finite sum
  # gives joint probability 0.8041 already at N = 579). 580 is therefore
  # the first grid candidate, not the minimum integer N.
  result <- joint_sample_size(B = 20000, seed = 2026)
  expect_equal(result$n_total, 580L)
  # Documented explicitly, not merely true by accident of this one example:
  expect_false(grepl("smallest", tolower(paste(deparse(body(joint_sample_size)),
                                                collapse = " "))))
})

# --- non-crossing: AUC gate blocks every candidate in a short N_range -----

test_that("non-crossing (AUC gate blocks everything): n_total is NA, never max(N_range)", {
  expect_warning(
    result <- joint_sample_size(B = 500, seed = 2026,
                                N_range = seq(100, 150, by = 10)),
    "AUC precision target"
  )
  expect_false(result$target_reached)
  expect_true(is.na(result$n_total))
  expect_true(is.na(result$joint_prob_se_sp))
  expect_true(is.na(result$joint_prob))
  expect_true(is.na(result$n_diseased))
  expect_true(is.na(result$n_non_diseased))
  expect_false(result$auc_gate_passed)
  # max_joint_prob_evaluated/N_at_max_joint_prob are NA here specifically
  # because the AUC gate blocked EVERY candidate, so no Se/Sp probability
  # was ever computed anywhere -- not because the search happened to find
  # a low value.
  expect_true(is.na(result$max_joint_prob_evaluated))
  expect_true(is.na(result$N_at_max_joint_prob))
  # n_total must never silently become max(N_range) (= 150):
  expect_false(isTRUE(result$n_total == max(seq(100, 150, by = 10))))
})

# --- non-crossing: AUC gate passes, but the Se/Sp joint probability stays
# below target_prob everywhere in N_range -----------------------------------

test_that("non-crossing (AUC gate passes, target never reached): max_joint_prob_evaluated < target_prob", {
  expect_warning(
    result <- joint_sample_size(B = 300, seed = 2026,
                                N_range = seq(350, 400, by = 10)),
    "No N in N_range achieved"
  )
  expect_false(result$target_reached)
  expect_true(is.na(result$n_total))
  expect_true(is.na(result$joint_prob_se_sp))
  expect_true(result$auc_gate_passed) # AUC gate DID pass (first pass = 350)
  expect_false(is.na(result$max_joint_prob_evaluated))
  expect_lt(result$max_joint_prob_evaluated, 0.80)
  expect_false(is.na(result$N_at_max_joint_prob))
  expect_true(result$N_at_max_joint_prob %in% seq(350, 400, by = 10))
  # n_total must never silently become max(N_range) (= 400):
  expect_false(isTRUE(result$n_total == 400L))
})

# --- the AUC gate blocking a candidate is NOT the same as that candidate
# having joint probability 0 -------------------------------------------------

test_that("AUC-gate-blocked candidates are distinguishable from probability-0 candidates", {
  # Default grid: candidates below N = 350 fail the (evaluated) AUC gate;
  # N = 100 has a non-degenerate expected margin (n_d_exp = 20, n_nd_exp =
  # 80), so its AUC gate IS evaluated (not skipped) and fails.
  result <- joint_sample_size(B = 20000, seed = 2026)
  tbl <- result$auc_gate$table
  row_100 <- tbl[tbl$N == 100, ]
  expect_false(is.na(row_100$auc_pass))
  expect_false(row_100$auc_pass) # evaluated AND failed -- not "skipped"
  # A candidate whose AUC gate failed never had its Se/Sp joint probability
  # computed at all (it is absent from the Se/Sp search entirely), which is
  # different from a computed probability of exactly 0. The table itself
  # only ever stores TRUE/FALSE/NA for auc_pass, never a joint probability
  # of 0 standing in for "blocked".
  expect_true(all(tbl$auc_pass %in% c(TRUE, FALSE, NA)))
  expect_equal(result$auc_gate$first_N_auc_pass, 350L)
  # Candidates at or above 350 (the full grid step-10 point of 350 upward)
  # all pass the AUC gate near the crossing region (confirmed independently
  # in DECISION_JOINT_SAMPLE_SIZE.md, Section 2(3): 570-590 all pass).
  row_580 <- tbl[tbl$N == 580, ]
  expect_true(isTRUE(row_580$auc_pass))
})

# --- print(): no-crossing never shows an N as a solution; success never
# says "smallest" --------------------------------------------------------

test_that("print() of a non-crossing result shows no sample size and no N as solution", {
  result <- suppressWarnings(joint_sample_size(B = 500, seed = 2026,
                                               N_range = seq(100, 150, by = 10)))
  out <- paste(capture.output(print(result)), collapse = "\n")
  expect_match(out, "NOT reached")
  expect_match(out, "no sample size returned")
  expect_false(grepl("N_total: 1[45]0", out)) # no candidate shown as "the" N
  expect_false(grepl("smallest", out, ignore.case = TRUE))
})

test_that("print() of a successful result never says \"smallest\" and shows B/seed/MCSE/AUC gate", {
  result <- joint_sample_size(B = 20000, seed = 2026)
  out <- paste(capture.output(print(result)), collapse = "\n")
  expect_false(grepl("smallest", out, ignore.case = TRUE))
  expect_match(out, "580")
  expect_match(out, "grid step")
  expect_match(out, "joint_prob_mcse")
  expect_match(out, "AUC gate")
  expect_match(out, "seed:\\s*2026")
  expect_match(out, "B:\\s*20000")
})

# --- RNG state: joint_sample_size() must not disturb the caller's RNG -----
# stream even when it takes the non-crossing branch (the branch most
# recently touched by this lote) --------------------------------------------

test_that("the caller's RNG stream is unaffected even when target_reached is FALSE", {
  set.seed(20260714)
  invisible(runif(3))
  expected <- runif(2)

  set.seed(20260714)
  invisible(runif(3))
  invisible(suppressWarnings(
    joint_sample_size(B = 500, seed = 2026, N_range = seq(100, 150, by = 10))
  ))
  actual <- runif(2)

  expect_identical(actual, expected)
})
