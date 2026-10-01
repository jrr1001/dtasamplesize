# G04 -- exact-mode invariance to B, seed and N_range
# --------------------------------------------------------------------------
# bam_sample_size(method = "exact") computes its headline N_total and
# joint_assurance from a closed-form Beta-Binomial enumeration
# (.bam_exact_joint_assurance(), see R/bam_sample_size.R) with NO Monte
# Carlo draws anywhere on that path. Consequently N_total and
# joint_assurance for a fixed set of priors/deltas/target must be
# IDENTICAL no matter what B or seed the caller passes (B/seed only affect
# legacy per-arm diagnostic fields that are not part of the headline
# result), and no matter whether N_range is left NULL (automatic integer
# scan) or supplied explicitly as a range that contains the crossing. This
# is the G04 gate in CORRECCION_02_COMPUERTAS_VERIFICACION.md section 7.
#
# Two published cases are checked:
#   - "678" / vague Sp prior (Figure 2 worked example): prior_se=c(17,3),
#     prior_sp=c(2,2), prior_prev=c(4,16), delta_se=0.14, delta_sp=0.10,
#     target_assurance=0.80 -> N_total=678, joint_assurance=0.8003489948.
#   - "672" / harmonized Sp prior (Table 4): prior_se=c(17,3),
#     prior_sp=c(18,2), prior_prev=c(4,16), delta_se=0.14, delta_sp=0.10,
#     target_assurance=0.80 -> N_total=672, joint_assurance=0.8002692084.
# (Both verbatim from validation/reproduce_manuscript.R sections 2 and 3.)
#
# Sweep: B in {5, 10, 50, 5000} x seed in {1, 2, 3, 4, 2026} x N_range in
# {NULL, an explicit range containing the crossing}. All combinations for a
# given case must return identical N_total and joint_assurance, and must
# leave the caller's .Random.seed and RNGkind() exactly as they were before
# the call.
#
# A separate block checks the no-crossing contract: an explicit N_range
# that does NOT contain the crossing (e.g. entirely below it) must return
# target_reached = FALSE and N_total = NA_integer_ -- never
# max(N_range) (the pre-0.6.6 defect this gate exists to catch).
#
# VERSION-CHECK OVERRIDE (for demonstrating the gate against v0.6.5,
# where it is expected to FAIL): pass --version-check=off on the command
# line, e.g.
#   Rscript validation/exact_invariance.R --version-check=off
# This does not weaken the normal guard; it only lets the script run
# against a non-0.6.6 installation so G04 can show the check actually
# detects the pre-0.6.6 defect rather than passing vacuously everywhere.
#
# Run from the package root (against the 0.6.6 private library):
#   Rscript validation/exact_invariance.R 2>&1 | tee validation/logs/exact_invariance.log
# Run against the 0.6.5 library (expected to FAIL; log kept separately):
#   R_LIBS=<path to lib065> Rscript validation/exact_invariance.R --version-check=off \
#     2>&1 | tee validation/logs/exact_invariance_v065.log

library(dtasamplesize)
options(width = 110, digits = 10)

args <- commandArgs(trailingOnly = TRUE)
version_check_off <- any(grepl("^--version-check=off$", args))

installed_version <- utils::packageVersion("dtasamplesize")
if (!version_check_off && installed_version != "0.6.6") {
  stop("dtasamplesize ", installed_version, " is installed, but this script ",
       "checks the exact-mode invariance contract introduced in 0.6.6. Pass ",
       "--version-check=off to run anyway (e.g. to demonstrate the gate ",
       "FAILS against an earlier version).", call. = FALSE)
}
cat("dtasamplesize installed version:", as.character(installed_version), "\n")
cat("--version-check=off:", version_check_off, "\n\n")

sep <- function(t) cat("\n==========", t, "==========\n")
ok  <- function(x) if (isTRUE(x)) "PASS" else "FAIL"

Bs    <- c(5, 10, 50, 5000)
seeds <- c(1, 2, 3, 4, 2026)

cases <- list(
  "678" = list(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    N_total_published = 678L, joint_assurance_published = 0.8003489948,
    N_range_explicit = 600:700
  ),
  "672" = list(
    prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    N_total_published = 672L, joint_assurance_published = 0.8002692084,
    N_range_explicit = 600:700
  )
)

## ---- RNG snapshot / restore helpers (caller state must be untouched) ----
rng_snapshot <- function() {
  list(
    kind = RNGkind(),
    seed_present = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
    seed = if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      get(".Random.seed", envir = .GlobalEnv) else NULL
  )
}
rng_identical <- function(a, b) {
  identical(a$kind, b$kind) &&
    identical(a$seed_present, b$seed_present) &&
    identical(a$seed, b$seed)
}
restore_rng <- function(snap) {
  suppressWarnings(RNGkind(kind = snap$kind[1], normal.kind = snap$kind[2],
                            sample.kind = snap$kind[3]))
  if (snap$seed_present) {
    assign(".Random.seed", snap$seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}
original_rng <- rng_snapshot()
on.exit(restore_rng(original_rng), add = TRUE)

## ---- 1. Sweep over B x seed x N_range for each published case ----
all_pass <- TRUE
mismatch_msgs <- character(0)

for (case_name in names(cases)) {
  cs <- cases[[case_name]]
  sep(paste0("Case ", case_name, ": B x seed x N_range sweep"))

  results <- list()
  for (B in Bs) {
    for (seed in seeds) {
      for (nr_label in c("NULL", "explicit")) {
        key <- paste(B, seed, nr_label, sep = "|")
        nr <- if (nr_label == "NULL") NULL else cs$N_range_explicit

        before <- rng_snapshot()
        set.seed(123456)  # arbitrary caller state, same for every call
        caller_before <- rng_snapshot()

        r <- suppressWarnings(bam_sample_size(
          prior_se = cs$prior_se, prior_sp = cs$prior_sp,
          prior_prev = cs$prior_prev, delta_se = cs$delta_se,
          delta_sp = cs$delta_sp, target_assurance = cs$target_assurance,
          method = "exact", N_range = nr, B = B, seed = seed))

        caller_after <- rng_snapshot()
        rng_ok <- rng_identical(caller_before, caller_after)

        results[[key]] <- list(
          B = B, seed = seed, N_range = nr_label,
          N_total = r$N_total, joint_assurance = r$joint_assurance,
          rng_preserved = rng_ok
        )
        restore_rng(before)
      }
    }
  }

  N_vals <- vapply(results, function(x) x$N_total, integer(1))
  A_vals <- vapply(results, function(x) x$joint_assurance, double(1))
  rng_vals <- vapply(results, function(x) x$rng_preserved, logical(1))

  N_identical <- length(unique(N_vals)) == 1L
  A_identical <- length(unique(A_vals)) == 1L
  N_matches_published <- N_identical && unique(N_vals) == cs$N_total_published
  A_matches_published <- A_identical &&
    isTRUE(all.equal(unique(A_vals), cs$joint_assurance_published, tolerance = 1e-8))
  rng_all_preserved <- all(rng_vals)

  cat(sprintf("  combinations tested: %d (B in {%s} x seed in {%s} x N_range in {NULL, explicit})\n",
              length(results), paste(Bs, collapse = ","), paste(seeds, collapse = ",")))
  cat(sprintf("  N_total identical across all combinations: %s (value: %s, published: %d)   %s\n",
              N_identical, if (N_identical) unique(N_vals) else "VARIES", cs$N_total_published,
              ok(N_identical)))
  cat(sprintf("  N_total matches published value: %s\n", ok(N_matches_published)))
  cat(sprintf("  joint_assurance identical across all combinations: %s (value: %s, published: %.10f)   %s\n",
              A_identical, if (A_identical) format(unique(A_vals), digits = 12) else "VARIES",
              cs$joint_assurance_published, ok(A_identical)))
  cat(sprintf("  joint_assurance matches published value: %s\n", ok(A_matches_published)))
  cat(sprintf("  caller RNG state (.Random.seed + RNGkind) preserved in all %d calls: %s\n",
              length(results), ok(rng_all_preserved)))

  case_pass <- N_identical && A_identical && N_matches_published &&
    A_matches_published && rng_all_preserved
  if (!case_pass) {
    all_pass <- FALSE
    mismatch_msgs <- c(mismatch_msgs, sprintf(
      "Case %s: N_identical=%s A_identical=%s N_matches_published=%s A_matches_published=%s rng_preserved=%s",
      case_name, N_identical, A_identical, N_matches_published, A_matches_published, rng_all_preserved))
  }
}

## ---- 2. No-crossing contract: insufficient N_range -> NA, not max(N_range) ----
sep("No-crossing contract: N_range below the crossing must return NA, not max(N_range)")

cs <- cases[["678"]]
insufficient_range <- 10:20  # far below any crossing for these priors/deltas
r_insufficient <- suppressWarnings(bam_sample_size(
  prior_se = cs$prior_se, prior_sp = cs$prior_sp, prior_prev = cs$prior_prev,
  delta_se = cs$delta_se, delta_sp = cs$delta_sp,
  target_assurance = cs$target_assurance, method = "exact",
  N_range = insufficient_range, B = 50, seed = 1))

target_reached_false <- isTRUE(identical(r_insufficient$target_reached, FALSE))
n_total_is_na <- is.na(r_insufficient$N_total)
not_max_n_range <- !identical(r_insufficient$N_total, max(insufficient_range))

cat(sprintf("  target_reached == FALSE: %s (value: %s)   %s\n",
            target_reached_false, r_insufficient$target_reached, ok(target_reached_false)))
cat(sprintf("  N_total is NA: %s (value: %s)   %s\n",
            n_total_is_na, format(r_insufficient$N_total), ok(n_total_is_na)))
cat(sprintf("  N_total is NOT max(N_range) (%d): %s\n",
            max(insufficient_range), ok(not_max_n_range)))

no_crossing_pass <- target_reached_false && n_total_is_na && not_max_n_range
if (!no_crossing_pass) {
  all_pass <- FALSE
  mismatch_msgs <- c(mismatch_msgs, "No-crossing contract failed (see above).")
}

## ---- Verdict ----
sep("Verdict")
if (all_pass) {
  cat("RESULTADO G04: PASS\n")
} else {
  cat("Failures:\n")
  for (m in mismatch_msgs) cat("  - ", m, "\n", sep = "")
  cat("RESULTADO G04: FAIL\n")
}

if (!all_pass) quit(status = 1, save = "no")
