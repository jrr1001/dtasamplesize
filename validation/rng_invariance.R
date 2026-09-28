# RNG-kind invariance check for the manuscript's headline numbers
# --------------------------------------------------------------
# The manuscript claims: "every published number in this article was
# regenerated under all five random-number generators available in R
# (Mersenne-Twister, L'Ecuyer-CMRG, Wichmann-Hill, Marsaglia-Multicarry, and
# Knuth-TAOCP-2002) and returned identical results in every case". No script
# in the repository previously generated this evidence; this script does.
#
# NOTE ON "available in R": R's documented uniform-generator kinds (see
# ?RNGkind) are Wichmann-Hill, Marsaglia-Multicarry, Super-Duper,
# Mersenne-Twister, Knuth-TAOCP, Knuth-TAOCP-2002, L'Ecuyer-CMRG, and
# user-supplied -- eight kinds, not five. The five named in the manuscript
# (and tested below) are therefore a SELECTION of the kinds R documents, not
# literally "all" of them. This script does not silently relabel the claim;
# it tests exactly the five named, and the header above says so explicitly
# so the discrepancy between "all five ... available in R" and "R documents
# eight" is visible to the reader rather than papered over.
#
# WHY THIS SHOULD PASS TRIVIALLY (recorded here BEFORE running, so the
# result cannot be read as more informative than it is): every function
# exercised below reseeds its own RNG stream internally with
# set.seed(seed, kind = "Mersenne-Twister", ...) before drawing anything,
# and saves/restores the CALLER's RNG kind and seed around that internal
# reseed (see save_rng_state()/restore_rng_state() in R/helpers.R, used by
# each function -- file:line citations in the header comment blocks below).
# One function (buderer_n()) contains no random draws at all: it is a
# closed-form ceiling() of a formula, so it cannot depend on the RNG kind
# in any way. And bam_sample_size(method = "exact") computes its headline
# N_total / joint_assurance from a closed-form Beta-Binomial enumeration
# with no sampling step at all (the function's only internal set.seed()
# calls feed legacy per-arm diagnostic fields that are not part of the
# published N = 678 / N = 672 numbers). Consequently, changing the CALLING
# session's RNGkind() before invoking these functions is not expected to
# change any published number -- not because the underlying computation is
# somehow robust to different random streams, but because the functions do
# not consume the caller's stream for the numbers being checked at all.
# This is worth stating plainly: the check below demonstrates exactly what
# it can demonstrate (caller-RNG-state independence, and caller-state
# preservation), not "the simulation is insensitive to which generator
# produced its draws" in some stronger sense.
#
# Internal-RNG-handling table (function : file:line of the relevant code):
#   buderer_n()                     -- R/helpers.R:139-142
#     No RNG use whatsoever: n <- ceiling(z^2 * p * (1 - p) / d^2).
#   mc_validate_buderer()           -- R/mc_validate_buderer.R
#     save_rng_state()/on.exit(restore_rng_state()): lines 63-64.
#     set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
#       sample.kind = "Rejection"): line 83.
#   bam_sample_size()               -- R/bam_sample_size.R
#     save_rng_state()/on.exit(restore_rng_state()): lines 329-330.
#     Legacy per-arm set.seed(..., kind = "Mersenne-Twister", ...) calls
#       (n_se / n_sp / N_total legacy-heuristic searches, NOT the headline
#       exact N_total/joint_assurance): lines 384, 405, 420, 439, 454, 490.
#     method = "exact" headline branch (N_total, joint_assurance for the
#       published N = 678 and N = 672 figures): lines 542-593 -- closed-form
#       Beta-Binomial enumeration via .bam_exact_joint_assurance(), no
#       set.seed() call anywhere in this branch, no B/seed dependence
#       (confirmed by in-code comment at lines 542-544).
#   joint_sample_size()             -- R/joint_sample_size.R
#     save_rng_state()/on.exit(restore_rng_state()): lines 194-195.
#     set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion",
#       sample.kind = "Rejection"): line 261 (inside the per-candidate-N
#       Monte Carlo search over Se/Sp coverage).
#
# Run from the package root:
#   Rscript validation/rng_invariance.R 2>&1 | tee validation/logs/rng_invariance.log

library(dtasamplesize)
options(width = 110, digits = 10)

if (utils::packageVersion("dtasamplesize") != "0.6.5") {
  stop("dtasamplesize ", utils::packageVersion("dtasamplesize"), " is installed, ",
       "but this script checks version 0.6.5's published numbers exactly.",
       call. = FALSE)
}

SEED <- 2026

cat("R.version.string:", R.version.string, "\n")
cat("RNGkind() documented kinds (per ?RNGkind, 'all five ... available in R' in\n",
    "the manuscript is a selection of these eight, not literally all of them):\n",
    "  Wichmann-Hill, Marsaglia-Multicarry, Super-Duper, Mersenne-Twister,\n",
    "  Knuth-TAOCP, Knuth-TAOCP-2002, L'Ecuyer-CMRG, user-supplied\n", sep = "")
cat("Kinds tested here (the five the manuscript names):\n",
    "  Mersenne-Twister, L'Ecuyer-CMRG, Wichmann-Hill, Marsaglia-Multicarry,\n",
    "  Knuth-TAOCP-2002\n\n", sep = "")

kinds <- c("Mersenne-Twister", "L'Ecuyer-CMRG", "Wichmann-Hill",
           "Marsaglia-Multicarry", "Knuth-TAOCP-2002")

# Restore the session's original RNG kind (and seed presence) no matter how
# this script exits.
original_rng_state <- RNGkind()
original_seed_present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
original_seed <- if (original_seed_present) get(".Random.seed", envir = .GlobalEnv) else NULL
restore_original_rng <- function() {
  suppressWarnings(RNGkind(kind = original_rng_state[1],
                            normal.kind = original_rng_state[2],
                            sample.kind = original_rng_state[3]))
  if (original_seed_present) {
    assign(".Random.seed", original_seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
}
on.exit(restore_original_rng(), add = TRUE)

## ---- Section 3 constants (harmonized priors), verbatim from
## validation/reproduce_manuscript.R ----
PRIOR_SE <- c(17, 3); PRIOR_SP_HARMONIZED <- c(18, 2); PRIOR_PREV <- c(4, 16)
DELTA_SE <- 0.07; DELTA_SP <- 0.05; B_MAIN <- 20000
E_SE <- PRIOR_SE[1] / sum(PRIOR_SE)
E_SP <- PRIOR_SP_HARMONIZED[1] / sum(PRIOR_SP_HARMONIZED)
E_PREV <- PRIOR_PREV[1] / sum(PRIOR_PREV)

compute_all <- function() {
  n_bud_085 <- buderer_n(0.85, 0.07)
  n_bud_090 <- buderer_n(0.90, 0.05)

  mv <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = 100,
                             B = 4000, ci_method = "wald", seed = SEED)$results

  bam_example <- bam_sample_size(
    prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
    delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
    method = "exact", N_range = 600:700, B = 5000, seed = SEED)

  N_buderer <- ceiling(max(buderer_n(0.85, 0.07) / 0.20,
                            buderer_n(0.90, 0.05) / 0.80))

  bam_harmonized <- suppressWarnings(bam_sample_size(
    prior_se = PRIOR_SE, prior_sp = PRIOR_SP_HARMONIZED, prior_prev = PRIOR_PREV,
    delta_se = 2 * DELTA_SE, delta_sp = 2 * DELTA_SP, target_assurance = 0.80,
    method = "exact", B = B_MAIN, seed = SEED))

  joint_res <- suppressWarnings(joint_sample_size(
    Se = E_SE, Sp = E_SP, delta_se = DELTA_SE, delta_sp = DELTA_SP, prev = E_PREV,
    design = "cohort", target_prob = 0.80, B = B_MAIN, seed = SEED))

  list(
    n_bud_085              = n_bud_085,
    n_bud_090               = n_bud_090,
    mc_P_width_target       = mv$P_width_target[1],
    N_buderer_combined      = N_buderer,
    bam_example_N_total     = bam_example$N_total,
    bam_example_assurance   = bam_example$joint_assurance,
    bam_harmonized_N_total  = bam_harmonized$N_total,
    bam_harmonized_assurance = bam_harmonized$joint_assurance,
    joint_N_total           = joint_res$n_total,
    joint_assurance         = joint_res$joint_prob_se_sp
  )
}

published <- list(
  n_bud_085                = 100L,
  n_bud_090                = 139L,
  mc_P_width_target         = 0.5730,
  N_buderer_combined        = 500L,
  bam_example_N_total       = 678L,
  bam_example_assurance     = 0.8003489948,
  bam_harmonized_N_total    = 672L,
  bam_harmonized_assurance  = 0.8002692084,
  joint_N_total             = 580L,
  joint_assurance           = 0.8041
)

results <- list()
rng_preserved <- list()

for (k in kinds) {
  # Set the CALLER's RNG kind and seed the generator under it, exactly as
  # instructed: the point is to prove the package's published numbers do
  # not depend on what generator the calling session happens to be using.
  suppressWarnings(RNGkind(kind = k))
  set.seed(1)

  results[[k]] <- compute_all()

  # Caller-state-preservation check: after each package call above, is the
  # CALLING session's RNGkind() still the kind we set it to, i.e. did the
  # package's internal save_rng_state()/restore_rng_state() correctly put
  # it back rather than leaking its internal Mersenne-Twister state out?
  rng_preserved[[k]] <- identical(RNGkind()[1], k)
}

## ---- Report: per-kind results table ----
sep <- function(t) cat("\n==========", t, "==========\n")
sep("Per-kind results (all numbers should be identical across rows)")

metric_names <- names(published)
tab <- do.call(rbind, lapply(kinds, function(k) {
  unlist(results[[k]][metric_names])
}))
rownames(tab) <- kinds
print(tab)

cat("\nPublished (manuscript) values, for reference:\n")
print(unlist(published))

## ---- Report: caller RNG-state preservation ----
sep("Caller RNG-state preservation (RNGkind()[1] after all calls, under each kind)")
preserved_tab <- data.frame(
  kind_set     = kinds,
  RNGkind_after = vapply(kinds, function(k) RNGkind_after <- {
    # recompute for display purposes is unnecessary; use stored logical
    if (isTRUE(rng_preserved[[k]])) k else paste0("MISMATCH (", k, ")")
  }, character(1)),
  preserved    = unlist(rng_preserved),
  stringsAsFactors = FALSE
)
print(preserved_tab, row.names = FALSE)

## ---- Compare every kind's results against Mersenne-Twister with identical() ----
sep("identical() comparison against Mersenne-Twister")

reference <- results[["Mersenne-Twister"]]
all_identical <- TRUE
mismatches <- character(0)

for (k in kinds) {
  if (k == "Mersenne-Twister") next
  for (m in metric_names) {
    is_id <- identical(reference[[m]], results[[k]][[m]])
    cat(sprintf("  %-22s vs Mersenne-Twister | %-26s | %s\n",
                k, m, if (is_id) "identical" else "**DIFFERS**"))
    if (!is_id) {
      all_identical <- FALSE
      mismatches <- c(mismatches, sprintf(
        "%s / %s: Mersenne-Twister = %s, %s = %s",
        k, m, format(reference[[m]], digits = 12),
        k, format(results[[k]][[m]], digits = 12)))
    }
  }
}

## ---- Compare every kind's results against the published manuscript values ----
sep("Comparison against published manuscript values")

published_mismatches <- character(0)
for (k in kinds) {
  for (m in metric_names) {
    got <- results[[k]][[m]]
    want <- published[[m]]
    ok_num <- if (is.integer(want) || (is.numeric(want) && want == round(want) && m %in% c("n_bud_085","n_bud_090","N_buderer_combined","bam_example_N_total","bam_harmonized_N_total","joint_N_total"))) {
      got == want
    } else {
      isTRUE(abs(got - want) <= 5e-4)
    }
    if (!isTRUE(ok_num)) {
      published_mismatches <- c(published_mismatches, sprintf(
        "%s / %s: computed = %s, published = %s",
        k, m, format(got, digits = 12), format(want, digits = 12)))
    }
  }
}
if (length(published_mismatches) == 0) {
  cat("All computed values match the published manuscript values under every kind tested.\n")
} else {
  cat("Computed values differing from the published manuscript values:\n")
  for (msg in published_mismatches) cat("  - ", msg, "\n", sep = "")
}

if (!all(unlist(rng_preserved))) {
  cat("\nCaller RNG state was NOT preserved for at least one kind (see table above).\n")
}

## ---- Verdict ----
sep("Verdict")
if (all_identical && all(unlist(rng_preserved)) && length(published_mismatches) == 0) {
  cat("PASS: all published numbers identical under 5 RNG kinds\n")
} else {
  msg <- character(0)
  if (!all_identical) {
    msg <- c(msg, "Numbers differed across RNG kinds:", paste0("  - ", mismatches))
  }
  if (!all(unlist(rng_preserved))) {
    msg <- c(msg, "Caller RNG state (RNGkind) was not preserved for at least one kind.")
  }
  if (length(published_mismatches) > 0) {
    msg <- c(msg, "Computed values differed from the published manuscript values:",
              paste0("  - ", published_mismatches))
  }
  stop(paste(msg, collapse = "\n"), call. = FALSE)
}
