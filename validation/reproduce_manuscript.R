# Reproduce every headline number reported in the dtasamplesize manuscript
# ------------------------------------------------------------------------
# This script regenerates, from scratch, each quantity quoted in the article,
# using the same parameters and the same random seed. Every number the article
# reports is produced here; nothing is read from a stored dataset, because the
# study analyses no empirical data: all data are simulated (or, for the exact
# Bayesian assurance calculation, computed in closed form) inside the
# package's functions from the stated parameters and are fully determined by
# the seed.
#
# SCOPE (article now reduced to a single joint estimand). The article
# describes only mc_validate_buderer(), the exact joint assurance evaluated
# by bam_sample_size(method = "exact"), and joint_sample_size(); ss_unified(),
# the AUC gate inside it, ss_net_benefit(), ss_imperfect_ref(),
# ss_adaptive_prevalence() and ss_time_dependent_roc() remain in the package
# but are no longer described in this article, so this script no longer
# checks their numbers (the nested-sequence Figure 4 of v0.6.9, the net-benefit
# feasibility ceiling, the imperfect-reference apparent/corrected estimates,
# and the fixed-margin-vs-cohort variance illustration built around net
# benefit have all been removed from this script accordingly).
#
# WHAT THIS SCRIPT DOES NOT DO. Every check below is a REPRODUCIBILITY
# check: it confirms that the installed package, run today, prints the same
# numbers the article prints. It is NOT, and cannot be, a VALIDITY check.
# A defect that biases a result reproduces just as cleanly as a result that
# does not -- a PASS below only means the code and the article agree with
# each other, not that either is statistically correct. This distinction is
# not academic: the numbers this script currently checks against were
# themselves revised twice (0.3.0, then 0.5.0 and 0.6.0) after earlier
# releases passed their own reproducibility checks against defective
# numbers. Validity is established separately, release by release, in
# NEWS.md and by independent review of the methods -- never by this script
# agreeing with itself.
#
# Run after installing the package:
#   Rscript validation/reproduce_manuscript.R [output_dir]
# The output directory is the first command-line argument; if none is given,
# validation/output/ inside the package tree (created if needed). The script
# writes ONLY there (a CSV with the result of every check, and, if DUMP_DATA
# is TRUE, a simulated_data/ sub-directory) and never to folders above the
# package.
#
# Optional: set DUMP_DATA <- TRUE below to also write the underlying simulated
# samples to CSV files in <output_dir>/simulated_data/ (see "Simulated data"
# at the end). The CSVs are a convenience: they are regenerated exactly by
# this script, so the script itself is the authoritative record of the data.

library(dtasamplesize)
options(width = 110, digits = 6)

# The numbers checked below are those of the R1 manuscript. The computational
# code that produces them is unchanged from 0.6.9 to 0.7.0 for the default
# arguments (0.7.0 adds arguments and functions -- bam_assurance_curve(),
# bam_sensitivity(), analysis priors -- that leave the default numerical
# paths bit-identical; see NEWS.md), so this script requires dtasamplesize
# 0.7.0 or later. Earlier releases (historical only; not presented here as
# current validation evidence) contain the validity errors that 0.3.0, 0.5.0
# and 0.6.0 correct in turn (see NEWS.md) -- notably a net-benefit sample
# size of 110 instead of 240, a BAM search that reported the marginal per-arm
# assurance instead of the joint one, an adaptive-prevalence design that
# discarded its own stage-1 pilot, and an imperfect-reference sample size
# inflated by a formula that belongs to a different estimand -- so earlier
# installations return different values, and 0.6.9 lacks the functions this
# script calls (the exact mode of mc_validate_buderer() among them). Stop early
# rather than let a stale installation look like a discrepancy in the article.
if (utils::packageVersion("dtasamplesize") < "0.7.0") {
  stop("dtasamplesize ", utils::packageVersion("dtasamplesize"), " is installed, ",
       "but this script needs version 0.7.0 or later (the R1 manuscript's ",
       "numbers and the exact mode of mc_validate_buderer()).\n",
       "  Install it from a local copy of the repository:\n",
       "    install.packages('.', repos = NULL, type = 'source')",
       call. = FALSE)
}

## ---- output directory: first argument, default validation/output/ ----
locate_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (basename(getwd()) == "validation") return(normalizePath(getwd()))
  if (dir.exists("validation")) return(normalizePath(file.path(getwd(), "validation")))
  stop("Cannot determine the location of reproduce_manuscript.R. Run it via\n",
       "  Rscript validation/reproduce_manuscript.R [output_dir]",
       call. = FALSE)
}
script_dir <- locate_script_dir()
cli_args   <- commandArgs(trailingOnly = TRUE)
out_dir    <- if (length(cli_args) >= 1 && nzchar(cli_args[1])) cli_args[1] else
  file.path(script_dir, "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir    <- normalizePath(out_dir, winslash = "/")
cat("R version:          ", R.version.string, "\n")
cat("dtasamplesize:      ", as.character(utils::packageVersion("dtasamplesize")), "\n")
cat("Output directory:   ", out_dir, "\n")

DUMP_DATA <- FALSE          # set to TRUE to also export the simulated samples
SEED      <- 2026           # the seed cited throughout the manuscript

check_log <- data.frame(check = character(0), result = character(0),
                        stringsAsFactors = FALSE)
# The helper ok returns "PASS" / "**CHECK**" and records the labelled result. It
# compares VALUES (callers pass the logical result of a numeric comparison),
# never storage types: identical() on an integer and a double is a false
# failure.
ok  <- function(label, x) {
  res <- if (isTRUE(x)) "PASS" else "**CHECK**"
  check_log[nrow(check_log) + 1L, ] <<- list(label, res)
  res
}
sep <- function(t) cat("\n==========", t, "==========\n")
near <- function(a, b, tol) isTRUE(abs(a - b) <= tol)

## ---- 1. The motivating result (Results; Figure 1) ----
sep("1. Assurance gap of the classical Buderer sample size")
n_bud <- buderer_n(0.85, 0.07)
cat(sprintf("buderer_n(Se=0.85, d=0.07) = %d   (article: 100)   %s\n",
            n_bud, ok("buderer_n(0.85, 0.07) = 100", n_bud == 100)))
cat(sprintf("buderer_n(Se=0.90, d=0.05) = %d   (article: 139)   %s\n",
            buderer_n(0.90, 0.05), ok("buderer_n(0.90, 0.05) = 139", buderer_n(0.90, 0.05) == 139)))

mv <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = 100,
                          B = 4000, ci_method = "wald", seed = SEED)$results
cat(sprintf("P(CI width <= target) = %.4f   (Additional file 2, Table S2.7: 0.5730)   %s\n",
            mv$P_width_target[1], ok("MC P(width) n=100 = 0.5730", near(mv$P_width_target[1], 0.5730, 5e-4))))
cat(sprintf("Wald coverage         = %.4f   (Additional file 2, Table S2.7: 0.9380)   %s\n",
            mv$coverage[1], ok("MC Wald coverage n=100 = 0.9380", near(mv$coverage[1], 0.9380, 5e-4))))

# Exact counterpart (R1 Figure 1): the same quantity by binomial enumeration,
# no simulation. At the Buderer n both intervals give P = 0.5683.
mv_ex <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = 100,
                             ci_method = "both", method = "exact")$results
for (m in c("wald", "wilson")) {
  pm <- mv_ex$P_width_target[mv_ex$ci_method == m]
  cat(sprintf("exact P(CI width <= 0.14), %-6s = %.4f   (article: 0.5683)   %s\n",
              m, pm, ok(paste("exact P(100),", m, "= 0.5683"),
                        isTRUE(round(pm, 4) == 0.5683))))
}

## ---- 2. Joint assurance for Se and Sp, exact Beta-Binomial (Figure 2) ----
sep("2. Joint assurance for Se and Sp: exact crossing at N = 678 (Figure 2)")
# Article's worked-example priors -- Se ~ Beta(17,3) and Sp ~ Beta(2,2) (the
# package defaults) with prevalence ~ Beta(4,16) (NOT the package default
# prior_prev, which is c(6, 14)) -- and full-width targets delta_se = 0.14,
# delta_sp = 0.10, method = "exact" (closed-form Beta-Binomial, no Monte
# Carlo error, no dependence on B or seed for this headline result; see
# ?bam_sample_size). This is the article's central calculation (the
# "minimal session" code block in Methods) and the curve plotted in Figure 2.
bam_example <- bam_sample_size(
  prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
  delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
  method = "exact", N_range = 600:700, B = 5000, seed = SEED)
cat(sprintf("N_total = %d   (article: 678)   %s\n",
            bam_example$N_total, ok("BAM exact N_total = 678", bam_example$N_total == 678)))
cat(sprintf("joint assurance at N=678 = %.10f   (article: 0.8003489948)   %s\n",
            bam_example$joint_assurance,
            ok("BAM exact A(678) = 0.8003489948", near(bam_example$joint_assurance, 0.8003489948, 1e-8))))

# N = 677 is one below the crossing and does not reach the target on its
# own; bam_sample_size() warns that the search "did not achieve the target
# JOINT assurance" for a single-value N_range like this one -- expected and
# suppressed here. Since target_reached = FALSE for this call (by design --
# see ?bam_sample_size, @return), N_total/n_total/joint_assurance are all
# NA_real_/NA_integer_: joint_assurance is NOT the exact value at N = 677 in
# this case. The exact value at N = 677 is instead reported in
# max_assurance_evaluated (the highest exact/Monte-Carlo-estimated assurance
# seen over the N_range swept, together with N_at_max_assurance identifying
# which N achieved it) -- see ?bam_sample_size, @return. This check confirms
# both that target_reached is indeed FALSE here (so reading
# max_assurance_evaluated instead of joint_assurance is the correct call,
# not an oversight) and that max_assurance_evaluated at N = 677 matches the
# article.
bam_677 <- suppressWarnings(bam_sample_size(
  prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
  delta_se = 0.14, delta_sp = 0.10, target_assurance = 0.80,
  method = "exact", N_range = 677, B = 5000, seed = SEED))
cat(sprintf("target_reached at N=677 = %s   (expected: FALSE)   %s\n",
            bam_677$target_reached, ok("BAM exact target_reached(677) = FALSE", isTRUE(bam_677$target_reached == FALSE))))
cat(sprintf("max_assurance_evaluated at N=677 = %.10f   (article: 0.7996848824)   %s\n",
            bam_677$max_assurance_evaluated,
            ok("BAM exact A(677) = 0.7996848824", near(bam_677$max_assurance_evaluated, 0.7996848824, 1e-8))))

## ---- 3. Method comparison (Table 2 of v0.6.9 and of the R1): three surviving methods under common assumptions ----
sep("3. Method comparison (Table 2 of v0.6.9 and of the R1): Buderer, BAM (exact, harmonized priors), joint Se/Sp")
# Unlike section 2 above (the package's own default, vague Sp prior), the method comparison
# harmonizes the Bayesian method's priors against the classical and Monte
# Carlo methods' point assumptions: Se ~ Beta(17,3) (mean 0.85), Sp ~
# Beta(18,2) (mean 0.90, matching Se/Sp = 0.85/0.90 used elsewhere in the
# article), prevalence ~ Beta(4,16) (mean 0.20). This is why the same
# closed-form calculation returns N = 672 here against N = 678 in section 2:
# the informative Sp prior narrows the Sp credible interval faster than the
# vague default, so fewer subjects are needed for Sp to reach its target
# width, and the joint requirement drops accordingly.
PRIOR_SE <- c(17, 3); PRIOR_SP_HARMONIZED <- c(18, 2); PRIOR_PREV <- c(4, 16)
DELTA_SE <- 0.07; DELTA_SP <- 0.05; B_MAIN <- 20000
E_SE <- PRIOR_SE[1] / sum(PRIOR_SE)
E_SP <- PRIOR_SP_HARMONIZED[1] / sum(PRIOR_SP_HARMONIZED)
E_PREV <- PRIOR_PREV[1] / sum(PRIOR_PREV)

N_buderer <- ceiling(max(buderer_n(0.85, 0.07) / 0.20,
                          buderer_n(0.90, 0.05) / 0.80))
cat(sprintf("Buderer (classical), combined N = %d   (article: 500)   %s\n",
            N_buderer, ok("Buderer combined N = 500", N_buderer == 500)))

bam_harmonized <- suppressWarnings(bam_sample_size(
  prior_se = PRIOR_SE, prior_sp = PRIOR_SP_HARMONIZED, prior_prev = PRIOR_PREV,
  delta_se = 2 * DELTA_SE, delta_sp = 2 * DELTA_SP, target_assurance = 0.80,
  method = "exact", B = B_MAIN, seed = SEED))
cat(sprintf("BAM (exact, harmonized priors), N = %d   (article: 672)   %s\n",
            bam_harmonized$N_total, ok("BAM harmonized N = 672", bam_harmonized$N_total == 672)))
cat(sprintf("  joint assurance at N=672  = %.10f   (article: 0.8002692084)   %s\n",
            bam_harmonized$joint_assurance,
            ok("BAM harmonized A(672) = 0.8002692084", near(bam_harmonized$joint_assurance, 0.8002692084, 1e-8))))

# joint_sample_size() performs a grid search over N_range = seq(100, 800, by
# = 10) (the function's own default grid, unchanged by the 0.6.6 correction)
# and reports the FIRST CANDIDATE OF THAT 10-STEP GRID that reaches
# target_prob under B Monte Carlo replicates at the given seed -- it is not,
# and is not described here as, the smallest integer N that would reach the
# target; a finer grid or a different seed could return a different
# candidate with a different joint_prob_se_sp. This was decided in
# version 0.6.6 (see NEWS.md): N_range, B and seed are unchanged, and N = 580 is reported as "first candidate of the 10-step
# grid", never as "the minimum N" or "the smallest N".
joint_res <- suppressWarnings(joint_sample_size(
  Se = E_SE, Sp = E_SP, delta_se = DELTA_SE, delta_sp = DELTA_SP, prev = E_PREV,
  design = "cohort", target_prob = 0.80, B = B_MAIN, seed = SEED))
cat(sprintf("Joint Se/Sp (joint_sample_size), N = %d   (article: 580, first candidate of the 10-step grid)   %s\n",
            joint_res$n_total, ok("joint_sample_size N = 580", joint_res$n_total == 580)))
cat(sprintf("  joint assurance at N=580  = %.4f   (article: 0.8041)   %s\n",
            joint_res$joint_prob_se_sp,
            ok("joint_sample_size P = 0.8041", near(joint_res$joint_prob_se_sp, 0.8041, 5e-4))))
cat(sprintf("  joint_prob_mcse at N=580  = %.6f   (article/log: ~0.002806)   %s\n",
            joint_res$joint_prob_mcse,
            ok("joint_sample_size MCSE = 0.002806", near(joint_res$joint_prob_mcse, 0.002806, 5e-4))))
cat(sprintf("  B = %d, seed = %d, search_type = %s, N_range_used equal in value to seq(100,800,10): %s\n",
            B_MAIN, SEED, joint_res$search_type,
            ok("joint_sample_size N_range_used = seq(100,800,10) (values)",
               isTRUE(all.equal(as.numeric(joint_res$N_range_used),
                                seq(100, 800, by = 10))))))
cat(sprintf("  AUC gate: first_N_auc_pass = %s\n",
            format(joint_res$auc_gate$first_N_auc_pass)))

## ---- Simulated data (optional export) ----
# The article analyses no empirical data. Section 1 above is computed from
# samples drawn inside the package from the stated parameters under seed
# 2026, so the parameters plus the seed fully determine the data. If a
# tangible file is required, the block below writes those underlying samples
# to CSV; rerunning this script reproduces them byte for byte. (Sections 2
# and 3 have no underlying "samples" to dump: the exact Beta-Binomial
# calculation is closed-form, and joint_sample_size()'s Monte Carlo draws
# are per-candidate-N intermediate counts, not a single fixed-N dataset like
# section 1's.)
if (isTRUE(DUMP_DATA)) {
  sep("Exporting the simulated samples to CSV")
  sim_dir <- file.path(out_dir, "simulated_data")
  dir.create(sim_dir, showWarnings = FALSE, recursive = TRUE)

  # Buderer validation: B binomial samples of successes among 100 diseased
  set.seed(SEED)
  buderer_draws <- rbinom(4000, 100, 0.85)
  write.csv(data.frame(replicate = seq_along(buderer_draws),
                       n_diseased = 100, true_Se = 0.85,
                       observed_positives = buderer_draws),
            file.path(sim_dir, "buderer_validation_samples.csv"),
            row.names = FALSE)

  cat("Written to", sim_dir, ":\n",
      " - buderer_validation_samples.csv  (4,000 replicates)\n", sep = "")
}

## ---- Summary: one line per check, and a CSV in the output directory ----
write.csv(check_log, file.path(out_dir, "reproduce_manuscript_checks.csv"),
          row.names = FALSE)
n_pass <- sum(check_log$result == "PASS")
cat(sprintf("\nSUMMARY: %d of %d checks PASS, %d **CHECK**\n",
            n_pass, nrow(check_log), nrow(check_log) - n_pass))
cat(paste0("  ", check_log$check, ": ", check_log$result), sep = "\n")
cat("\n=== reproduce_manuscript done ===\n")
