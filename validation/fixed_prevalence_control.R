#!/usr/bin/env Rscript
# =============================================================================
# fixed_prevalence_control.R
#
# Control for the decomposition of the 678 -> 533 reduction of the negative
# control in the article (where the number of diseased participants is fixed
# at round(0.20 N) instead of being Beta-Binomial). That control fixes two
# things at once: the prevalence (at its prior mean 0.20) and the number of
# diseased participants (at its expectation). This script separates them: it
# keeps the number of diseased participants binomial and fixes only the
# prevalence at 0.20.
#
# bam_sample_size() has no argument for a fixed prevalence, so the script
# approximates it with a very concentrated prior, prevalence ~ Beta(2e6, 8e6)
# (mean 0.20, standard deviation 1.3e-4). With such a prior the
# Beta-Binomial distribution of the number of diseased participants is, to
# practical precision, Binomial(N, 0.20) with a known prevalence, and the
# exact method (a finite sum) needs no change. A prior this concentrated is a
# computational device, not a model of the prevalence: it stands for "the
# prevalence is known to be 0.20".
#
# Worked example of the article (full-width targets 0.14 for Se and 0.10 for
# Sp, 95% equal-tailed intervals, target assurance 0.80, exact method):
#   P1: Se ~ Beta(17, 3), Sp ~ Beta(2, 2)
#   P2: Se ~ Beta(17, 3), Sp ~ Beta(18, 2)
# With the Beta(4, 16) prevalence prior the exact sample sizes are 678 (P1)
# and 672 (P2). With the number of diseased participants fixed at
# round(0.20 N) they are 533 (P1 and P2; variant V5 of
# validation/python/conf_PY_05_shared_assumptions.py, which is a registered
# file and is not edited here). Expected here, with the prevalence fixed and
# the number of diseased participants binomial: 540 for both priors.
#
# Usage
#   Rscript validation/fixed_prevalence_control.R <output_dir>
#   Writes ONLY into <output_dir>: fixed_prevalence_control.log and
#   fixed_prevalence_control.csv.
#   Libraries: if the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB
#   exist they are prepended to .libPaths() (in that order). Single core; each
#   design takes about one minute.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript fixed_prevalence_control.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
paths <- paths[nzchar(paths)]
if (length(paths)) .libPaths(c(paths, .libPaths()))
suppressPackageStartupMessages(library(dtasamplesize))

sink(file.path(outdir, "fixed_prevalence_control.log"), split = TRUE)
on.exit(sink(), add = TRUE)
cat("fixed_prevalence_control.R\n")
cat("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n", sep = "")
cat("R: ", R.version.string, "\n", sep = "")
cat("dtasamplesize ", as.character(packageVersion("dtasamplesize")),
    " (", find.package("dtasamplesize"), ")\n", sep = "")
cat("Prevalence prior used to fix the prevalence at 0.20: Beta(2e6, 8e6)\n\n")

prev_fixed <- c(2e6, 8e6)
run <- function(design, sp) {
  t0 <- Sys.time()
  r <- bam_sample_size(prior_se = c(17, 3), prior_sp = sp, prior_prev = prev_fixed, method = "exact")
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("%s: N_total = %s, joint assurance = %.10f, %.1f s\n", design, r$N_total, r$joint_assurance, secs))
  data.frame(design = design, prior_se = "Beta(17,3)", prior_sp = sprintf("Beta(%g,%g)", sp[1], sp[2]),
             prior_prev = "Beta(2e6,8e6)", N_total = r$N_total, joint_assurance = r$joint_assurance,
             seconds = round(secs, 1), stringsAsFactors = FALSE)
}
tab <- rbind(run("P1", c(2, 2)), run("P2", c(18, 2)))
write.csv(tab, file.path(outdir, "fixed_prevalence_control.csv"), row.names = FALSE)

cat("\nReading: with the prevalence fixed at 0.20 and the number of diseased participants\n")
cat("binomial, the sample size is the value above for each design. Compare with 678 and 672\n")
cat("(Beta(4,16) prevalence prior) and with 533 (number of diseased participants fixed at\n")
cat("round(0.20 N)).\n")
ok <- all(tab$N_total == 540L)
cat(sprintf("\nExpected 540 for both designs: %s\n", if (ok) "PASS" else "**CHECK**"))
