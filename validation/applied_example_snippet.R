#!/usr/bin/env Rscript
# =============================================================================
# applied_example_snippet.R
#
# Short batch fragment for the manuscript. The code between the
# BEGIN and END markers (25 lines or fewer) is the fragment: it runs as is and
# reproduces the main table of validation/applied_example.R (published Panbio
# pilot, Albert et al. 2021, doi:10.1016/j.cmi.2020.11.004: TP 43, FN 11,
# FP 0, TN 358; see that script for the data sources and conventions).
# Everything outside the markers is the verification harness, not part of the
# fragment.
#
# Usage
#   Rscript applied_example_snippet.R <output_dir> <applied_example_main_table.csv>
#   Writes ONLY to <output_dir>: log and the reproduced table (CSV).
#   Libraries: DTA_PKG_LIB and DTA_EXTRA_LIB, if set, are prepended to
#   .libPaths(). Single core; exact calculation (no Monte Carlo).
#
# Conventions: one-sided 95% lower semi-widths 0.10 (Se) and 0.05 (Sp), or
# two-sided 95% full widths 0.20 and 0.10; joint objective unless "se";
# assurance target 0.80; total N; base prior Beta(1,1); a0 = 0.5; analysis
# prior = base ("design") or the updated prior ("both"). Comparison with the
# reference CSV: N identical, A identical to 1e-12 (TOL). The informative package
# warning for use_in = "both" (analysis prior already within the target width)
# is not suppressed: it is printed once, as a single message, when the fragment runs.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript applied_example_snippet.R <output_dir> <main_table.csv>")
outdir <- args[1L]; ref_csv <- args[2L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
libs <- c(Sys.getenv("DTA_PKG_LIB"), Sys.getenv("DTA_EXTRA_LIB")); libs <- libs[nzchar(libs)]
if (length(libs)) .libPaths(c(libs, .libPaths()))
logcon <- file(file.path(outdir, "applied_example_snippet.log"), open = "wt", encoding = "UTF-8")
sink(logcon, split = TRUE); on.exit({ sink(); close(logcon) }, add = TRUE)
t_start <- Sys.time()

# --- BEGIN MANUSCRIPT FRAGMENT -----------------------------------------------
library(dtasamplesize)
pilot <- list(tp = 43, fn = 11, tn = 358, fp = 0)            # Albert et al. 2021
pri <- function(use_in) do.call(bam_prior_from_pilot, c(pilot, list(
  base_prior_se = c(1, 1), base_prior_sp = c(1, 1), base_prior_prev = c(1, 1),
  a0 = 0.5, use_in = use_in)))
p_des <- pri("design"); p_both <- pri("both")
two <- list(width_type = "full_width", delta_se = 0.20, delta_sp = 0.10)
scen <- list(
  design_onesided_joint = list(priors = p_des),
  both_onesided_joint   = list(priors = p_both),
  design_twosided_joint = c(list(priors = p_des), two),
  both_twosided_joint   = c(list(priors = p_both), two),
  design_onesided_se    = list(priors = p_des, assurance_for = "se"),
  design_twosided_se    = c(list(priors = p_des, assurance_for = "se"), two[1:2]))
tab <- bam_sensitivity(scen, delta_se = 0.10, delta_sp = 0.05,
  width_type = "lower_half_width", alpha_ci = 0.95, target_assurance = 0.80)
tab[, c("scenario", "N_total", "A_N", "A_N_minus_1")]
# --- END MANUSCRIPT FRAGMENT -------------------------------------------------

cat("\n== harness ==\n")
cat("R:", R.version.string, "| dtasamplesize", as.character(packageVersion("dtasamplesize")),
    "|", find.package("dtasamplesize"), "\n")
print(tab[, c("scenario", "N_total", "A_N", "A_N_minus_1")], digits = 10)
write.csv(tab, file.path(outdir, "applied_example_snippet_table.csv"), row.names = FALSE)
TOL <- 1e-12
ref <- read.csv(ref_csv, stringsAsFactors = FALSE)
stopifnot(identical(sort(ref$scenario), sort(tab$scenario)))
ref <- ref[match(tab$scenario, ref$scenario), ]
dN <- tab$N_total - ref$N_total
dA <- abs(tab$A_N - ref$A_N); dA1 <- abs(tab$A_N_minus_1 - ref$A_N_minus_1)
cmp <- data.frame(scenario = tab$scenario, N_snippet = tab$N_total, N_reference = ref$N_total,
                  diff_N = dN, A_snippet = tab$A_N, A_reference = ref$A_N, absdiff_A = dA,
                  absdiff_A_Nminus1 = dA1)
print(cmp, digits = 12, row.names = FALSE)
ok <- all(dN == 0L) && all(dA <= TOL) && all(dA1 <= TOL)
cat(sprintf("\nCOMPARISON WITH %s: N identical = %s; max|dA| = %.3e; max|dA(N-1)| = %.3e; TOL = %g -> %s\n",
            basename(ref_csv), all(dN == 0L), max(dA), max(dA1), TOL, if (ok) "MATCH" else "MISMATCH"))
cat(sprintf("elapsed: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
if (!ok) quit(status = 1)
