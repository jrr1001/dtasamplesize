#!/usr/bin/env Rscript
# =============================================================================
# negative_control_halfwidth.R
#
# Direction of the "negative control" in the Python cross-check
# (conf_PY_05_shared_assumptions.py): which misreading of the target
# convention reproduces the reported assurance of ~0.0032 at N = 678?
#
# Uses dtasamplesize::bam_sample_size(method = "exact") at a SINGLE N
# (N_range = 678L for P1, 672L for P2), so the scan cannot run past that N:
# with N_range = NULL the scan extends to N_max and max_assurance_evaluated
# would refer to a different N.
#
# Scenarios (worked-example priors; target_assurance = 0.80, level 0.95)
#   P1: prior_se = c(17,3), prior_sp = c(2,2),  prior_prev = c(4,16), N = 678
#   P2: prior_se = c(17,3), prior_sp = c(18,2), prior_prev = c(4,16), N = 672
#
# Three readings of the target (delta_* arguments are FULL credible-interval
# widths in bam_sample_size()):
#   (i)   published targets:           (0.14, 0.10)  [full widths]
#   (ii)  delta/2 used AS the full width: (0.07, 0.05)  ("treating the
#         half-width delta/2 as the full-width target": a STRICTER criterion)
#   (iii) 2*delta used as the full width: (0.28, 0.20)  (a LOOSER criterion)
# When the objective is not reached at the single N, joint_assurance is NA
# and the value is read from max_assurance_evaluated, after checking that
# N_at_max_assurance equals the single N.
#
# Usage: Rscript negative_control_halfwidth.R <output_dir>
#   Writes negative_control_halfwidth.csv and negative_control_halfwidth.log
#   to <output_dir> only. If the environment variables DTA_PKG_LIB and
#   DTA_EXTRA_LIB exist they are prepended to .libPaths() (in that order).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript negative_control_halfwidth.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
paths <- paths[nzchar(paths)]
if (length(paths)) .libPaths(c(paths, .libPaths()))

sink(file.path(outdir, "negative_control_halfwidth.log"), split = TRUE)
on.exit(sink(), add = TRUE)

library(dtasamplesize)
TOL_AGREE <- 1e-6   # declared tolerance for agreement with the expected values

cat("negative_control_halfwidth.R\n")
cat("date:          ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:     ", R.version.string, "\n")
cat("dtasamplesize: ", as.character(packageVersion("dtasamplesize")), "\n")
cat("width tolerance: none used in this script (package keeps its exact width <= delta comparison);\n")
cat("agreement with expected values: |difference| <= ", format(TOL_AGREE), "\n", sep = "")
cat("target_assurance = 0.80; credible level 0.95 (alpha_ci default); method = 'exact'\n\n")

scen <- list(
  P1 = list(N = 678L, prior_se = c(17, 3), prior_sp = c(2, 2),  prior_prev = c(4, 16)),
  P2 = list(N = 672L, prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16))
)
readings <- list(
  "(i) published full-width targets"           = c(0.14, 0.10),
  "(ii) delta/2 used as full-width target"     = c(0.07, 0.05),
  "(iii) 2*delta used as full-width target"    = c(0.28, 0.20)
)
expected <- list(
  "P1" = c("(i) published full-width targets" = 0.8003489948,
           "(ii) delta/2 used as full-width target" = 0.003178),
  "P2" = c("(i) published full-width targets" = 0.8002692084,
           "(ii) delta/2 used as full-width target" = 0.038115)
)

rows <- list()
for (sn in names(scen)) for (rn in names(readings)) {
  s <- scen[[sn]]; dl <- readings[[rn]]
  r <- suppressWarnings(bam_sample_size(
    prior_se = s$prior_se, prior_sp = s$prior_sp, prior_prev = s$prior_prev,
    delta_se = dl[1], delta_sp = dl[2], target_assurance = 0.80,
    N_range = s$N, method = "exact"))
  ja   <- r$joint_assurance
  mx   <- r$max_assurance_evaluated
  nmx  <- r$N_at_max_assurance
  used <- if (!is.na(ja)) ja else mx
  exp_v <- if (!is.null(expected[[sn]][rn]) && !is.na(expected[[sn]][rn])) expected[[sn]][[rn]] else NA_real_
  rows[[length(rows) + 1L]] <- data.frame(
    scenario = sn, N = s$N, reading = rn, delta_se_full = dl[1], delta_sp_full = dl[2],
    target_reached = isTRUE(r$target_reached), joint_assurance = ja,
    max_assurance_evaluated = mx, N_at_max_assurance = nmx,
    N_at_max_equals_single_N = identical(as.integer(nmx), s$N),
    assurance_read = used, expected = exp_v,
    abs_diff_vs_expected = if (is.na(exp_v)) NA_real_ else abs(used - exp_v),
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows)
write.csv(tab, file.path(outdir, "negative_control_halfwidth.csv"), row.names = FALSE)

cat("--- Results (assurance_read = joint_assurance if not NA, else max_assurance_evaluated) ---\n")
print(transform(tab, joint_assurance = round(joint_assurance, 10),
                max_assurance_evaluated = round(max_assurance_evaluated, 10),
                assurance_read = round(assurance_read, 10),
                abs_diff_vs_expected = signif(abs_diff_vs_expected, 3)), row.names = FALSE)

cat("\n--- Agreement with expected values ---\n")
for (i in seq_len(nrow(tab))) {
  if (is.na(tab$expected[i])) {
    cat(sprintf("%s %-42s no expected value; read %.10f\n", tab$scenario[i], tab$reading[i], tab$assurance_read[i]))
  } else {
    ok <- tab$abs_diff_vs_expected[i] <= if (grepl("\\(ii\\)", tab$reading[i])) 5e-7 else TOL_AGREE
    cat(sprintf("%s %-42s read %.10f expected %.10f -> %s\n", tab$scenario[i], tab$reading[i],
                tab$assurance_read[i], tab$expected[i], if (ok) "AGREE" else "DIFFER"))
  }
}
cat("(expected (ii) values are given to 6 decimals, so the tolerance for those is 5e-7)\n")
cat("N_at_max_assurance equals the single N in every (ii) row: ",
    all(tab$N_at_max_equals_single_N[grepl("\\(ii\\)", tab$reading)]), "\n", sep = "")

cat("\n--- Direction ---\n")
cat("Only reading (ii) reproduces the negative-control value of about 0.0032 (P1, N = 678):\n")
cat("it corresponds to TREATING THE HALF-WIDTH delta/2 AS THE FULL-WIDTH TARGET, i.e. a\n")
cat("target twice as strict as the published one. Reading (iii), 2*delta as full width, is\n")
cat("the opposite (looser) misreading and gives the assurance shown in the table.\n")
cat("\nDone.\n")
