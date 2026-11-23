#!/usr/bin/env Rscript
# =============================================================================
# edge_case_empty_arm.R
#
# Empty-arm edge case of the exact joint assurance: what happens when the
# number of diseased subjects K (Se arm size) is 0, or N (so the Sp arm is
# empty), and how much the convention used for that case matters.
#
# Implementation: SELF-CONTAINED exact beta-binomial enumeration (base R only;
# no Monte Carlo), cross-checked against dtasamplesize::bam_sample_size(
# method = "exact") at the same single N when the package is installed.
#
# Model (as in bam_sample_size): K ~ BetaBinomial(N, a_prev, b_prev);
# given K, the Se arm has size K and the Sp arm has size N - K. For an arm of
# size m with a Beta(a, b) prior and x successes the equal-tailed credible
# interval is Beta(a + x, b + m - x); P_arm(m) = P(width <= delta) with
# x ~ BetaBinomial(m, a, b). Joint assurance
#   A(N) = sum_K P(K) * P_se(K) * P_sp(N - K).
#
# Operational definition of "suppressing the convention"
#   Convention in the package ("empty arm = failure"): P_arm(0) = 0 for BOTH
#   arms, i.e. any replication with K = 0 or K = N contributes 0 to A(N).
#   Suppressing it means: remove that override and evaluate the empty arm
#   (m = 0, x = 0) with the SAME formula as every other arm, i.e. by its
#   prior-only credible interval Beta(a, b):
#       P_arm(0) = 1 if width(prior) <= delta + TOL, else 0.
#   All terms with m >= 1 are unchanged. Hence the assurance can change only
#   if the prior-only width of the empty arm is <= its target delta, and the
#   change is bounded by P(K in {0, N}).
#
# Declared tolerances: width <= delta uses TOL_W = 1e-12 (the package keeps its
# exact comparison; no value here lies within 1e-12 of its target, which the
# script checks and prints). Agreement between this implementation and the
# package: |difference| <= TOL_IMPL = 1e-10 on A(N).
#
# Scenarios
#   C10 (LOCKED_confirmatory_grid.json): Se Beta(8,2), Sp Beta(2,8), prevalence
#       Beta(3,7), full-width targets 0.52 (Se) and 0.48 (Sp), level 0.90, N = 3.
#   Worked example: Se Beta(17,3), Sp Beta(2,2), prevalence Beta(4,16), full-width
#       targets 0.14 (Se) and 0.10 (Sp), level 0.95, at N = 3 and N = 678.
#
# Usage: Rscript edge_case_empty_arm.R <output_dir>
#   Writes edge_case_empty_arm.csv and edge_case_empty_arm.log to <output_dir>
#   only. If the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB exist
#   they are prepended to .libPaths() (in that order).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript edge_case_empty_arm.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
paths <- paths[nzchar(paths)]
if (length(paths)) .libPaths(c(paths, .libPaths()))

sink(file.path(outdir, "edge_case_empty_arm.log"), split = TRUE)
on.exit(sink(), add = TRUE)

TOL_W    <- 1e-12
TOL_IMPL <- 1e-10
pkg_version <- tryCatch(as.character(utils::packageVersion("dtasamplesize")), error = function(e) "not installed")
pkg_ok <- requireNamespace("dtasamplesize", quietly = TRUE)

cat("edge_case_empty_arm.R\n")
cat("date:          ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:     ", R.version.string, "\n")
cat("dtasamplesize: ", pkg_version, "(cross-check only; the computation is self-contained beta-binomial enumeration)\n")
cat("tolerance: width <= delta uses TOL_W = ", format(TOL_W), "; implementation agreement |dA| <= ", format(TOL_IMPL), "\n", sep = "")
cat("\nOperational definition of 'suppressing the convention':\n")
cat("  package convention = 'empty arm = failure': P_arm(0) = 0 for both arms, so any replication\n")
cat("  with K = 0 or K = N contributes 0 to the joint assurance.\n")
cat("  suppressed convention = remove that override and evaluate the empty arm (m = 0, x = 0)\n")
cat("  by its prior-only credible interval Beta(a, b) with the same formula as every other arm:\n")
cat("  P_arm(0) = 1 if prior width <= delta (+TOL_W), else 0. Terms with m >= 1 are unchanged.\n")
cat("  The assurance can therefore change only if the prior-only width of an empty arm is <= its\n")
cat("  target, and by at most P(K in {0, N}).\n\n")

# ---- helpers ------------------------------------------------------------------
bb_pmf <- function(k, n, a, b) exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))
cred_width <- function(a, b, level) qbeta(1 - (1 - level) / 2, a, b) - qbeta((1 - level) / 2, a, b)

# P_arm(m) for m = 0..M; empty_ok = FALSE -> P_arm(0) = 0 (package convention)
p_arm <- function(M, a, b, delta, level, empty_ok) {
  P <- numeric(M + 1L)
  P[1L] <- if (empty_ok) as.numeric(cred_width(a, b, level) <= delta + TOL_W) else 0
  for (m in seq_len(M)) {
    x <- 0:m
    w <- cred_width(a + x, b + m - x, level)
    P[m + 1L] <- sum(bb_pmf(x, m, a, b)[w <= delta + TOL_W])
  }
  P
}
assurance <- function(N, Pse, Psp, ap, bp) {
  k <- 0:N
  sum(bb_pmf(k, N, ap, bp) * Pse[k + 1L] * Psp[N - k + 1L])
}
min_margin <- function(M, a, b, delta, level) {   # smallest |width - delta| over all (m, x)
  mm <- Inf
  for (m in seq_len(M)) { x <- 0:m; mm <- min(mm, abs(cred_width(a + x, b + m - x, level) - delta)) }
  mm
}

run_case <- function(label, prior_se, prior_sp, prior_prev, dse, dsp, level, N, target) {
  ap <- prior_prev[1]; bp <- prior_prev[2]
  p_deg <- bb_pmf(0, N, ap, bp) + bb_pmf(N, N, ap, bp)
  w_se0 <- cred_width(prior_se[1], prior_se[2], level)
  w_sp0 <- cred_width(prior_sp[1], prior_sp[2], level)
  out <- list()
  for (eo in c(FALSE, TRUE)) {
    Pse <- p_arm(N, prior_se[1], prior_se[2], dse, level, eo)
    Psp <- p_arm(N, prior_sp[1], prior_sp[2], dsp, level, eo)
    out[[as.character(eo)]] <- assurance(N, Pse, Psp, ap, bp)
  }
  margin <- min(min_margin(min(N, 60L), prior_se[1], prior_se[2], dse, level),
                min_margin(min(N, 60L), prior_sp[1], prior_sp[2], dsp, level))
  # package cross-check at the same single N (convention 'empty arm = failure')
  pkg_A <- NA_real_
  if (pkg_ok && N <= 700L) {
    r <- suppressWarnings(dtasamplesize::bam_sample_size(
      prior_se = prior_se, prior_sp = prior_sp, prior_prev = prior_prev,
      delta_se = dse, delta_sp = dsp, target_assurance = target, alpha_ci = level,
      N_range = as.integer(N), method = "exact"))
    pkg_A <- if (!is.na(r$joint_assurance)) r$joint_assurance else r$max_assurance_evaluated
  }
  data.frame(
    case = label, N = N, level = level, delta_se = dse, delta_sp = dsp,
    P_K_in_0_N = p_deg, prior_width_se = w_se0, prior_width_sp = w_sp0,
    A_empty_fail = out[["FALSE"]], A_empty_prior_interval = out[["TRUE"]],
    shift = out[["TRUE"]] - out[["FALSE"]],
    min_abs_width_minus_delta = margin,
    A_package_exact = pkg_A,
    abs_diff_package_vs_script = abs(pkg_A - out[["FALSE"]]),
    stringsAsFactors = FALSE)
}

res <- rbind(
  run_case("C10_N3",            c(8, 2),  c(2, 8), c(3, 7),  0.52, 0.48, 0.90, 3L,   0.50),
  run_case("worked_N3",         c(17, 3), c(2, 2), c(4, 16), 0.14, 0.10, 0.95, 3L,   0.80),
  run_case("worked_N678",       c(17, 3), c(2, 2), c(4, 16), 0.14, 0.10, 0.95, 678L, 0.80)
)
write.csv(res, file.path(outdir, "edge_case_empty_arm.csv"), row.names = FALSE)

old <- options(width = 250, digits = 10); on.exit(options(old), add = TRUE)
cat("--- Results ---\n")
print(res, row.names = FALSE)
cat("\nWorked example (P1) prior widths at 95%: Se Beta(17,3) = ",
    format(cred_width(17, 3, 0.95), digits = 6), "; Sp Beta(2,2) = ", format(cred_width(2, 2, 0.95), digits = 6),
    " (both above their targets 0.14 and 0.10, so the empty-arm convention cannot change A).\n", sep = "")
cat("C10 prior widths at 90%: Se Beta(8,2) = ", format(cred_width(8, 2, 0.90), digits = 6),
    "; Sp Beta(2,8) = ", format(cred_width(2, 8, 0.90), digits = 6),
    " (both below 0.52 and 0.48: the prior-only interval already meets the target).\n", sep = "")

# ---- agreement with expected values ------------------------------------------
cat("\n--- Agreement with expected values from planning evidence ---\n")
chk <- function(label, got, exp, digits) {
  ok <- isTRUE(all.equal(signif(got, digits), signif(exp, digits)))
  cat(sprintf("%-58s got %-14s expected %-14s -> %s\n", label, format(signif(got, digits)),
              format(signif(exp, digits)), if (ok) "AGREE" else "DIFFER"))
}
g <- function(case) res[res$case == case, ]
chk("C10 N=3  P(K in {0,3})",                     g("C10_N3")$P_K_in_0_N, 0.4273, 4)
chk("C10 N=3  A, empty arm = failure",            g("C10_N3")$A_empty_fail, 0.5727, 4)
chk("C10 N=3  A, empty arm by prior-only interval", g("C10_N3")$A_empty_prior_interval, 1.0000, 4)
chk("C10 prior width Se",                         g("C10_N3")$prior_width_se, 0.388, 3)
chk("C10 prior width Sp",                         g("C10_N3")$prior_width_sp, 0.388, 3)
chk("Worked N=3  P(K in {0,3})",                  g("worked_N3")$P_K_in_0_N, 0.5429, 4)
chk("Worked N=3  A, empty arm = failure (unchanged by convention)", g("worked_N3")$A_empty_fail, g("worked_N3")$A_empty_prior_interval, 6)
chk("Worked prior width Se",                      g("worked_N3")$prior_width_se, 0.298, 3)
chk("Worked prior width Sp",                      g("worked_N3")$prior_width_sp, 0.811, 3)
chk("Worked N=678  P(K in {0,678})",              g("worked_N678")$P_K_in_0_N, 3.976e-07, 4)
chk("Worked N=678  A(678), empty arm = failure (published P1 0.8003489948)",
    g("worked_N678")$A_empty_fail, 0.8003489948, 10)
cat("shift in A due to the convention: worked N=3 = ", format(g("worked_N3")$shift), "; worked N=678 = ",
    format(g("worked_N678")$shift), "\n", sep = "")
if (pkg_ok) {
  d <- res$abs_diff_package_vs_script
  cat("Package (bam_sample_size, method='exact', N_range = single N) vs this script, max |dA| = ",
      format(max(d, na.rm = TRUE)), " -> ", if (all(d <= TOL_IMPL, na.rm = TRUE)) "AGREE" else "DIFFER",
      " (tolerance ", format(TOL_IMPL), ")\n", sep = "")
}
cat("Smallest |width - delta| over arm sizes 1..min(N,60) and all x (should be >> TOL_W): ",
    format(min(res$min_abs_width_minus_delta)), "\n", sep = "")
cat("\nDone.\n")
