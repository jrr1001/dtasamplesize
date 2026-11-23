#!/usr/bin/env Rscript
# =============================================================================
# mc_validate_buderer_exact_vs_reference.R
#
# Compares mc_validate_buderer(method = "exact") with the separate (base R,
# same team)
# reference curve produced by validation/width_attainment_exact.R
# (file curve_Se085_d007.csv: P(width <= 0.14) vs n, n = 40..200, Se = 0.85,
# d = 0.07, Wald truncated to [0, 1] and Wilson).
#
# Usage:
#   Rscript mc_validate_buderer_exact_vs_reference.R <curve_Se085_d007.csv>
#
# The reference CSV is read, never written. Output goes to the console only.
# Tolerance: 1e-12 (absolute) on every P(width <= 0.14) value.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript <script> <curve_Se085_d007.csv>")
ref <- utils::read.csv(args[1L])
stopifnot(identical(names(ref), c("n", "P_wald_trunc", "P_wilson")))

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(dtasamplesize)
}

TOL <- 1e-12
got <- t(vapply(ref$n, function(n) {
  r <- mc_validate_buderer(Se = 0.85, d = 0.07, n_diseased = n, method = "exact")$results
  c(r$P_width_target[r$ci_method == "wald"], r$P_width_target[r$ci_method == "wilson"])
}, numeric(2)))
dw <- max(abs(got[, 1] - ref$P_wald_trunc))
dl <- max(abs(got[, 2] - ref$P_wilson))
cat("rows compared:", nrow(ref), "\n")
cat("max |difference| Wald truncated:", format(dw), "\n")
cat("max |difference| Wilson:        ", format(dl), "\n")
cat("all within", format(TOL), ":", dw <= TOL && dl <= TOL, "\n")
if (!(dw <= TOL && dl <= TOL)) quit(status = 1)
