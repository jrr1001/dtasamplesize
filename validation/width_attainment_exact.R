#!/usr/bin/env Rscript
# =============================================================================
# width_attainment_exact.R
#
# EXACT (binomial enumeration, no Monte Carlo) probability that the realised
# 95% confidence-interval width for a sensitivity (or specificity) estimate
# is <= the target full width 2d, for sample sizes chosen by four sizing rules.
#
# Purpose
#   Replaces the Monte-Carlo sweep (B = 4000) used for the motivating result
#   and Figure 1. Every number here is a finite sum of binomial probabilities.
#
# Usage
#   Rscript width_attainment_exact.R <output_dir>
#   (the output directory is the ONLY thing written to; it is created if
#   missing. Nothing is written anywhere else.)
#   Libraries: if the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB
#   exist they are prepended to .libPaths() (in that order). The package
#   dtasamplesize is only needed for one cross-check; the script runs
#   without it.
#
# Conventions (named explicitly)
#   * target is the FULL width 2d (d is the half-width), two-sided 95% level,
#     z = qnorm(0.975).
#   * All comparisons "width <= target" use the declared tolerance
#     TOL = 1e-12:  width <= 2d + TOL. The same TOL is used for
#     "p_hat(1-p_hat) <= p0(1-p0)" and for "P >= 0.80" in the curve section.
#   * Wald (truncated): p_hat +/- z*sqrt(p_hat(1-p_hat)/n), each limit
#     truncated to [0,1] (as R/mc_validate_buderer.R L93-94 of the package).
#   * Wald (untruncated): same, with no truncation (width = 2 z se).
#   * Wilson score: as in the package (R/mc_validate_buderer.R, "else" branch),
#     limits clamped to [0,1].
#   * Exact coverage: P(lower <= p0 <= upper) with the package's exact
#     comparison (no tolerance), and, as a check, with TOL on both sides.
#
# Sizing rules (per cell Se = p0, d)
#   (a) Buderer:  n = ceiling(z^2 p0 (1-p0) / d^2).
#   (b) Wilson plug-in: smallest n whose Wilson width AT p_hat = p0 is <= 2d.
#   (c) Expected width, Wald (truncated): smallest n such that
#       E[width(n')] <= 2d for EVERY n' in n..N_CAP.
#   (d) Expected width, Wilson: same with the Wilson width.
#   N_CAP (declared) = 5000. E[width] is not monotone in n and, for extreme
#   Se, is small at very small n because degenerate intervals (x in {0,n})
#   have zero width (Wald); the "for every n'" definition avoids that.
#   The first n with E[width] <= 2d is also reported (n_first).
#
# Cell flags (declared thresholds)
#   "borde"     : boundary-dominated: P(x in {0,n}) >= BORDE_P (0.25), or the
#                 expected-width criterion is met only because of the
#                 degenerate points (E[width] <= 2d but the mean width over
#                 1 <= x <= n-1 is > 2d).
#   "discrecion": small-n discreteness: n < DISCRETE_N (30) and not "borde".
#   "ninguna"   : otherwise.
#
# Outputs (in <output_dir>)
#   width_attainment_cells.csv      one row per cell x rule x interval
#   width_attainment_summary.csv    min/Q1/median/Q3/max by interval and rule
#   curve_Se085_d007.csv            P(width <= 0.14) vs n, n = 40..200
#   curve_Se085_d007.png            300 dpi base-graphics figure
#   curve_Se085_d007_summary.csv    first/stable n for P >= 0.80; coverage n=100
#   width_attainment_exact.log      log of the run (also printed on console)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript width_attainment_exact.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
paths <- paths[nzchar(paths)]
if (length(paths)) .libPaths(c(paths, .libPaths()))

logfile <- file.path(outdir, "width_attainment_exact.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)

# ---- declared constants ------------------------------------------------------
TOL        <- 1e-12
Z          <- qnorm(0.975)
N_CAP      <- 5000L
BORDE_P    <- 0.25
DISCRETE_N <- 30L
SE_GRID    <- c(0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.98)
D_GRID     <- c(0.03, 0.05, 0.07, 0.10)

pkg_version <- tryCatch(as.character(utils::packageVersion("dtasamplesize")),
                        error = function(e) "not installed")
cat("width_attainment_exact.R\n")
cat("date:            ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:       ", R.version.string, "\n")
cat("dtasamplesize:   ", pkg_version, "(used only for the Buderer-n cross-check)\n")
cat("tolerance (width <= target, p_hat(1-p_hat) <= p0(1-p0), P >= 0.80): ", format(TOL), "\n")
cat("z = qnorm(0.975) = ", format(Z, digits = 15), "\n")
cat("N_CAP = ", N_CAP, "; BORDE_P = ", BORDE_P, "; DISCRETE_N = ", DISCRETE_N, "\n", sep = "")
cat("output dir:      ", normalizePath(outdir, winslash = "/"), "\n\n")

# ---- interval widths as functions of p_hat and n -----------------------------
# p may be a vector (x/n) or a scalar plug-in value p0; n is a scalar.
w_wald_trunc <- function(p, n) {
  h <- Z * sqrt(p * (1 - p) / n)
  pmin(1, p + h) - pmax(0, p - h)
}
w_wald_untr <- function(p, n) 2 * Z * sqrt(p * (1 - p) / n)
w_wilson <- function(p, n) {
  den <- 1 + Z^2 / n
  ctr <- (p + Z^2 / (2 * n)) / den
  mar <- Z * sqrt((p * (1 - p) + Z^2 / (4 * n)) / n) / den
  pmin(1, ctr + mar) - pmax(0, ctr - mar)
}
W <- list(wald_trunc = w_wald_trunc, wald_untr = w_wald_untr, wilson = w_wilson)

# ---- exact quantities for one (n, p0, d, interval) ---------------------------
cell_stats <- function(n, p0, d, interval) {
  x   <- 0:n
  ph  <- x / n
  pr  <- dbinom(x, n, p0)
  w   <- W[[interval]](ph, n)
  ok  <- w <= 2 * d + TOL
  interior <- x >= 1L & x <= n - 1L
  p_b <- pr[1L] + pr[n + 1L]
  Ew  <- sum(pr * w)
  Ew_int <- if (sum(pr[interior]) > 0) sum(pr[interior] * w[interior]) / sum(pr[interior]) else NA_real_
  list(
    P_width   = sum(pr[ok]),
    P_q       = sum(pr[ph * (1 - ph) <= p0 * (1 - p0) + TOL]),
    P_x0      = pr[1L],
    P_xn      = pr[n + 1L],
    P_x0n     = p_b,
    E_width   = Ew,
    E_width_interior = Ew_int
  )
}
width_at_p0 <- function(n, p0) c(
  wald_trunc = w_wald_trunc(p0, n),
  wald_untr  = w_wald_untr(p0, n),
  wilson     = w_wilson(p0, n)
)
flag_cell <- function(n, P_x0n, E_width, E_int, d, rule_is_E) {
  e_only_degenerate <- rule_is_E && is.finite(E_width) && is.finite(E_int) &&
    E_width <= 2 * d + TOL && E_int > 2 * d + TOL
  if (P_x0n >= BORDE_P || e_only_degenerate) "borde"
  else if (n < DISCRETE_N) "discrecion"
  else "ninguna"
}

# ---- expected-width tables E[width(n)], n = 1..N_CAP, per Se and interval ----
# E[width] does not depend on d, so it is computed once per Se.
cat("Computing E[width(n)] for n = 1..", N_CAP, " for ", length(SE_GRID),
    " Se values x 2 intervals (exact enumeration)...\n", sep = "")
EW <- list(wald_trunc = matrix(NA_real_, length(SE_GRID), N_CAP),
           wilson     = matrix(NA_real_, length(SE_GRID), N_CAP))
for (i in seq_along(SE_GRID)) {
  p0 <- SE_GRID[i]
  for (n in seq_len(N_CAP)) {
    x  <- 0:n
    pr <- dbinom(x, n, p0)
    ph <- x / n
    EW$wald_trunc[i, n] <- sum(pr * w_wald_trunc(ph, n))
    EW$wilson[i, n]     <- sum(pr * w_wilson(ph, n))
  }
}
# n_stable: smallest n with E[width(n')] <= 2d + TOL for all n' in n..N_CAP
n_stable_E <- function(ew, d) {
  bad <- which(ew > 2 * d + TOL)
  if (length(bad) == 0L) return(1L)
  if (max(bad) >= N_CAP) return(NA_integer_)
  max(bad) + 1L
}
n_first_E <- function(ew, d) {
  ok <- which(ew <= 2 * d + TOL)
  if (length(ok) == 0L) NA_integer_ else ok[1L]
}

# ---- Buderer-n cross-check against the package -------------------------------
buderer_formula <- function(p, d) ceiling(Z^2 * p * (1 - p) / d^2)
pkg_ok <- requireNamespace("dtasamplesize", quietly = TRUE)

# ---- main loop over the 28 cells ---------------------------------------------
rows <- list()
add_row <- function(cell, p0, d, rule, interval, n, n_first = NA_integer_) {
  if (is.na(n)) {
    rows[[length(rows) + 1L]] <<- data.frame(
      cell = cell, Se = p0, d = d, rule = rule, interval = interval, n = NA_integer_,
      n_first = n_first, P_width = NA_real_, P_q = NA_real_,
      w_p0_wald_trunc = NA_real_, w_p0_wald_untr = NA_real_, w_p0_wilson = NA_real_,
      P_x0 = NA_real_, P_xn = NA_real_, P_x0n = NA_real_,
      E_width = NA_real_, E_width_interior = NA_real_, flag = NA_character_,
      stringsAsFactors = FALSE)
    return(invisible())
  }
  s  <- cell_stats(n, p0, d, interval)
  wp <- width_at_p0(n, p0)
  fl <- flag_cell(n, s$P_x0n, s$E_width, s$E_width_interior, d, rule %in% c("c_wald_E", "d_wilson_E"))
  rows[[length(rows) + 1L]] <<- data.frame(
    cell = cell, Se = p0, d = d, rule = rule, interval = interval, n = as.integer(n),
    n_first = n_first, P_width = s$P_width, P_q = s$P_q,
    w_p0_wald_trunc = wp[["wald_trunc"]], w_p0_wald_untr = wp[["wald_untr"]],
    w_p0_wilson = wp[["wilson"]],
    P_x0 = s$P_x0, P_xn = s$P_xn, P_x0n = s$P_x0n,
    E_width = s$E_width, E_width_interior = s$E_width_interior, flag = fl,
    stringsAsFactors = FALSE)
}

buderer_mismatch <- 0L
plugin_nonmono <- 0L
n_grid_plug <- seq_len(N_CAP)
for (i in seq_along(SE_GRID)) for (d in D_GRID) {
  p0   <- SE_GRID[i]
  cell <- sprintf("Se%.2f_d%.2f", p0, d)

  # (a) Buderer
  nB <- buderer_formula(p0, d)
  if (pkg_ok && !identical(as.integer(dtasamplesize::buderer_n(p0, d)), as.integer(nB)))
    buderer_mismatch <- buderer_mismatch + 1L
  for (iv in c("wald_trunc", "wald_untr", "wilson")) add_row(cell, p0, d, "a_buderer", iv, nB)

  # (b) Wilson plug-in at p0
  wil_p0 <- vapply(n_grid_plug, function(n) w_wilson(p0, n), numeric(1))
  nW <- which(wil_p0 <= 2 * d + TOL)[1L]
  if (any(wil_p0[nW:N_CAP] > 2 * d + TOL)) plugin_nonmono <- plugin_nonmono + 1L
  for (iv in c("wald_trunc", "wald_untr", "wilson")) add_row(cell, p0, d, "b_wilson_p0", iv, nW)

  # (c), (d) expected width
  nc <- n_stable_E(EW$wald_trunc[i, ], d); ncf <- n_first_E(EW$wald_trunc[i, ], d)
  nd <- n_stable_E(EW$wilson[i, ], d);     ndf <- n_first_E(EW$wilson[i, ], d)
  add_row(cell, p0, d, "c_wald_E",   "wald_trunc", nc, ncf)
  add_row(cell, p0, d, "d_wilson_E", "wilson",     nd, ndf)
}
cells <- do.call(rbind, rows)
write.csv(cells, file.path(outdir, "width_attainment_cells.csv"), row.names = FALSE)

cat("\n--- Cross-checks internal to the script ---\n")
if (pkg_ok) {
  cat("Buderer n: formula vs dtasamplesize::buderer_n() over 28 cells -> mismatches: ",
      buderer_mismatch, "\n", sep = "")
} else cat("Buderer n: package dtasamplesize not available; formula only.\n")
cat("Wilson width at p0 stays <= 2d for all n' >= n_b up to N_CAP in all cells (non-monotone cells count): ",
    plugin_nonmono, "\n", sep = "")
cat("Rows written: ", nrow(cells), " (", length(unique(cells$cell)), " cells)\n", sep = "")

# ---- summary by interval and rule --------------------------------------------
sm <- do.call(rbind, lapply(split(cells, list(cells$rule, cells$interval), drop = TRUE), function(g) {
  g <- g[!is.na(g$P_width), ]
  q <- quantile(g$P_width, c(0, 0.25, 0.5, 0.75, 1), type = 7, names = FALSE)
  data.frame(rule = g$rule[1], interval = g$interval[1], n_cells = nrow(g),
             min = q[1], q1 = q[2], median = q[3], q3 = q[4], max = q[5],
             n_borde = sum(g$flag == "borde"), n_discrecion = sum(g$flag == "discrecion"),
             stringsAsFactors = FALSE)
}))
rownames(sm) <- NULL
sm <- sm[order(sm$rule, sm$interval), ]
write.csv(sm, file.path(outdir, "width_attainment_summary.csv"), row.names = FALSE)
cat("\n--- Summary of P(width <= 2d) by rule and interval (quantile type 7) ---\n")
print(transform(sm, min = round(min, 4), q1 = round(q1, 4), median = round(median, 4),
                q3 = round(q3, 4), max = round(max, 4)), row.names = FALSE)

cat("\n--- Cells (rounded) ---\n")
show <- cells
for (v in c("P_width", "P_q", "P_x0", "P_xn", "P_x0n", "E_width", "E_width_interior"))
  show[[v]] <- round(show[[v]], 4)
for (v in c("w_p0_wald_trunc", "w_p0_wald_untr", "w_p0_wilson")) show[[v]] <- round(show[[v]], 4)
options(width = 250)
print(show, row.names = FALSE)

# ---- curve for Se = 0.85, d = 0.07 -------------------------------------------
cat("\n--- Curve: Se 0.85, d 0.07, target full width 0.14 ---\n")
p0c <- 0.85; dc <- 0.07
ns_curve  <- 40:200
ns_long   <- 40:N_CAP   # longer horizon used to verify the stable n
P_fun <- function(iv, n) cell_stats(n, p0c, dc, iv)$P_width
Pw_long <- vapply(ns_long, function(n) P_fun("wald_trunc", n), numeric(1))
Pl_long <- vapply(ns_long, function(n) P_fun("wilson", n), numeric(1))
curve <- data.frame(n = ns_curve,
                    P_wald_trunc = Pw_long[match(ns_curve, ns_long)],
                    P_wilson     = Pl_long[match(ns_curve, ns_long)])
write.csv(curve, file.path(outdir, "curve_Se085_d007.csv"), row.names = FALSE)

first_ge <- function(P, ns, thr = 0.80) { i <- which(P >= thr - TOL)[1L]; if (is.na(i)) NA_integer_ else ns[i] }
stable_ge <- function(P, ns, thr = 0.80) {
  bad <- which(P < thr - TOL)
  if (length(bad) == 0L) return(ns[1L])
  if (max(bad) == length(ns)) return(NA_integer_)
  ns[max(bad) + 1L]
}
sel200 <- ns_long <= 200L
res_curve <- data.frame(
  interval = c("wald_trunc", "wilson"),
  first_n_P_ge_0.80 = c(first_ge(Pw_long, ns_long), first_ge(Pl_long, ns_long)),
  stable_n_to_200   = c(stable_ge(Pw_long[sel200], ns_long[sel200]), stable_ge(Pl_long[sel200], ns_long[sel200])),
  stable_n_to_N_CAP = c(stable_ge(Pw_long, ns_long), stable_ge(Pl_long, ns_long)),
  n_nonmonotone_steps_40_200 = c(sum(diff(Pw_long[sel200]) < -TOL), sum(diff(Pl_long[sel200]) < -TOL)),
  stringsAsFactors = FALSE)

# exact coverage at n = 100
cov_fun <- function(iv, n, p0, tol = 0) {
  x <- 0:n; ph <- x / n
  if (iv == "wald_trunc") {
    h <- Z * sqrt(ph * (1 - ph) / n); lo <- pmax(0, ph - h); hi <- pmin(1, ph + h)
  } else {
    den <- 1 + Z^2 / n; ctr <- (ph + Z^2 / (2 * n)) / den
    mar <- Z * sqrt((ph * (1 - ph) + Z^2 / (4 * n)) / n) / den
    lo <- pmax(0, ctr - mar); hi <- pmin(1, ctr + mar)
  }
  sum(dbinom(x, n, p0)[(lo <= p0 + tol) & (p0 <= hi + tol)])
}
res_curve$coverage_n100_exact_comparison <- c(cov_fun("wald_trunc", 100, p0c), cov_fun("wilson", 100, p0c))
res_curve$coverage_n100_with_TOL         <- c(cov_fun("wald_trunc", 100, p0c, TOL), cov_fun("wilson", 100, p0c, TOL))
res_curve$P_width_n100 <- c(Pw_long[ns_long == 100L], Pl_long[ns_long == 100L])
write.csv(res_curve, file.path(outdir, "curve_Se085_d007_summary.csv"), row.names = FALSE)
print(res_curve, row.names = FALSE, digits = 7)

png(file.path(outdir, "curve_Se085_d007.png"), width = 6.5, height = 4.5, units = "in", res = 300)
par(mar = c(4.2, 4.4, 1.2, 0.8))
plot(curve$n, curve$P_wald_trunc, type = "l", lwd = 2, col = "black", ylim = c(0, 1),
     xlab = "n (number of diseased subjects)",
     ylab = "Exact P(CI full width <= 0.14)", las = 1)
lines(curve$n, curve$P_wilson, lwd = 2, col = "grey45", lty = 2)
abline(h = 0.80, lty = 3)
abline(v = 100, lty = 3, col = "grey60")
legend("bottomright", c("Wald (truncated to [0,1])", "Wilson score", "P = 0.80", "n = 100 (Buderer)"),
       lwd = c(2, 2, 1, 1), lty = c(1, 2, 3, 3), col = c("black", "grey45", "black", "grey60"), bty = "n")
dev.off()

# ---- agreement with the expected (preliminary) values ------------------------
cat("\n--- Agreement with expected values from planning evidence ---\n")
chk <- function(label, got, exp, digits) {
  ok <- isTRUE(all.equal(round(got, digits), round(exp, digits)))
  cat(sprintf("%-62s got %-12s expected %-12s -> %s\n", label,
              paste(format(round(got, digits), nsmall = digits), collapse = "/"),
              paste(format(round(exp, digits), nsmall = digits), collapse = "/"),
              if (ok) "AGREE" else "DIFFER"))
}
g <- function(cell, rule, iv) cells[cells$cell == cell & cells$rule == rule & cells$interval == iv, ]
chk("Se.85 d.07 Buderer n", g("Se0.85_d0.07", "a_buderer", "wald_trunc")$n, 100, 0)
chk("Se.85 d.07 n=100 P(width<=.14) Wald trunc", g("Se0.85_d0.07", "a_buderer", "wald_trunc")$P_width, 0.5683, 4)
chk("Se.85 d.07 n=100 P(width<=.14) Wilson", g("Se0.85_d0.07", "a_buderer", "wilson")$P_width, 0.5683, 4)
wb <- cells[cells$rule == "b_wilson_p0" & cells$interval == "wilson", ]
chk("Wilson own plug-in rule: range over 28 cells (min, max)", range(wb$P_width), c(0.482, 0.721), 3)
r98 <- g("Se0.98_d0.05", "a_buderer", "wald_trunc")
chk("Se.98 d.05 Buderer n", r98$n, 31, 0)
chk("Se.98 d.05 Wald trunc P at n=31", r98$P_width, 0.8728, 4)
chk("Se.98 d.05 Wald trunc width at p0, n=31", r98$w_p0_wald_trunc, 0.0693, 4)
chk("First n with P>=0.80, Wald trunc / Wilson", res_curve$first_n_P_ge_0.80, c(114, 114), 0)
chk("Stable n (to 200) with P>=0.80, Wald trunc / Wilson", res_curve$stable_n_to_200, c(117, 114), 0)
chk("Exact coverage n=100: Wald trunc / Wilson", res_curve$coverage_n100_exact_comparison, c(0.9330, 0.9332), 4)
cat("Stable n with horizon N_CAP (should equal horizon-200 values): ",
    paste(res_curve$stable_n_to_N_CAP, collapse = " / "), "\n")
cat("\nDone.\n")
