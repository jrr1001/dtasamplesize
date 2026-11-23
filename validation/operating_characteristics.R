#!/usr/bin/env Rscript
# =============================================================================
# operating_characteristics.R
#
# EXACT operating characteristics (finite sums, no Monte Carlo) of the five
# design/analysis combinations of the worked example of the manuscript
# (v0.6.9), each evaluated with its OWN analysis method, at FIXED true values
# of Se, Sp and prevalence:
#   (1) Buderer N = 500, Wald interval         (also reported with Wilson, (2))
#   (3) Joint frequentist N = 580, Wilson interval   [joint_sample_size()]
#   (4) BAM exact N = 678 (P1: Se ~ Beta(17,3), Sp ~ Beta(2,2), vague Sp)
#   (5) BAM exact N = 672 (P2: Se ~ Beta(17,3), Sp ~ Beta(18,2), informative Sp)
# For every design and every true (Se, Sp, prev) on the declared grid it gives
#   * P(both intervals meet their target width)  ("success"),
#   * the frequentist coverage of the Se interval and of the Sp interval,
#   * the bias of the point estimator (MLE in the frequentist designs,
#     posterior mean in the Bayesian ones).
#
# Usage
#   Rscript operating_characteristics.R <output_dir>
#   Writes ONLY into <output_dir> (created if missing): log, CSV tables, PNG.
#   If the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB exist they are
#   prepended to .libPaths() (in that order). Single core; no Monte Carlo.
#
# INPUTS OF THE DESIGNS (all from the 0.6.9 article)
#   Planning values Se = 0.85, Sp = 0.90, prevalence = 0.20, 95% level
#       (Table 2 of the article; Se = 0.85, d = 0.07).
#   Targets: HALF-widths 0.07 (Se) and 0.05 (Sp) for the classical and Monte
#       Carlo methods; equivalently FULL widths 0.14 and 0.10 for the Bayesian
#       method (delta_se = 0.14, delta_sp = 0.10).
#       In this script every success event is "full width <= 0.14 (Se) AND
#       full width <= 0.10 (Sp)", the same event for all five designs.
#   Buderer N = 500: Table 2 of the article; = buderer_total_N(0.85, 0.90,
#       0.07, 0.05, 0.20) (checked in the log).
#   Joint frequentist N = 580 (first candidate of seq(100,800,10); Wilson,
#       cohort design, K ~ Binomial(N, prev)): Table 2 of the article (579 on
#       the Monte Carlo integer grid; 578 by exact enumeration).
#   BAM N = 678 under P1 (target assurance 0.80): priors Se ~ Beta(17,3),
#       Sp ~ Beta(2,2), prev ~ Beta(4,16), delta_se = 0.14, delta_sp = 0.10
#       (also in validation/compat_v069_check.R, calls P1 and P2).
#   BAM N = 672 under P2: Sp ~ Beta(18,2), the rest as P1;
#       validation/compat_v069_check.R (call P2).
#   Published values compared in the log: A(678) = 0.8003489948, A(677) =
#       0.7996848824, A(672) = 0.8002692084, joint Wilson MC 0.8041 at N = 580
#       with MCSE 0.002806 (Table 2 of the article), exact enumeration
#       A(578) = 0.801244683 and A(577) = 0.7983424783
#       (validation/comparator_calls.R), exact Wald/Wilson P(width <= 0.14) at
#       n = 100, Se = 0.85: 0.5683 (validation/width_attainment_exact.R).
#
# CONVENTIONS (named)
#   * Level 0.95, two-sided, equal tails, z = qnorm(0.975) in Wald and Wilson.
#     Targets are FULL widths; a width meets the target when
#     width <= target + TOL, TOL = 1e-12 (rule 11).
#   * Wald: p +/- z sqrt(p(1-p)/n), limits truncated to [0,1] (as
#     wald_ci() of the package). A degenerate Wald interval (x = 0 or x = n) has
#     width 0 and therefore meets the target (as in mc_validate_buderer()); it
#     does not cover a true value in (0,1). Wilson: as wilson_ci() of the
#     package (limits clamped to [0,1]). Bayesian: equal-tailed interval
#     [qbeta(0.025), qbeta(0.975)] of the posterior Beta(a + x, b + n - x) with
#     the ANALYSIS prior of the design (equal to the design prior in P1 and P2).
#   * Number of diseased K ~ Binomial(N, prev) (cohort sampling); Se arm: x ~
#     Binomial(K, Se); Sp arm: y ~ Binomial(N - K, Sp), independent given K.
#   * Empty arm (K = 0 or N - K = 0): counts as FAILURE of the design for
#     "success" (as the package: P_arm(0) = 0). For coverage and bias an
#     empty arm has no estimate: coverage/bias of the Se interval are computed
#     CONDITIONAL on K >= 1 and those of the Sp interval conditional on
#     N - K >= 1. P(K = 0) and P(N - K = 0) are reported at every grid point
#     (columns/magnitudes p_empty_se_arm, p_empty_sp_arm).
#   * Coverage = P(lower <= truth <= upper), closed interval, no tolerance.
#   * Bias = E[estimate] - truth (conditional on a non-empty arm).
#   * Grid of true values (declared): Se = 0.65, 0.675, ..., 0.975 (14 values;
#     planning value 0.85 and 0.025 steps on both sides, upper end 0.975);
#     Sp = 0.70, 0.725, ..., 0.975 (12 values; planning 0.90); prevalence =
#     0.10, 0.20, 0.30 (planning 0.20). Full factorial, 504 points per design.
#
# OUTPUTS (in <output_dir>)
#   operating_characteristics.log     log of the run
#   oc_cells_long.csv                 (design, Se, Sp, prevalence, quantity, value)
#   oc_cells_wide.csv                 same, one row per (design, Se, Sp, prev)
#   oc_summary_planning_point.csv     the five designs at (0.85, 0.90, 0.20)
#   oc_summary_by_design.csv          range of each magnitude over the grid
#   oc_prior_average_check.csv        average over the design prior vs A(N)
#   oc_success_vs_Se.png              success probability vs true Se
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript operating_characteristics.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pkg_lib   <- Sys.getenv("DTA_PKG_LIB",   unset = "")
extra_lib <- Sys.getenv("DTA_EXTRA_LIB", unset = "")
newlibs <- c(if (nzchar(pkg_lib)) pkg_lib, if (nzchar(extra_lib)) extra_lib)
if (length(newlibs)) .libPaths(c(newlibs, .libPaths()))

logfile <- file.path(outdir, "operating_characteristics.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)
t_start <- Sys.time()

suppressMessages(library(dtasamplesize))
cat("operating_characteristics.R\n")
cat("date:                 ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R:                    ", R.version.string, "\n")
cat("dtasamplesize:        ", as.character(packageVersion("dtasamplesize")), "\n")
cat("find.package:         ", normalizePath(find.package("dtasamplesize"), winslash = "/"), "\n")
cat("output dir:           ", normalizePath(outdir, winslash = "/"), "\n")

# ---- declared constants ------------------------------------------------------
TOL    <- 1e-12
ALPHA  <- 0.05
Z      <- qnorm(1 - ALPHA / 2)
DELTA_SE <- 0.14       # full width target, Se
DELTA_SP <- 0.10       # full width target, Sp
PLAN <- c(Se = 0.85, Sp = 0.90, prev = 0.20)
SE_GRID   <- round(seq(0.65, 0.975, by = 0.025), 3)
SP_GRID   <- round(seq(0.70, 0.975, by = 0.025), 3)
PREV_GRID <- c(0.10, 0.20, 0.30)
cat("TOL (width <= target + TOL): ", format(TOL), "; level 0.95; z = ", format(Z, digits = 15), "\n", sep = "")
cat("targets (full widths): Se ", DELTA_SE, ", Sp ", DELTA_SP, "\n", sep = "")
cat("grid: Se ", length(SE_GRID), " values (", min(SE_GRID), "..", max(SE_GRID),
    "), Sp ", length(SP_GRID), " values (", min(SP_GRID), "..", max(SP_GRID),
    "), prev ", paste(PREV_GRID, collapse = "/"), "\n", sep = "")

# ---- interval tables per arm size n: lo, hi, width, point estimate for x = 0..n
tab_wald <- function(n) {
  x <- 0:n; p <- x / n
  se <- sqrt(p * (1 - p) / n)
  lo <- pmax(0, p - Z * se); hi <- pmin(1, p + Z * se)
  list(lo = lo, hi = hi, w = hi - lo, est = p)
}
tab_wilson <- function(n) {
  x <- 0:n; p <- x / n
  den <- 1 + Z^2 / n
  ctr <- (p + Z^2 / (2 * n)) / den
  mar <- Z * sqrt((p * (1 - p) + Z^2 / (4 * n)) / n) / den
  lo <- pmax(0, ctr - mar); hi <- pmin(1, ctr + mar)
  list(lo = lo, hi = hi, w = hi - lo, est = p)
}
tab_bayes <- function(a, b) function(n) {
  x <- 0:n
  pa <- a + x; pb <- b + n - x
  lo <- qbeta(ALPHA / 2, pa, pb); hi <- qbeta(1 - ALPHA / 2, pa, pb)
  list(lo = lo, hi = hi, w = hi - lo, est = pa / (pa + pb))
}
make_tables <- function(fun, n_max) lapply(1:n_max, fun)   # element n = table for arm size n

# For an arm with tables `tabs` (sizes 1..n_max), target full width `target`,
# and a vector of true values `p`, return (n_max + 1) x length(p) matrices:
#   pw  = P(width <= target + TOL)   (row for n = 0 is 0: empty arm = failure)
#   cov = P(lo <= p <= hi)           (row n = 0: NA)
#   mu  = E[estimate]                (row n = 0: NA)
arm_eval <- function(tabs, target, p) {
  n_max <- length(tabs)
  pw <- cov <- mu <- matrix(NA_real_, n_max + 1L, length(p))
  pw[1L, ] <- 0
  for (n in seq_len(n_max)) {
    tb <- tabs[[n]]; x <- 0:n
    ok <- tb$w <= target + TOL
    for (j in seq_along(p)) {
      pm <- dbinom(x, n, p[j])
      pw[n + 1L, j]  <- sum(pm[ok])
      cov[n + 1L, j] <- sum(pm[tb$lo <= p[j] & p[j] <= tb$hi])
      mu[n + 1L, j]  <- sum(pm * tb$est)
    }
  }
  list(pw = pw, cov = cov, mu = mu)
}
# Same as arm_eval() but only P(width <= target + TOL) (used by the quadrature control).
arm_pw <- function(tabs, target, p) {
  n_max <- length(tabs)
  pw <- matrix(NA_real_, n_max + 1L, length(p)); pw[1L, ] <- 0
  for (n in seq_len(n_max)) {
    x <- 0:n; ok <- which(tabs[[n]]$w <= target + TOL)
    for (j in seq_along(p)) pw[n + 1L, j] <- sum(dbinom(ok - 1L, n, p[j]))
  }
  list(pw = pw)
}


# Operating characteristics of a design over the (Se, Sp, prev) grid.
# Returns a data.frame in wide form.
oc_design <- function(name, N, tab_se, tab_sp, se_grid = SE_GRID, sp_grid = SP_GRID,
                      prev_grid = PREV_GRID) {
  ev_se <- arm_eval(tab_se, DELTA_SE, se_grid)   # rows n = 0..N (arm size of Se arm = K)
  ev_sp <- arm_eval(tab_sp, DELTA_SP, sp_grid)   # arm size of Sp arm = N - K
  K <- 0:N
  out <- list()
  for (pr in prev_grid) {
    pK <- dbinom(K, N, pr)
    # success: sum_K pK * PWse[K, i] * PWsp[N-K, j]
    PWsp_rev <- ev_sp$pw[N - K + 1L, , drop = FALSE]        # row K -> arm size N-K
    succ <- crossprod(ev_se$pw * pK, PWsp_rev)                # length(se_grid) x length(sp_grid)
    # Se-arm quantities conditional on K >= 1
    k1 <- K >= 1L
    p_empty_se <- pK[1L]
    cov_se  <- colSums(ev_se$cov[k1, , drop = FALSE] * pK[k1]) / sum(pK[k1])
    bias_se <- colSums(ev_se$mu[k1, , drop = FALSE]  * pK[k1]) / sum(pK[k1]) - se_grid
    pw_se   <- colSums(ev_se$pw * pK)                         # marginal P(Se width ok); empty = 0
    # Sp-arm quantities conditional on N - K >= 1
    kN <- K <= N - 1L
    p_empty_sp <- pK[N + 1L]
    cov_sp  <- colSums(ev_sp$cov[N - K[kN] + 1L, , drop = FALSE] * pK[kN]) / sum(pK[kN])
    bias_sp <- colSums(ev_sp$mu[N - K[kN] + 1L, , drop = FALSE]  * pK[kN]) / sum(pK[kN]) - sp_grid
    pw_sp   <- colSums(ev_sp$pw[N - K + 1L, , drop = FALSE] * pK)
    for (i in seq_along(se_grid)) for (j in seq_along(sp_grid)) {
      out[[length(out) + 1L]] <- data.frame(
        design = name, N = N, Se = se_grid[i], Sp = sp_grid[j], prev = pr,
        p_success = succ[i, j], p_width_se = pw_se[i], p_width_sp = pw_sp[j],
        cov_se = cov_se[i], cov_sp = cov_sp[j], bias_se = bias_se[i], bias_sp = bias_sp[j],
        p_empty_se_arm = p_empty_se, p_empty_sp_arm = p_empty_sp,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

# ---- 0. Code controls --------------------------------------------------------
cat("\n========== 0. Controls of the interval code ==========\n")
chk_int <- function(n) {
  x <- 0:n
  d1 <- max(abs(tab_wald(n)$lo - pmax(0, x / n - Z * sqrt((x / n) * (1 - x / n) / n))))
  ww <- vapply(x, function(xx) wald_ci(xx, n)[["width"]], numeric(1))
  wl <- vapply(x, function(xx) wilson_ci(xx, n)[["width"]], numeric(1))
  c(wald = max(abs(tab_wald(n)$w - ww)), wilson = max(abs(tab_wilson(n)$w - wl)), d1 = d1)
}
for (n in c(7L, 100L, 116L, 543L)) {
  r <- chk_int(n)
  cat(sprintf("n = %3d: max|width mine - package| Wald %.2e, Wilson %.2e\n", n, r[["wald"]], r[["wilson"]]))
}
N_BUD <- 500L; N_JOINT <- 580L; N_P1 <- 678L; N_P2 <- 672L
bt <- dtasamplesize:::buderer_total_N(0.85, 0.90, 0.07, 0.05, 0.20)
cat("buderer_total_N(0.85, 0.90, 0.07, 0.05, 0.20) = ", bt, " (manuscript Table 2: 500) ",
    if (bt == N_BUD) "PASS" else "**CHECK**", "\n", sep = "")

# ---- 1. Arm tables per design ------------------------------------------------
cat("\n========== 1. Building interval tables ==========\n")
t0 <- Sys.time()
P1 <- list(se = c(17, 3), sp = c(2, 2),  prev = c(4, 16))
P2 <- list(se = c(17, 3), sp = c(18, 2), prev = c(4, 16))
tabs <- list(
  wald_500_se   = make_tables(tab_wald, N_BUD),   wald_500_sp   = make_tables(tab_wald, N_BUD),
  wilson_500_se = make_tables(tab_wilson, N_BUD), wilson_500_sp = make_tables(tab_wilson, N_BUD),
  wilson_580_se = make_tables(tab_wilson, N_JOINT), wilson_580_sp = make_tables(tab_wilson, N_JOINT),
  p1_se = make_tables(tab_bayes(P1$se[1], P1$se[2]), N_P1),
  p1_sp = make_tables(tab_bayes(P1$sp[1], P1$sp[2]), N_P1),
  p2_se = make_tables(tab_bayes(P2$se[1], P2$se[2]), N_P2),
  p2_sp = make_tables(tab_bayes(P2$sp[1], P2$sp[2]), N_P2)
)
cat("tables built in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s\n", sep = "")

# Package convention for an empty arm: P_arm(0) = 0. Check that the posterior-only
# width of the prior meets no target in P1/P2 (so the convention coincides with
# evaluating the width at n = 0).
for (nm in c("P1 Se", "P1 Sp", "P2 Se", "P2 Sp")) {
  pr <- switch(nm, "P1 Se" = P1$se, "P1 Sp" = P1$sp, "P2 Se" = P2$se, "P2 Sp" = P2$sp)
  w0 <- qbeta(1 - ALPHA / 2, pr[1], pr[2]) - qbeta(ALPHA / 2, pr[1], pr[2])
  tg <- if (grepl("Se", nm)) DELTA_SE else DELTA_SP
  cat(sprintf("prior-only width %s Beta(%g,%g) = %.4f vs target %.2f -> meets target: %s (P_arm(0) = 0 coincides)\n",
              nm, pr[1], pr[2], w0, tg, w0 <= tg + TOL))
}

# ---- 2. Operating characteristics over the grid -------------------------------
cat("\n========== 2. Operating characteristics on the grid ==========\n")
t0 <- Sys.time()
res <- rbind(
  oc_design("Buderer_N500_Wald",   N_BUD,   tabs$wald_500_se,   tabs$wald_500_sp),
  oc_design("Buderer_N500_Wilson", N_BUD,   tabs$wilson_500_se, tabs$wilson_500_sp),
  oc_design("Joint_N580_Wilson",   N_JOINT, tabs$wilson_580_se, tabs$wilson_580_sp),
  oc_design("BAM_P1_N678",         N_P1,    tabs$p1_se,         tabs$p1_sp),
  oc_design("BAM_P2_N672",         N_P2,    tabs$p2_se,         tabs$p2_sp)
)
cat("grid evaluated in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
    " s; rows: ", nrow(res), " (", length(unique(res$design)), " designs x ",
    nrow(res) / length(unique(res$design)), " points)\n", sep = "")
stopifnot(all(is.finite(res$p_success)), all(res$p_success >= -1e-15 & res$p_success <= 1 + 1e-12))
write.csv(res, file.path(outdir, "oc_cells_wide.csv"), row.names = FALSE)
mags <- c("p_success", "p_width_se", "p_width_sp", "cov_se", "cov_sp", "bias_se", "bias_sp",
          "p_empty_se_arm", "p_empty_sp_arm")
long <- do.call(rbind, lapply(mags, function(m)
  data.frame(design = res$design, Se = res$Se, Sp = res$Sp, prevalence = res$prev,
             quantity = m, value = res[[m]], stringsAsFactors = FALSE)))
write.csv(long, file.path(outdir, "oc_cells_long.csv"), row.names = FALSE)

# ---- 3. Summaries --------------------------------------------------------------
cat("\n========== 3. Planning point (Se 0.85, Sp 0.90, prev 0.20) ==========\n")
at_plan <- res$Se == PLAN[["Se"]] & res$Sp == PLAN[["Sp"]] & res$prev == PLAN[["prev"]]
pp <- res[at_plan, c("design", "N", mags)]
write.csv(pp, file.path(outdir, "oc_summary_planning_point.csv"), row.names = FALSE)
options(width = 250)
print(transform(pp, p_success = round(p_success, 4), p_width_se = round(p_width_se, 4),
                p_width_sp = round(p_width_sp, 4), cov_se = round(cov_se, 4), cov_sp = round(cov_sp, 4),
                bias_se = round(bias_se, 5), bias_sp = round(bias_sp, 5),
                p_empty_se_arm = signif(p_empty_se_arm, 3), p_empty_sp_arm = signif(p_empty_sp_arm, 3)),
      row.names = FALSE)

cat("\n========== 4. Range over the grid, by design ==========\n")
by_design <- do.call(rbind, lapply(split(res, res$design), function(g) {
  data.frame(design = g$design[1], N = g$N[1], points = nrow(g),
             success_min = min(g$p_success), success_median = median(g$p_success),
             success_max = max(g$p_success),
             frac_success_ge_080 = mean(g$p_success >= 0.80 - TOL),
             cov_se_min = min(g$cov_se), cov_se_max = max(g$cov_se),
             cov_sp_min = min(g$cov_sp), cov_sp_max = max(g$cov_sp),
             bias_se_min = min(g$bias_se), bias_se_max = max(g$bias_se),
             bias_sp_min = min(g$bias_sp), bias_sp_max = max(g$bias_sp),
             max_p_empty_se_arm = max(g$p_empty_se_arm), max_p_empty_sp_arm = max(g$p_empty_sp_arm),
             stringsAsFactors = FALSE)
}))
rownames(by_design) <- NULL
write.csv(by_design, file.path(outdir, "oc_summary_by_design.csv"), row.names = FALSE)
print(transform(by_design, success_min = round(success_min, 4), success_median = round(success_median, 4),
                success_max = round(success_max, 4), frac_success_ge_080 = round(frac_success_ge_080, 3),
                cov_se_min = round(cov_se_min, 4), cov_se_max = round(cov_se_max, 4),
                cov_sp_min = round(cov_sp_min, 4), cov_sp_max = round(cov_sp_max, 4),
                bias_se_min = round(bias_se_min, 4), bias_se_max = round(bias_se_max, 4),
                bias_sp_min = round(bias_sp_min, 4), bias_sp_max = round(bias_sp_max, 4),
                max_p_empty_se_arm = signif(max_p_empty_se_arm, 3),
                max_p_empty_sp_arm = signif(max_p_empty_sp_arm, 3)), row.names = FALSE)

cat("\n--- P(success) vs true Se (Sp 0.90, prev 0.20) ---\n")
sl <- res[res$Sp == 0.90 & res$prev == 0.20, c("design", "Se", "p_success")]
wide_se <- reshape(sl, idvar = "Se", timevar = "design", direction = "wide")
names(wide_se) <- sub("^p_success\\.", "", names(wide_se))
print(round(wide_se, 4), row.names = FALSE)
cat("\n--- P(success) vs true Sp (Se 0.85, prev 0.20) ---\n")
sl <- res[res$Se == 0.85 & res$prev == 0.20, c("design", "Sp", "p_success")]
wide_sp <- reshape(sl, idvar = "Sp", timevar = "design", direction = "wide")
names(wide_sp) <- sub("^p_success\\.", "", names(wide_sp))
print(round(wide_sp, 4), row.names = FALSE)
cat("\n--- P(success) vs true prevalence (Se 0.85, Sp 0.90) ---\n")
sl <- res[res$Se == 0.85 & res$Sp == 0.90, c("design", "prev", "p_success")]
wide_pr <- reshape(sl, idvar = "prev", timevar = "design", direction = "wide")
names(wide_pr) <- sub("^p_success\\.", "", names(wide_pr))
print(round(wide_pr, 4), row.names = FALSE)

# ---- 5. Comparison with published / preliminary values --------------------------
cat("\n========== 5. Comparison with published and preliminary values ==========\n")
cmp <- list()
add_cmp <- function(item, source, published, obtained, tol, note = "") {
  cmp[[length(cmp) + 1L]] <<- data.frame(item = item, source = source, published = published,
    obtained = obtained, diff = obtained - published, tolerance = tol,
    agrees = abs(obtained - published) <= tol, note = note, stringsAsFactors = FALSE)
}
get_pp <- function(d, m) pp[[m]][pp$design == d]
# exact Wald / Wilson P(width <= 0.14) at fixed n = 100, Se = 0.85 (conditioning on K = 100)
ev100 <- arm_eval(tabs$wald_500_se[1:100], DELTA_SE, 0.85)$pw[101, 1]
ev100w <- arm_eval(tabs$wilson_500_se[1:100], DELTA_SE, 0.85)$pw[101, 1]
add_cmp("P(width<=0.14) Wald, n=100 fixed, Se=0.85", "width_attainment_exact.R / plan evidence", 0.5683, ev100, 5e-5,
        "published value rounded to 4 decimals")
add_cmp("P(width<=0.14) Wilson, n=100 fixed, Se=0.85", "width_attainment_exact.R / plan evidence", 0.5683, ev100w, 5e-5,
        "published value rounded to 4 decimals")
# Joint Wilson exact at N = 578 and 577 and 580 (point values)
joint_A <- function(N, se = 0.85, sp = 0.90, prev = 0.20) {
  K <- 0:N; pK <- dbinom(K, N, prev)
  ok_se <- vapply(K, function(k) if (k == 0) 0 else
    { tb <- tab_wilson(k); sum(dbinom(0:k, k, se)[tb$w <= DELTA_SE + TOL]) }, numeric(1))
  ok_sp <- vapply(N - K, function(k) if (k == 0) 0 else
    { tb <- tab_wilson(k); sum(dbinom(0:k, k, sp)[tb$w <= DELTA_SP + TOL]) }, numeric(1))
  sum(pK * ok_se * ok_sp)
}
A578 <- joint_A(578L); A577 <- joint_A(577L); A580 <- joint_A(580L)
add_cmp("Joint Wilson exact A(578)", "comparator_calls.log L62", 0.801244683, A578, 5e-10)
add_cmp("Joint Wilson exact A(577)", "comparator_calls.log L63", 0.7983424783, A577, 5e-11)
add_cmp("Joint Wilson exact A(580) vs this script's planning point", "internal", A580,
        get_pp("Joint_N580_Wilson", "p_success"), 1e-12, "separate code path (scalar loop vs matrix)")
mc_z <- (A580 - 0.8041) / 0.002806
add_cmp("Joint Wilson MC 0.8041 (B = 20000, seed 2026) vs exact A(580)", "article, Table 2", 0.8041, A580, 3 * 0.002806,
        sprintf("(MC - exact)/MCSE = %.2f", (0.8041 - A580) / 0.002806))
# Bayesian A(N) from the package
cv <- bam_assurance_curve(prior_se = P1$se, prior_sp = P1$sp, prior_prev = P1$prev,
                          delta_se = DELTA_SE, delta_sp = DELTA_SP, N_range = c(677L, 678L))
cv2 <- bam_assurance_curve(prior_se = P2$se, prior_sp = P2$sp, prior_prev = P2$prev,
                           delta_se = DELTA_SE, delta_sp = DELTA_SP, N_range = 672L)
A678 <- cv$curve$assurance[cv$curve$N == 678]; A677 <- cv$curve$assurance[cv$curve$N == 677]
A672 <- cv2$curve$assurance[1]
add_cmp("BAM P1 A(678), package", "article, worked example", 0.8003489948, A678, 5e-11)
add_cmp("BAM P1 A(677), package", "article, worked example", 0.7996848824, A677, 5e-11)
add_cmp("BAM P2 A(672), package", "article, Table 2", 0.8002692084, A672, 5e-11)
cmp <- do.call(rbind, cmp)
write.csv(cmp, file.path(outdir, "oc_comparison_published.csv"), row.names = FALSE)
print(transform(cmp, published = signif(published, 10), obtained = signif(obtained, 10),
                diff = signif(diff, 3), tolerance = signif(tolerance, 3)), row.names = FALSE)
cat(sprintf("Exact Joint Wilson A(580) = %.10f ; MC 0.8041 +/- 0.002806 differs by %.2f MCSE\n",
            A580, (0.8041 - A580) / 0.002806))

# ---- 6. Mandatory check: average of P(success) over the design prior = A(N) -------
cat("\n========== 6. Average of P(success) over the DESIGN prior vs A(N) of the package ==========\n")
# (a) Closed form (Beta-Binomial sums): tautological for the Bayesian designs; a code control.
bb_pmf <- function(x, n, a, b) exp(lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b))
avg_bb <- function(N, tab_se, tab_sp, pri) {
  K <- 0:N
  pse <- vapply(K, function(k) if (k == 0) 0 else { tb <- tab_se[[k]]; sum(bb_pmf(0:k, k, pri$se[1], pri$se[2])[tb$w <= DELTA_SE + TOL]) }, numeric(1))
  psp <- vapply(N - K, function(k) if (k == 0) 0 else { tb <- tab_sp[[k]]; sum(bb_pmf(0:k, k, pri$sp[1], pri$sp[2])[tb$w <= DELTA_SP + TOL]) }, numeric(1))
  sum(bb_pmf(K, N, pri$prev[1], pri$prev[2]) * pse * psp)
}
# (b) Gauss-Legendre quadrature of the SAME fixed-truth functions used on the grid
#     (arm_eval): integrates P(success | Se, Sp, prev) against the Beta priors.
#     It exercises the binomial layer of the code, which (a) does not.
gl_nodes <- function(m) {
  i <- seq_len(m - 1L); b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, m, m); J[cbind(i, i + 1L)] <- b; J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = (e$values[o] + 1) / 2, w = e$vectors[1, o]^2)   # nodes/weights on [0,1] (weights sum to 1)
}
GL_M <- 400L
gl <- gl_nodes(GL_M)
cat("Gauss-Legendre nodes: ", GL_M, " (exact for polynomial degree <= ", 2 * GL_M - 1,
    "; integrands have degree <= N + 18 = 696)\n", sep = "")
avg_quad <- function(N, tab_se, tab_sp, pri) {
  ev_se <- arm_pw(tab_se, DELTA_SE, gl$x)
  ev_sp <- arm_pw(tab_sp, DELTA_SP, gl$x)
  w_se <- gl$w * dbeta(gl$x, pri$se[1], pri$se[2])
  w_sp <- gl$w * dbeta(gl$x, pri$sp[1], pri$sp[2])
  w_pr <- gl$w * dbeta(gl$x, pri$prev[1], pri$prev[2])
  K <- 0:N
  Sebar <- as.vector(ev_se$pw %*% w_se)               # by arm size 0..N
  Spbar <- as.vector(ev_sp$pw %*% w_sp)
  Wprev <- vapply(K, function(k) sum(w_pr * dbinom(k, N, gl$x)), numeric(1))
  sum(Wprev * Sebar[K + 1L] * Spbar[N - K + 1L])
}
chk_avg <- data.frame()
for (cfg in list(list(nm = "BAM P1 N=678", N = N_P1, ts = tabs$p1_se, tp = tabs$p1_sp, pri = P1, A = A678),
                 list(nm = "BAM P2 N=672", N = N_P2, ts = tabs$p2_se, tp = tabs$p2_sp, pri = P2, A = A672))) {
  a_bb <- avg_bb(cfg$N, cfg$ts, cfg$tp, cfg$pri)
  a_q  <- avg_quad(cfg$N, cfg$ts, cfg$tp, cfg$pri)
  chk_avg <- rbind(chk_avg, data.frame(design = cfg$nm, A_package = cfg$A,
    avg_betabinomial = a_bb, diff_betabinomial = a_bb - cfg$A,
    avg_quadrature = a_q, diff_quadrature = a_q - cfg$A, stringsAsFactors = FALSE))
}
write.csv(chk_avg, file.path(outdir, "oc_prior_average_check.csv"), row.names = FALSE)
print(chk_avg, row.names = FALSE, digits = 12)
cat("max |difference| over both designs and both methods: ",
    format(max(abs(c(chk_avg$diff_betabinomial, chk_avg$diff_quadrature))), digits = 3), "\n", sep = "")
cat("Check passes if <= 1e-10: ",
    if (max(abs(c(chk_avg$diff_betabinomial, chk_avg$diff_quadrature))) <= 1e-10) "PASS" else "**CHECK**", "\n", sep = "")

# ---- 7. Figure -------------------------------------------------------------------
cat("\n========== 7. Figure ==========\n")
fig <- res[res$Sp == 0.90 & res$prev == 0.20, ]
des <- c("Buderer_N500_Wald", "Buderer_N500_Wilson", "Joint_N580_Wilson", "BAM_P1_N678", "BAM_P2_N672")
lab <- c("Buderer N = 500, Wald", "Buderer N = 500, Wilson", "Joint N = 580, Wilson",
         "BAM N = 678 (P1)", "BAM N = 672 (P2)")
col <- c("grey40", "grey60", "black", "#0072B2", "#D55E00")
lty <- c(2, 3, 1, 1, 1); pch <- c(1, 2, 15, 16, 17)
png(file.path(outdir, "oc_success_vs_Se.png"), width = 7, height = 4.8, units = "in", res = 300)
par(mar = c(4.4, 4.6, 1.0, 0.8))
plot(NA, xlim = range(SE_GRID), ylim = c(0, 1), xlab = "True sensitivity (Sp = 0.90, prevalence = 0.20)",
     ylab = "Exact P(both interval widths meet the targets)", las = 1)
abline(h = 0.80, lty = 3, col = "grey70"); abline(v = 0.85, lty = 3, col = "grey70")
for (k in seq_along(des)) {
  g <- fig[fig$design == des[k], ]
  lines(g$Se, g$p_success, col = col[k], lty = lty[k], lwd = 2)
  points(g$Se, g$p_success, col = col[k], pch = pch[k], cex = 0.8)
}
legend("bottomleft", lab, col = col, lty = lty, pch = pch, lwd = 2, bty = "n", cex = 0.85)
dev.off()
cat("figure written: oc_success_vs_Se.png\n")

cat("\nElapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s\n", sep = "")
cat("Done.\n")
