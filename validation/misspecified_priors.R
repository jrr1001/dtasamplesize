#!/usr/bin/env Rscript
# =============================================================================
# misspecified_priors.R
#
# EXACT (finite-sum) behaviour of the Bayesian designs of the manuscript's
# worked example when the TRUTH lies away from the centre of the design prior,
# and of the post-study prior-data conflict check of bam_prior_data_conflict().
# For each design and each declared scenario of true (Se, Sp, prev):
#   * probability that both credible intervals meet their full-width targets,
#   * frequentist coverage of the Se and Sp credible intervals, and bias of
#     the posterior mean,
#   * probability that the conflict check FLAGS a conflict (detection
#     probability): the study data (K diseased of N; TP of K; TN of N - K) are
#     compared with the Beta-Binomial prior predictive of the ANALYSIS prior
#     (inferential reference) and of the DESIGN prior (calibration).
# The two frequentist comparators (Buderer N = 500 Wald, joint N = 580 Wilson)
# are evaluated at the same truths for success and coverage ("when applicable":
# they have no prior, so no conflict check).
#
# Usage
#   Rscript misspecified_priors.R <output_dir>
#   Writes ONLY into <output_dir>: log, CSV tables. If DTA_PKG_LIB and
#   DTA_EXTRA_LIB exist they are prepended to .libPaths() (in that order).
#   Single core; no Monte Carlo (every number is a finite sum).
#
# INPUTS (from the 0.6.9 article)
#   Planning values Se 0.85, Sp 0.90, prevalence 0.20, level 0.95.
#   Full-width targets 0.14 (Se) and 0.10 (Sp).
#   P1 (N = 678): Se ~ Beta(17,3), Sp ~ Beta(2,2), prev ~ Beta(4,16);
#       validation/compat_v069_check.R.
#   P2 (N = 672): Sp ~ Beta(18,2), the rest as P1.
#   Frequentist N = 500 (Buderer, Wald) and N = 580 (joint, Wilson): Table 2 of the article.
#   In P1 and P2 the analysis prior equals the design prior (as in the
#   manuscript), so the "analysis" and "design" references of the check
#   coincide there. To exercise the two references, one ILLUSTRATIVE variant
#   is added: P2w = P2 design priors and N = 672, with analysis priors of half
#   the strength and the same means (Se ~ Beta(8.5,1.5), Sp ~ Beta(9,1));
#   these analysis priors are an assumption of this script, not a published
#   value, and N = 672 is not re-optimised for them.
#
# SCENARIOS OF THE TRUTH (declared; centre of the prior: Se 0.85, Sp 0.90 in P2
# and 0.50 in P1, prevalence 0.20)
#   S00 (0.85, 0.90, 0.20)  planning values (well specified for Se, prev, P2)
#   S01 (0.75, 0.90, 0.20)  Se 0.10 below the Se prior mean
#   S02 (0.65, 0.90, 0.20)  Se 0.20 below
#   S03 (0.85, 0.80, 0.20)  Sp 0.10 below the Sp prior mean of P2
#   S04 (0.85, 0.70, 0.20)  Sp 0.20 below
#   S05 (0.85, 0.90, 0.35)  prevalence above the prior centre (prior sd 0.087)
#   S06 (0.85, 0.90, 0.08)  prevalence below
#   S07 (0.95, 0.90, 0.20)  Se 0.10 above
#   S08 (0.75, 0.80, 0.35)  Se, Sp and prevalence all off
#   S09 (0.65, 0.70, 0.35)  all off, severe
#   PRIOR_PREDICTIVE        truth drawn from the DESIGN prior: the detection
#                           probability is the false-alarm rate of the check
#                           (calibration) and the success probability is A(N).
#   For P1 the Sp prior is vague (mean 0.5): the planned truth Sp = 0.90 lies
#   0.4 above its mean, which is the manuscript's own case. The log prints,
#   per design, the distance of each truth from the prior mean in prior-sd
#   units (columns z_se, z_sp, z_prev).
#
# CONFLICT CHECK AS EXACT SUMS
#   The default of bam_prior_data_conflict() is used: alpha = 0.05 with a
#   1e-12 tolerance. For x_obs of n: p_lower = P(X <= x_obs), p_upper =
#   P(X >= x_obs), p_two_sided = min(1, 2 min(p_lower, p_upper)), X ~
#   BetaBinomial(n, a, b); a tail is flagged if p <= 0.05 + 1e-12 and the
#   two-sided value if <= 0.05 + 1e-12. Quantities: Se (TP of K diseased),
#   Sp (TN of N - K non-diseased), prevalence (K of N against the prevalence
#   prior; design reference only). Declared decision rules:
#     "two_sided": the two-sided flag of the quantity (primary rule);
#     "any_tail" : flag_lower OR flag_upper (union of the one-tail rules).
#   Quantity-level detection = P(flag); "any" = P(at least one of Se, Sp,
#   prevalence is flagged). Given K the Se and Sp outcomes are independent, so
#   P(any | K) = 1 - (1 - qSe(K)) (1 - qSp(N-K)) (1 - 1[prev flagged at K]),
#   summed over K ~ Binomial(N, prev). An empty arm (K = 0 or N - K = 0) has
#   no row in the package ("row omitted"): it contributes no flag for that
#   quantity. Success of the design with an empty arm is failure (P_arm(0) =
#   0). Coverage and bias are conditional on a non-empty arm, as in
#   operating_characteristics.R.
#   The closed form is checked against bam_prior_data_conflict() on 20 random
#   outcomes (p-values to 1e-12, flags identical), on every outcome of two
#   arm sizes, and for the P2w two-reference output.
#
# Conventions as in operating_characteristics.R: full-width targets, Wald
# truncated to [0,1] (degenerate Wald width 0 meets the target), Wilson as
# wilson_ci(), equal-tailed Beta posterior intervals with the analysis prior,
# width <= target + TOL with TOL = 1e-12, coverage with closed intervals,
# K ~ Binomial(N, prev).
#
# OUTPUTS (in <output_dir>)
#   misspecified_priors.log
#   mp_cells_long.csv      (design, Se, Sp, prevalence, quantity, value, scenario)
#   mp_cells_wide.csv      one row per (design, scenario)
#   mp_summary.csv         compact table: success, coverage, detection
#   mp_prior_average_check.csv   average over the design prior vs A(N)
#   mp_conflict_closed_form_check.csv
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript misspecified_priors.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pkg_lib   <- Sys.getenv("DTA_PKG_LIB",   unset = "")
extra_lib <- Sys.getenv("DTA_EXTRA_LIB", unset = "")
newlibs <- c(if (nzchar(pkg_lib)) pkg_lib, if (nzchar(extra_lib)) extra_lib)
if (length(newlibs)) .libPaths(c(newlibs, .libPaths()))

logfile <- file.path(outdir, "misspecified_priors.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)
t_start <- Sys.time()

suppressMessages(library(dtasamplesize))
cat("misspecified_priors.R\n")
cat("date:                 ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R:                    ", R.version.string, "\n")
cat("dtasamplesize:        ", as.character(packageVersion("dtasamplesize")), "\n")
cat("find.package:         ", normalizePath(find.package("dtasamplesize"), winslash = "/"), "\n")
cat("output dir:           ", normalizePath(outdir, winslash = "/"), "\n")

TOL <- 1e-12
ALPHA <- 0.05                # level of the interval AND flagging threshold of the check
Z <- qnorm(1 - ALPHA / 2)
FLAG_ALPHA <- 0.05; FLAG_TOL <- 1e-12   # defaults of bam_prior_data_conflict()
DELTA_SE <- 0.14; DELTA_SP <- 0.10
cat("TOL = ", format(TOL), "; interval level 0.95; targets (full widths) Se ", DELTA_SE, ", Sp ", DELTA_SP,
    "; check alpha = ", FLAG_ALPHA, " (+ ", format(FLAG_TOL), ")\n", sep = "")

# ---- interval tables (as in operating_characteristics.R) ---------------------
tab_wald <- function(n) {
  x <- 0:n; p <- x / n; se <- sqrt(p * (1 - p) / n)
  lo <- pmax(0, p - Z * se); hi <- pmin(1, p + Z * se)
  list(lo = lo, hi = hi, w = hi - lo, est = p)
}
tab_wilson <- function(n) {
  x <- 0:n; p <- x / n; den <- 1 + Z^2 / n
  ctr <- (p + Z^2 / (2 * n)) / den
  mar <- Z * sqrt((p * (1 - p) + Z^2 / (4 * n)) / n) / den
  lo <- pmax(0, ctr - mar); hi <- pmin(1, ctr + mar)
  list(lo = lo, hi = hi, w = hi - lo, est = p)
}
tab_bayes <- function(a, b) function(n) {
  x <- 0:n; pa <- a + x; pb <- b + n - x
  lo <- qbeta(ALPHA / 2, pa, pb); hi <- qbeta(1 - ALPHA / 2, pa, pb)
  list(lo = lo, hi = hi, w = hi - lo, est = pa / (pa + pb))
}
make_tables <- function(fun, n_max) lapply(1:n_max, fun)
arm_eval <- function(tabs, target, p) {
  n_max <- length(tabs)
  pw <- cov <- mu <- matrix(NA_real_, n_max + 1L, length(p))
  pw[1L, ] <- 0
  for (n in seq_len(n_max)) {
    tb <- tabs[[n]]; x <- 0:n; ok <- tb$w <= target + TOL
    for (j in seq_along(p)) {
      pm <- dbinom(x, n, p[j])
      pw[n + 1L, j]  <- sum(pm[ok])
      cov[n + 1L, j] <- sum(pm[tb$lo <= p[j] & p[j] <= tb$hi])
      mu[n + 1L, j]  <- sum(pm * tb$est)
    }
  }
  list(pw = pw, cov = cov, mu = mu)
}

# ---- Beta-Binomial predictive tails and flags (exact) --------------------------
bb_pmf <- function(x, n, a, b) exp(lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b))
# Flags for every x_obs = 0..n against BetaBinomial(n, a, b). Same arithmetic as
# the package: pl = min(1, sum of first x+1 masses), pu = min(1, sum of last masses).
bb_flags <- function(n, a, b) {
  pm <- bb_pmf(0:n, n, a, b)
  pl <- pmin(1, cumsum(pm)); pu <- pmin(1, rev(cumsum(rev(pm))))
  p2 <- pmin(1, 2 * pmin(pl, pu))
  lim <- FLAG_ALPHA + FLAG_TOL
  list(pl = pl, pu = pu, p2 = p2, two = p2 <= lim, tail = (pl <= lim) | (pu <= lim))
}
# list over arm sizes 1..n_max of flag vectors for the prior (a, b)
make_flags <- function(n_max, ab) lapply(1:n_max, function(n) bb_flags(n, ab[1], ab[2]))
# P(flag) for arm size n and truth p: matrix (n_max+1) x length(p); row n = 0 -> 0 (row omitted)
arm_flagprob <- function(fl, which, p) {
  n_max <- length(fl)
  out <- matrix(0, n_max + 1L, length(p))
  for (n in seq_len(n_max)) {
    f <- fl[[n]][[which]]; x <- 0:n
    for (j in seq_along(p)) out[n + 1L, j] <- sum(dbinom(x, n, p[j])[f])
  }
  out
}

# ---- designs --------------------------------------------------------------------
P1  <- list(se = c(17, 3), sp = c(2, 2),  prev = c(4, 16), an_se = c(17, 3),   an_sp = c(2, 2))
P2  <- list(se = c(17, 3), sp = c(18, 2), prev = c(4, 16), an_se = c(17, 3),   an_sp = c(18, 2))
P2w <- list(se = c(17, 3), sp = c(18, 2), prev = c(4, 16), an_se = c(8.5, 1.5), an_sp = c(9, 1))
N_BUD <- 500L; N_JOINT <- 580L; N_P1 <- 678L; N_P2 <- 672L

SC <- data.frame(
  scenario = sprintf("S%02d", 0:9),
  Se   = c(0.85, 0.75, 0.65, 0.85, 0.85, 0.85, 0.85, 0.95, 0.75, 0.65),
  Sp   = c(0.90, 0.90, 0.90, 0.80, 0.70, 0.90, 0.90, 0.90, 0.80, 0.70),
  prev = c(0.20, 0.20, 0.20, 0.20, 0.20, 0.35, 0.08, 0.20, 0.35, 0.35),
  desc = c("planning values", "Se -0.10", "Se -0.20", "Sp -0.10", "Sp -0.20", "prev 0.35 (above)",
           "prev 0.08 (below)", "Se +0.10", "Se, Sp, prev off", "Se, Sp, prev off (severe)"),
  stringsAsFactors = FALSE)
cat("\nScenarios of the truth:\n"); print(SC, row.names = FALSE)
beta_mean <- function(ab) ab[1] / sum(ab)
beta_sd   <- function(ab) sqrt(ab[1] * ab[2] / (sum(ab)^2 * (sum(ab) + 1)))
cat("\nPrior means and sd: Se Beta(17,3) mean ", beta_mean(P1$se), " sd ", round(beta_sd(P1$se), 4),
    "; Sp P1 Beta(2,2) mean ", beta_mean(P1$sp), " sd ", round(beta_sd(P1$sp), 4),
    "; Sp P2 Beta(18,2) mean ", beta_mean(P2$sp), " sd ", round(beta_sd(P2$sp), 4),
    "; prev Beta(4,16) mean ", beta_mean(P1$prev), " sd ", round(beta_sd(P1$prev), 4), "\n", sep = "")

cat("\n========== 1. Building tables ==========\n")
t0 <- Sys.time()
tb_wald_se <- make_tables(tab_wald, N_BUD);   tb_wald_sp <- tb_wald_se
tb_wil_se  <- make_tables(tab_wilson, N_JOINT); tb_wil_sp <- tb_wil_se
bayes_tabs <- function(des, N) list(se = make_tables(tab_bayes(des$an_se[1], des$an_se[2]), N),
                                    sp = make_tables(tab_bayes(des$an_sp[1], des$an_sp[2]), N))
bt_P1 <- bayes_tabs(P1, N_P1); bt_P2 <- bayes_tabs(P2, N_P2); bt_P2w <- bayes_tabs(P2w, N_P2)
# conflict flags: reference "analysis" (analysis prior) and "design" (design prior)
mk_fl <- function(des, N) list(
  an_se = make_flags(N, des$an_se), an_sp = make_flags(N, des$an_sp),
  de_se = make_flags(N, des$se),    de_sp = make_flags(N, des$sp),
  prev  = bb_flags(N, des$prev[1], des$prev[2]))
fl_P1 <- mk_fl(P1, N_P1); fl_P2 <- mk_fl(P2, N_P2); fl_P2w <- mk_fl(P2w, N_P2)
cat("tables built in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s\n", sep = "")

# ---- 2. Closed form vs bam_prior_data_conflict() ------------------------------------
cat("\n========== 2. Closed-form tails and flags vs bam_prior_data_conflict() ==========\n")
call_pkg <- function(tp, fn, tn, fp, des) {
  suppressWarnings(bam_prior_data_conflict(tp, fn, tn, fp, prior_se = des$se, prior_sp = des$sp,
    prior_prev = des$prev, analysis_prior_se = des$an_se, analysis_prior_sp = des$an_sp))
}
mine_rows <- function(tp, fn, tn, fp, des, fl, N) {
  K <- tp + fn; M <- fp + tn
  g <- function(f, x) c(pl = f$pl[x + 1L], pu = f$pu[x + 1L], p2 = f$p2[x + 1L],
                        two = f$two[x + 1L], tail_lo = NA, tail_up = NA)
  list(prevalence = g(fl$prev, K),
       se_analysis = g(fl$an_se[[K]], tp), se_design = g(fl$de_se[[K]], tp),
       sp_analysis = g(fl$an_sp[[M]], tn), sp_design = g(fl$de_sp[[M]], tn))
}
cmp_one <- function(tp, fn, tn, fp, des, fl, N, label) {
  pk <- call_pkg(tp, fn, tn, fp, des); mine <- mine_rows(tp, fn, tn, fp, des, fl, N)
  rows <- list()
  for (i in seq_len(nrow(pk))) {
    key <- if (pk$reference[i] == "prevalence") "prevalence" else
      paste0(if (grepl("^sens", pk$quantity[i])) "se_" else "sp_", pk$reference[i])
    m <- mine[[key]]
    # when the package returns a single row (priors coincide) it is the analysis row
    rows[[length(rows) + 1L]] <- data.frame(case = label, key = key,
      d_pl = abs(pk$p_lower[i] - m[["pl"]]), d_pu = abs(pk$p_upper[i] - m[["pu"]]),
      d_p2 = abs(pk$p_two_sided[i] - m[["p2"]]),
      flag_two_pkg = pk$flag_two_sided[i],
      flag_two_same = identical(pk$flag_two_sided[i], unname(m[["two"]] == 1)),
      flag_lo_same = identical(pk$flag_lower[i], as.logical(pk$p_lower[i] <= FLAG_ALPHA + FLAG_TOL)),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}
# helper that also checks my "tail" rule against the package flags
tail_same <- function(tp, fn, tn, fp, des, fl) {
  pk <- call_pkg(tp, fn, tn, fp, des); K <- tp + fn; M <- fp + tn
  ok <- TRUE
  for (i in seq_len(nrow(pk))) {
    f <- if (pk$reference[i] == "prevalence") fl$prev else {
      se <- grepl("^sens", pk$quantity[i])
      if (pk$reference[i] == "analysis") (if (se) fl$an_se[[K]] else fl$an_sp[[M]])
      else (if (se) fl$de_se[[K]] else fl$de_sp[[M]]) }
    x <- pk$x_obs[i]
    mine_tail <- f$tail[x + 1L]
    pkg_tail <- isTRUE(pk$flag_lower[i]) || isTRUE(pk$flag_upper[i])
    ok <- ok && identical(mine_tail, pkg_tail)
  }
  ok
}
set.seed(20261115)
chk <- list()
cat("20 random outcomes per design (seed 20261115; K, TP, TN uniform on their ranges):\n")
for (cfg in list(list(nm = "P1", des = P1, fl = fl_P1, N = N_P1),
                 list(nm = "P2", des = P2, fl = fl_P2, N = N_P2),
                 list(nm = "P2w", des = P2w, fl = fl_P2w, N = N_P2))) {
  N <- cfg$N
  for (r in 1:20) {
    K <- sample(1:(N - 1L), 1L); tp <- sample(0:K, 1L); tn <- sample(0:(N - K), 1L)
    c1 <- cmp_one(tp, K - tp, tn, N - K - tn, cfg$des, cfg$fl, N, paste0(cfg$nm, "_rand", r))
    c1$tail_same <- tail_same(tp, K - tp, tn, N - K - tn, cfg$des, cfg$fl)
    chk[[length(chk) + 1L]] <- c1
  }
}
# Draws from the neighbourhood of the tails, where flags change (exhaustive over x for fixed other counts)
for (cfg in list(list(nm = "P2", des = P2, fl = fl_P2, N = N_P2, K = 134L),
                 list(nm = "P1", des = P1, fl = fl_P1, N = N_P1, K = 136L),
                 list(nm = "P2w", des = P2w, fl = fl_P2w, N = N_P2, K = 134L))) {
  N <- cfg$N; K <- cfg$K
  tn0 <- round(0.9 * (N - K))
  for (tp in 0:K) {
    c1 <- cmp_one(tp, K - tp, tn0, N - K - tn0, cfg$des, cfg$fl, N, paste0(cfg$nm, "_allTP_K", K))
    c1$tail_same <- tail_same(tp, K - tp, tn0, N - K - tn0, cfg$des, cfg$fl)
    chk[[length(chk) + 1L]] <- c1
  }
  tp0 <- round(0.85 * K)
  for (tn in 0:(N - K)) {
    c1 <- cmp_one(tp0, K - tp0, tn, N - K - tn, cfg$des, cfg$fl, N, paste0(cfg$nm, "_allTN_K", K))
    c1$tail_same <- tail_same(tp0, K - tp0, tn, N - K - tn, cfg$des, cfg$fl)
    chk[[length(chk) + 1L]] <- c1
  }
}
chk <- do.call(rbind, chk)
write.csv(chk, file.path(outdir, "mp_conflict_closed_form_check.csv"), row.names = FALSE)
cat("rows compared (quantity-reference rows): ", nrow(chk), " from ", length(unique(chk$case)), " groups\n", sep = "")
cat("max |p_lower diff| = ", format(max(chk$d_pl), digits = 3), "; max |p_upper diff| = ", format(max(chk$d_pu), digits = 3),
    "; max |p_two_sided diff| = ", format(max(chk$d_p2), digits = 3), "\n", sep = "")
cat("two-sided flags identical in ", sum(chk$flag_two_same), "/", nrow(chk), "; lower-tail flags identical in ",
    sum(chk$flag_lo_same), "/", nrow(chk), "; any-tail rule identical in ", sum(chk$tail_same), "/", nrow(chk), "\n", sep = "")
cat("two-sided flag values in the compared rows (package): TRUE ", sum(chk$flag_two_pkg), ", FALSE ", sum(!chk$flag_two_pkg), " (shows the comparison is not vacuous)\n", sep = "")
cat("closed form agrees with bam_prior_data_conflict() to 1e-12 and flags identical: ",
    if (max(c(chk$d_pl, chk$d_pu, chk$d_p2)) <= 1e-12 && all(chk$flag_two_same) && all(chk$flag_lo_same) && all(chk$tail_same))
      "PASS" else "**CHECK**", "\n", sep = "")

# ---- 3. Operating characteristics per design and scenario -------------------------------
cat("\n========== 3. Success, coverage, bias and detection per design and scenario ==========\n")
oc_scen <- function(name, N, tab_se, tab_sp, prev_vec = SC$prev, se_vec = SC$Se, sp_vec = SC$Sp,
                    fl = NULL, des = NULL) {
  ev_se <- arm_eval(tab_se, DELTA_SE, se_vec)
  ev_sp <- arm_eval(tab_sp, DELTA_SP, sp_vec)
  K <- 0:N
  if (!is.null(fl)) {
    fq <- list(
      an_se_two = arm_flagprob(fl$an_se, "two", se_vec),  an_se_tail = arm_flagprob(fl$an_se, "tail", se_vec),
      de_se_two = arm_flagprob(fl$de_se, "two", se_vec),  de_se_tail = arm_flagprob(fl$de_se, "tail", se_vec),
      an_sp_two = arm_flagprob(fl$an_sp, "two", sp_vec),  an_sp_tail = arm_flagprob(fl$an_sp, "tail", sp_vec),
      de_sp_two = arm_flagprob(fl$de_sp, "two", sp_vec),  de_sp_tail = arm_flagprob(fl$de_sp, "tail", sp_vec))
    prev_two  <- fl$prev$two[K + 1L]; prev_tail <- fl$prev$tail[K + 1L]
  }
  out <- list()
  for (s in seq_len(nrow(SC))) {
    pK <- dbinom(K, N, prev_vec[s])
    k1 <- K >= 1L; kN <- K <= N - 1L
    row <- data.frame(design = name, N = N, scenario = SC$scenario[s], Se = se_vec[s], Sp = sp_vec[s], prev = prev_vec[s],
      p_success = sum(pK * ev_se$pw[K + 1L, s] * ev_sp$pw[N - K + 1L, s]),
      cov_se = sum(pK[k1] * ev_se$cov[K[k1] + 1L, s]) / sum(pK[k1]),
      cov_sp = sum(pK[kN] * ev_sp$cov[N - K[kN] + 1L, s]) / sum(pK[kN]),
      bias_se = sum(pK[k1] * ev_se$mu[K[k1] + 1L, s]) / sum(pK[k1]) - se_vec[s],
      bias_sp = sum(pK[kN] * ev_sp$mu[N - K[kN] + 1L, s]) / sum(pK[kN]) - sp_vec[s],
      stringsAsFactors = FALSE)
    if (!is.null(fl)) {
      # z distances from the prior means
      row$z_se <- (se_vec[s] - beta_mean(des$se)) / beta_sd(des$se)
      row$z_sp <- (sp_vec[s] - beta_mean(des$sp)) / beta_sd(des$sp)
      row$z_prev <- (prev_vec[s] - beta_mean(des$prev)) / beta_sd(des$prev)
      for (rf in c("an", "de")) for (rule in c("two", "tail")) {
        qse <- fq[[paste0(rf, "_se_", rule)]][K + 1L, s]            # K = 0 -> 0 (row omitted)
        qsp <- fq[[paste0(rf, "_sp_", rule)]][N - K + 1L, s]
        fpr <- as.numeric(if (rule == "two") prev_two else prev_tail)
        nm_rule <- if (rule == "two") "two_sided" else "any_tail"
        nm_ref <- if (rf == "an") "analysis" else "design"
        row[[paste0("det_se_", nm_rule, "_", nm_ref)]]  <- sum(pK * qse)
        row[[paste0("det_sp_", nm_rule, "_", nm_ref)]]  <- sum(pK * qsp)
        row[[paste0("det_any_", nm_rule, "_", nm_ref)]] <- sum(pK * (1 - (1 - qse) * (1 - qsp) * (1 - fpr)))
        if (rf == "an") row[[paste0("det_prev_", nm_rule)]] <- sum(pK * fpr)
      }
    }
    out[[s]] <- row
  }
  do.call(rbind, out)
}
t0 <- Sys.time()
r_P1  <- oc_scen("BAM_P1_N678",  N_P1, bt_P1$se,  bt_P1$sp,  fl = fl_P1,  des = P1)
r_P2  <- oc_scen("BAM_P2_N672",  N_P2, bt_P2$se,  bt_P2$sp,  fl = fl_P2,  des = P2)
r_P2w <- oc_scen("BAM_P2w_N672_weaker_analysis_prior", N_P2, bt_P2w$se, bt_P2w$sp, fl = fl_P2w, des = P2w)
r_bud <- oc_scen("Buderer_N500_Wald", N_BUD, tb_wald_se, tb_wald_sp)
r_jnt <- oc_scen("Joint_N580_Wilson", N_JOINT, tb_wil_se, tb_wil_sp)
cat("scenarios evaluated in ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s\n", sep = "")

# ---- 4. Prior-predictive (truth ~ design prior): false-alarm rate and A(N) -----------------
cat("\n========== 4. Truth drawn from the design prior (calibration) and check against A(N) ==========\n")
prior_pred <- function(name, N, des, bt, fl) {
  K <- 0:N
  pK <- bb_pmf(K, N, des$prev[1], des$prev[2])
  # arm-level P under the DESIGN prior: BB(n, a_d, b_d) over x
  arm_bb <- function(tabs, target, ab_d, flags) {
    n_max <- length(tabs)
    pw <- numeric(n_max + 1L)
    q <- lapply(flags, function(f) numeric(n_max + 1L))
    for (n in seq_len(n_max)) {
      pm <- bb_pmf(0:n, n, ab_d[1], ab_d[2])
      pw[n + 1L] <- sum(pm[tabs[[n]]$w <= (if (identical(target, "se")) DELTA_SE else DELTA_SP) + TOL])
      for (nm in names(flags)) q[[nm]][n + 1L] <- sum(pm[flags[[nm]][[n]]])
    }
    list(pw = pw, q = q)
  }
  fse <- list(an_two = lapply(fl$an_se, `[[`, "two"), an_tail = lapply(fl$an_se, `[[`, "tail"),
              de_two = lapply(fl$de_se, `[[`, "two"), de_tail = lapply(fl$de_se, `[[`, "tail"))
  fsp <- list(an_two = lapply(fl$an_sp, `[[`, "two"), an_tail = lapply(fl$an_sp, `[[`, "tail"),
              de_two = lapply(fl$de_sp, `[[`, "two"), de_tail = lapply(fl$de_sp, `[[`, "tail"))
  A_se <- arm_bb(bt$se, "se", des$se, fse); A_sp <- arm_bb(bt$sp, "sp", des$sp, fsp)
  row <- data.frame(design = name, N = N, scenario = "PRIOR_PREDICTIVE", Se = NA_real_, Sp = NA_real_, prev = NA_real_,
                    p_success = sum(pK * A_se$pw[K + 1L] * A_sp$pw[N - K + 1L]), stringsAsFactors = FALSE)
  for (rf in c("an", "de")) for (rule in c("two", "tail")) {
    qse <- A_se$q[[paste0(rf, "_", rule)]][K + 1L]; qsp <- A_sp$q[[paste0(rf, "_", rule)]][N - K + 1L]
    fpr <- as.numeric(if (rule == "two") fl$prev$two[K + 1L] else fl$prev$tail[K + 1L])
    nm_rule <- if (rule == "two") "two_sided" else "any_tail"; nm_ref <- if (rf == "an") "analysis" else "design"
    row[[paste0("det_se_", nm_rule, "_", nm_ref)]]  <- sum(pK * qse)
    row[[paste0("det_sp_", nm_rule, "_", nm_ref)]]  <- sum(pK * qsp)
    row[[paste0("det_any_", nm_rule, "_", nm_ref)]] <- sum(pK * (1 - (1 - qse) * (1 - qsp) * (1 - fpr)))
    if (rf == "an") row[[paste0("det_prev_", nm_rule)]] <- sum(pK * fpr)
  }
  row
}
pp_P1  <- prior_pred("BAM_P1_N678", N_P1, P1, bt_P1, fl_P1)
pp_P2  <- prior_pred("BAM_P2_N672", N_P2, P2, bt_P2, fl_P2)
pp_P2w <- prior_pred("BAM_P2w_N672_weaker_analysis_prior", N_P2, P2w, bt_P2w, fl_P2w)

cvA <- function(des, N) bam_assurance_curve(prior_se = des$se, prior_sp = des$sp, prior_prev = des$prev,
  delta_se = DELTA_SE, delta_sp = DELTA_SP, N_range = N,
  analysis_prior_se = if (identical(des$an_se, des$se)) NULL else des$an_se,
  analysis_prior_sp = if (identical(des$an_sp, des$sp)) NULL else des$an_sp)$curve$assurance[1]
avg_chk <- data.frame(
  design = c("BAM_P1_N678", "BAM_P2_N672", "BAM_P2w_N672_weaker_analysis_prior"),
  A_package = c(cvA(P1, N_P1), cvA(P2, N_P2), cvA(P2w, N_P2)),
  average_over_design_prior = c(pp_P1$p_success, pp_P2$p_success, pp_P2w$p_success),
  source = c("article, worked example: 0.8003489948", "article, Table 2: 0.8002692084", "no published value (illustrative)"),
  stringsAsFactors = FALSE)
avg_chk$difference <- avg_chk$average_over_design_prior - avg_chk$A_package
write.csv(avg_chk, file.path(outdir, "mp_prior_average_check.csv"), row.names = FALSE)
print(avg_chk, row.names = FALSE, digits = 12)
cat("max |difference| = ", format(max(abs(avg_chk$difference)), digits = 3),
    "  (control: must be <= 1e-10) ", if (max(abs(avg_chk$difference)) <= 1e-10) "PASS" else "**CHECK**", "\n", sep = "")
cat("published A(678) = 0.8003489948 and A(672) = 0.8002692084: package values above agree to ",
    format(max(abs(avg_chk$A_package[1] - 0.8003489948), abs(avg_chk$A_package[2] - 0.8002692084)), digits = 3), "\n", sep = "")

# ---- 5. Assemble and write ----------------------------------------------------------------
all_cols <- unique(unlist(lapply(list(r_P1, r_P2, r_P2w, r_bud, r_jnt, pp_P1, pp_P2, pp_P2w), names)))
fill <- function(d) { for (m in setdiff(all_cols, names(d))) d[[m]] <- NA_real_; d[, all_cols] }
wide <- rbind(fill(r_P1), fill(r_P2), fill(r_P2w), fill(r_bud), fill(r_jnt), fill(pp_P1), fill(pp_P2), fill(pp_P2w))
write.csv(wide, file.path(outdir, "mp_cells_wide.csv"), row.names = FALSE)
mag <- setdiff(names(wide), c("design", "N", "scenario", "Se", "Sp", "prev"))
long <- do.call(rbind, lapply(mag, function(m) {
  d <- data.frame(design = wide$design, Se = wide$Se, Sp = wide$Sp, prevalence = wide$prev,
                  quantity = m, value = wide[[m]], scenario = wide$scenario, stringsAsFactors = FALSE)
  d[!is.na(d$value), ]
}))
write.csv(long, file.path(outdir, "mp_cells_long.csv"), row.names = FALSE)
cat("\nwritten: mp_cells_wide.csv (", nrow(wide), " rows), mp_cells_long.csv (", nrow(long), " rows)\n", sep = "")

# ---- 6. Printed tables -----------------------------------------------------------------------
options(width = 260)
r4 <- function(d, cols, k = 4) { d[cols] <- lapply(d[cols], function(v) round(v, k)); d }
cat("\n========== 6a. Success probability, coverage, bias (all designs) ==========\n")
print(r4(wide[, c("design", "scenario", "Se", "Sp", "prev", "p_success", "cov_se", "cov_sp", "bias_se", "bias_sp")],
         c("p_success", "cov_se", "cov_sp", "bias_se", "bias_sp")), row.names = FALSE)
cat("\n========== 6b. Detection probability of the conflict check (two-sided rule; analysis = reference, design = calibration) ==========\n")
dd <- wide[grepl("^BAM", wide$design), c("design", "scenario", "z_se", "z_sp", "z_prev",
  "det_se_two_sided_analysis", "det_sp_two_sided_analysis", "det_prev_two_sided", "det_any_two_sided_analysis",
  "det_any_two_sided_design", "det_any_any_tail_analysis", "det_any_any_tail_design")]
print(r4(dd, setdiff(names(dd), c("design", "scenario")), 3), row.names = FALSE)
cat("\nNote: in P1 and P2 the analysis and design priors coincide, so *_analysis equals *_design; they differ only in P2w.\n")
sm <- wide[, c("design", "scenario", "Se", "Sp", "prev", "p_success", "cov_se", "cov_sp",
  "det_se_two_sided_analysis", "det_sp_two_sided_analysis", "det_prev_two_sided", "det_any_two_sided_analysis",
  "det_any_two_sided_design", "det_any_any_tail_analysis", "det_any_any_tail_design")]
write.csv(sm, file.path(outdir, "mp_summary.csv"), row.names = FALSE)

cat("\nElapsed: ", round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1), " s\n", sep = "")
cat("Done.\n")
