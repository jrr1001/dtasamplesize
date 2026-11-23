#!/usr/bin/env Rscript
# =============================================================================
# compare_williamson2023.R
#
# Reproduces the sample sizes of Williamson et al. (2023) with the package's
# EXACT Bayesian-assurance calculation, and compares the package with the
# outputs recorded from the authors' Shiny application. No Monte Carlo is used
# anywhere in this script (every number is a finite sum of Beta-Binomial
# probabilities); no random numbers are drawn.
#
# Usage
#   Rscript compare_williamson2023.R <output_dir>
#   The output directory is the ONLY thing written to (created if missing):
#   log (compare_williamson2023.log), CSV tables. Single core.
#   If the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB exist they are
#   prepended to .libPaths() (in that order).
#
# Published figures compared (source of each figure)
#   [M]  Williamson SF, Williams CJ, Lendrem BC, Wilson KJ. Diagn Progn Res
#        2023;7:17. doi:10.1186/s41512-023-00153-1.
#     - Illustrative example (section "Illustrative example"): required N = 321
#       at assurance 80 %; priors Se Beta(25,7), Sp Beta(30,2),
#       prevalence Beta(13.56,122.06); half-width targets 0.10 (Se) and 0.05
#       (Sp) from the "Simulation study: application to COVID-19" paragraph and
#       Table 4 (acceptable TPP). Criterion: one-sided 95 % credible interval,
#       half-width = posterior median minus lower bound. The joint reading is
#       inferred from "the resulting sensitivity and specificity intervals will
#       have sufficient precision with probability 80%".
#     - Figure 3 text (Results): n_c = 10 per group, Se only, assurance 0.8,
#       priors (after pessimistic filter) Beta(9,3), Beta(10,2), Beta(11,1) for
#       lab results n11 = 8, 9, 10 (derived: Beta(1,1) + data; the article does
#       not print them for Se). Published: prevalence 0.05 -> 1630, 1295, 783;
#       prevalence 0.30 -> 182, 154, 101.
#   [APP] Observations of the Shiny application of Williamson et al., 2026-11-14/15
#       (dated screenshots kept by the authors, available on request).
#
# Conventions (named explicitly)
#   * Width criterion: ONE-SIDED, HALF-width = posterior median minus the lower
#     bound qbeta(1 - 0.95, ., .) of the one-sided 95 % credible interval
#     (package: width_type = "lower_half_width", alpha_ci = 0.95). Targets are
#     half-widths. Comparison "width <= target + TOL", TOL = 1e-12 in this
#     script's own enumeration; the package compares exactly (width <= delta).
#   * Objective: JOINT (Se and Sp both) in part 1 and 3; MARGINAL Se
#     (assurance_for = "se") in part 2.
#   * Design prior = analysis prior (as in Williamson et al.): the package's
#     defaults (analysis_prior_* = NULL).
#   * Success rule: assurance A(N) >= 0.80 at the first N of an integer scan
#     (no tolerance, as in the package). Empty-arm convention of the package:
#     an arm with 0 subjects fails (P(0) = 0).
#   * Prevalence prior of the Figure 3 text (NOT published). Declared
#     assumption of this script: mean = rho_T and 95 % equal-tailed
#     interval width = 0.1, derived as Beta(mu*s, (1-mu)*s) with s solving
#     qbeta(.975, mu*s, (1-mu)*s) - qbeta(.025, mu*s, (1-mu)*s) = 0.1
#     (uniroot, tol 1e-12; an earlier exploratory script, not distributed, used the
#     default uniroot tolerance and printed Beta(3.43,65.23) and
#     Beta(96.34,224.79)). Same as Beta(13.56,122.06) of the illustrative
#     example, which has mean 0.1 and width 0.1.
#   * Own enumeration (an independent implementation of the same exact
#     assurance, used for A(N-1), A(N) and for the variants of part 2): per-arm
#     success probability P(m) for m = 0..M diseased (or non-diseased)
#     subjects, M = 400; beyond M, P(m) is taken equal to P(M), which is
#     checked to be 1 within 1e-9 for every scenario (logged), together with
#     the Beta-Binomial mass of m > M.
#
# Part 2 variants tried for the difference of 1-2 units (proposal P7):
#   V0 base; V1 strict success (width < target - TOL); V2a/V2b prevalence
#   hyperparameters rounded to 2 / 1 decimals; V3a/V3b quantiles by uniroot on
#   pbeta (default tolerance 1.2e-4 / tight 1e-12) instead of qbeta; V4 prior
#   moment-matched to mean rho_T and sd = 0.1 / (2 * 1.96); V5 point
#   prevalence (m ~ Bin(N, rho_T), no prior); V6 constant prior size
#   Beta(321 rho, 321 (1 - rho)) (the prior size implied at rho = 0.30).
#   V5 and V6 are exploratory hypotheses of this script, not P7 proposals.
#
# Outputs (in <output_dir>)
#   compare_williamson2023.log
#   f32a_illustrative.csv            part 1: A(320), A(321), A(322), package N
#   f32b_figure3.csv                 part 2: package vs published, per prior
#   f32b_figure3_variants.csv        part 2: variants V0..V6
#   f33_package_vs_app.csv           part 3: app observation vs package
#   f33_restrictions.csv             part 3: restrictions observed in captures
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript compare_williamson2023.R <output_dir>")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

pkg_lib <- Sys.getenv("DTA_PKG_LIB", unset = "")
ext_lib <- Sys.getenv("DTA_EXTRA_LIB", unset = "")
libs <- c(if (nzchar(pkg_lib)) pkg_lib, if (nzchar(ext_lib)) ext_lib)
if (length(libs)) .libPaths(c(libs, .libPaths()))
suppressPackageStartupMessages(library(dtasamplesize))

TOL <- 1e-12
logfile <- file.path(outdir, "compare_williamson2023.log")
logcon <- file(logfile, open = "wt", encoding = "UTF-8")
say <- function(...) {
  s <- paste0(..., collapse = "")
  cat(s, "\n", sep = "")
  cat(s, "\n", sep = "", file = logcon)
  flush(logcon)
}
fmt <- function(x, d = 7) formatC(x, digits = d, format = "f")
tstart <- Sys.time()
say("compare_williamson2023.R")
say("date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("packageVersion(dtasamplesize): ", as.character(packageVersion("dtasamplesize")))
say("find.package(dtasamplesize): ", find.package("dtasamplesize"))
say("R.version.string: ", R.version.string)
say("TOL = ", format(TOL), "; output_dir = ", outdir)
say("")

# ---- helpers ----------------------------------------------------------------
bbpmf <- function(k, n, a, b) exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))

q_alg <- function(p, A, B, alg) {
  if (alg == "qbeta") return(qbeta(p, A, B))
  tol <- if (alg == "uniroot_default") .Machine$double.eps^0.25 else 1e-12
  mapply(function(a, b) uniroot(function(q) pbeta(q, a, b) - p, c(0, 1), tol = tol)$root, A, B)
}

# per-arm success probability P(m), m = 0..M, one-sided 95 % half-width <= delta
p_arm <- function(M, a, b, delta, strict = FALSE, alg = "qbeta") {
  P <- numeric(M + 1L)
  for (m in 1:M) {
    x <- 0:m
    A <- a + x; B <- b + m - x
    w <- q_alg(0.5, A, B, alg) - q_alg(0.05, A, B, alg)
    ok <- if (strict) w < delta - TOL else w <= delta + TOL
    P[m + 1L] <- sum(bbpmf(x, m, a, b)[ok])
  }
  P
}

# marginal assurance A(N) = sum_m pmf(m; N) P(m); pmf Beta-Binomial or Binomial
a_marg <- function(N, P, ap = NULL, bp = NULL, rho = NULL) {
  M <- length(P) - 1L
  k <- 0:min(M, N)
  pm <- if (is.null(rho)) bbpmf(k, N, ap, bp) else dbinom(k, N, rho)
  sum(pm * P[k + 1L]) + (1 - sum(pm)) * P[M + 1L]
}
mass_above <- function(N, M, ap, bp, rho = NULL) {
  k <- 0:min(M, N)
  pm <- if (is.null(rho)) bbpmf(k, N, ap, bp) else dbinom(k, N, rho)
  1 - sum(pm)
}

# joint assurance A(N) = sum_m pmf(m) Pse(m) Psp(N - m)   (Beta-Binomial m)
a_joint <- function(N, Pse, Psp, ap, bp) {
  m <- 0:N
  pm <- bbpmf(m, N, ap, bp)
  sum(pm * Pse[m + 1L] * Psp[N - m + 1L])
}

first_N <- function(Ns, Afun) {
  for (N in Ns) {
    a <- Afun(N)
    if (a >= 0.80) return(c(N = N, A = a))
  }
  c(N = NA_real_, A = NA_real_)
}

# Beta(mu*s, (1-mu)*s) with equal-tailed 95 % interval width `wid`
fit_beta <- function(mu, wid = 0.1, tol = 1e-12) {
  f <- function(s) diff(qbeta(c(.025, .975), mu * s, (1 - mu) * s)) - wid
  s <- uniroot(f, c(2, 1e4), tol = tol)$root
  c(mu * s, (1 - mu) * s)
}
prev_ci_width <- function(ab) diff(qbeta(c(.025, .975), ab[1], ab[2]))

# package call, warnings captured and counted
pkg_call <- function(...) {
  w <- character()
  r <- withCallingHandlers(
    bam_sample_size(..., method = "exact", width_type = "lower_half_width",
                    alpha_ci = 0.95),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  list(res = r, warnings = w)
}

# =============================================================================
# PART 1  illustrative example, N = 321
# =============================================================================
say("=== PART 1: illustrative example [M, Illustrative example, Table 3(i), Table 4] ===")
t1 <- Sys.time()
ill <- list(prior_se = c(25, 7), prior_sp = c(30, 2), prior_prev = c(13.56, 122.06),
            delta_se = 0.10, delta_sp = 0.05)
say("inputs: Se Beta(25,7); Sp Beta(30,2); prevalence Beta(13.56,122.06); ",
    "one-sided 95% half-widths 0.10 (Se), 0.05 (Sp); assurance 0.80; joint objective")
pk <- pkg_call(prior_se = ill$prior_se, prior_sp = ill$prior_sp, prior_prev = ill$prior_prev,
               delta_se = ill$delta_se, delta_sp = ill$delta_sp, assurance_for = "joint")
say("package bam_sample_size(): N_total = ", pk$res$N_total,
    "; joint_assurance = ", fmt(pk$res$joint_assurance, 7),
    "; warnings issued: ", length(pk$warnings))
cv <- bam_assurance_curve(prior_se = ill$prior_se, prior_sp = ill$prior_sp,
                          prior_prev = ill$prior_prev, delta_se = ill$delta_se,
                          delta_sp = ill$delta_sp, width_type = "lower_half_width",
                          alpha_ci = 0.95, assurance_for = "joint", N_range = 250:400)
say("package bam_assurance_curve(250:400): N_first = ", cv$N_first,
    "; n_decreases = ", cv$n_decreases, "; monotone = ", cv$monotone)
# independent enumeration
Pse_i <- p_arm(400, 25, 7, 0.10); Psp_i <- p_arm(400, 30, 2, 0.05)
Aown <- sapply(320:322, a_joint, Pse = Pse_i, Psp = Psp_i, ap = 13.56, bp = 122.06)
Apkg <- cv$curve$assurance[cv$curve$N %in% 320:322]
f32a <- data.frame(
  N = 320:322, A_package = Apkg, A_own_enumeration = Aown,
  abs_diff = abs(Apkg - Aown),
  published_N = c(NA, 321, NA),
  source = c("", "[M] Illustrative example (required sample size 321)", ""),
  plan_preliminary_A = c(0.7985633, 0.8007552, 0.8029264),
  plan_source = "earlier exploratory calculation (not distributed)")
f32a$package_N_first <- cv$N_first
f32a$package_N_total_bam_sample_size <- pk$res$N_total
f32a$diff_N_vs_published <- pk$res$N_total - 321L
write.csv(f32a, file.path(outdir, "f32a_illustrative.csv"), row.names = FALSE)
for (i in 1:3) say(sprintf("  N = %d: A package = %s; A own enumeration = %s; |diff| = %.2e; plan preliminary = %s",
                           f32a$N[i], fmt(f32a$A_package[i], 10), fmt(f32a$A_own_enumeration[i], 10),
                           f32a$abs_diff[i], fmt(f32a$plan_preliminary_A[i], 7)))
say("  published N = 321; package N = ", pk$res$N_total, "; difference = ", pk$res$N_total - 321L)
say("  part 1 elapsed: ", round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), " s")
say("")

# =============================================================================
# PART 2  Figure 3 text, n_c = 10, Se only
# =============================================================================
say("=== PART 2: Figure 3 text [M, Results; Figure 3] ===")
say("inputs: Se only (assurance_for = 'se'); one-sided 95% half-width 0.10; assurance 0.80;")
say("design = analysis prior Beta(9,3), Beta(10,2), Beta(11,1); prevalence prior: declared assumption")
say("(mean rho_T, 95% equal-tailed width 0.1); NOT published in [M].")
t2 <- Sys.time()
pubs <- list(`0.05` = c(1630, 1295, 783), `0.3` = c(182, 154, 101))
priors_se <- list(c(9, 3), c(10, 2), c(11, 1))
M <- 400L
rows <- list(); vrows <- list()
for (rho in c(0.05, 0.30)) {
  pp <- fit_beta(rho)
  say(sprintf("rho_T = %.2f: prevalence prior Beta(%.10f, %.10f); mean = %.6f; 95%% width = %.8f",
              rho, pp[1], pp[2], pp[1] / sum(pp), prev_ci_width(pp)))
  pub <- pubs[[as.character(rho)]]
  Nmax <- if (rho < 0.1) 2600L else 700L
  for (j in 1:3) {
    pr <- priors_se[[j]]
    # package
    t0 <- Sys.time()
    pk2 <- pkg_call(prior_se = pr, prior_prev = pp, delta_se = 0.10, assurance_for = "se",
                    N_max = 4000L)
    tp <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    # own enumeration, base variant (V0)
    P0 <- p_arm(M, pr[1], pr[2], 0.10)
    f0 <- first_N(10:Nmax, function(N) a_marg(N, P0, pp[1], pp[2]))
    Nm1 <- f0["N"] - 1
    Am1 <- a_marg(Nm1, P0, pp[1], pp[2])
    say(sprintf("  Se Beta(%d,%d): package N = %s (A = %s, %.0f s); own enumeration N = %d (A(N-1) = %s, A(N) = %s); published %d",
                pr[1], pr[2], pk2$res$N_total, fmt(pk2$res$joint_assurance, 7), tp,
                f0["N"], fmt(Am1, 7), fmt(f0["A"], 7), pub[j]))
    say(sprintf("    P_se(M) = %.12f; Beta-Binomial mass of m > M at N = %d: %.3e; warnings: %d",
                P0[M + 1L], f0["N"], mass_above(f0["N"], M, pp[1], pp[2]), length(pk2$warnings)))
    expl <- if (pk2$res$N_total == pub[j]) "reproduced exactly" else
      "see f32b_figure3_variants.csv and the report"
    rows[[length(rows) + 1L]] <- data.frame(
      rho_T = rho, prior_se = sprintf("Beta(%d,%d)", pr[1], pr[2]),
      prev_prior = sprintf("Beta(%.6f,%.6f)", pp[1], pp[2]),
      published_N = pub[j], published_source = "[M] Results, text on Figure 3",
      plan_preliminary_N = if (rho < 0.1) c(1628, 1294, 782)[j] else pub[j],
      plan_source = "earlier exploratory calculation (not distributed)",
      package_N = pk2$res$N_total, package_A = pk2$res$joint_assurance,
      own_N = unname(f0["N"]), own_A_at_N = unname(f0["A"]), own_A_at_Nminus1 = unname(Am1),
      diff_package_minus_published = pk2$res$N_total - pub[j],
      P_se_at_M = P0[M + 1L], package_seconds = round(tp, 1), note = expl)
  }
}
f32b <- do.call(rbind, rows)
write.csv(f32b, file.path(outdir, "f32b_figure3.csv"), row.names = FALSE)
say("  part 2 (package + base enumeration) elapsed: ",
    round(as.numeric(difftime(Sys.time(), t2, units = "secs")), 1), " s")
say("")

# ---- variants for the 1-2 unit differences ------------------------------------
say("--- variants for the difference at rho_T = 0.05 (P7) and consistency check at rho_T = 0.30 ---")
t2v <- Sys.time()
variant_defs <- list(
  V0  = list(desc = "base: qbeta, non-strict, unrounded prevalence prior"),
  V1  = list(desc = "strict success (width < target)", strict = TRUE),
  V2a = list(desc = "prevalence hyperparameters rounded to 2 decimals", round = 2),
  V2b = list(desc = "prevalence hyperparameters rounded to 1 decimal", round = 1),
  V3a = list(desc = "quantiles by uniroot on pbeta, default tolerance (1.2e-4)", alg = "uniroot_default"),
  V3b = list(desc = "quantiles by uniroot on pbeta, tolerance 1e-12", alg = "uniroot_tight"),
  V4  = list(desc = "prevalence prior moment-matched, mean rho_T, sd = 0.1/(2*1.96)", mm = TRUE),
  V5  = list(desc = "point prevalence: m ~ Bin(N, rho_T), no prior", point = TRUE),
  V6  = list(desc = "constant prior size: Beta(321 rho_T, 321 (1 - rho_T))", n321 = TRUE))
Pcache <- list()
for (rho in c(0.05, 0.30)) {
  pp0 <- fit_beta(rho)
  pub <- pubs[[as.character(rho)]]
  Nmax <- if (rho < 0.1) 2600L else 700L
  for (vn in names(variant_defs)) {
    v <- variant_defs[[vn]]
    ppv <- pp0
    if (!is.null(v$round)) ppv <- round(pp0, v$round)
    if (isTRUE(v$mm)) {
      sd0 <- 0.1 / (2 * 1.96); s <- rho * (1 - rho) / sd0^2 - 1; ppv <- c(rho * s, (1 - rho) * s)
    }
    if (isTRUE(v$n321)) ppv <- c(321 * rho, 321 * (1 - rho))
    alg <- if (is.null(v$alg)) "qbeta" else v$alg
    strict <- isTRUE(v$strict)
    Ns <- integer(3); As <- numeric(3)
    for (j in 1:3) {
      pr <- priors_se[[j]]
      key <- paste(alg, strict, pr[1], pr[2])
      if (is.null(Pcache[[key]])) Pcache[[key]] <- p_arm(M, pr[1], pr[2], 0.10, strict = strict, alg = alg)
      P <- Pcache[[key]]
      f <- if (isTRUE(v$point)) first_N(10:Nmax, function(N) a_marg(N, P, rho = rho))
           else first_N(10:Nmax, function(N) a_marg(N, P, ppv[1], ppv[2]))
      Ns[j] <- f["N"]; As[j] <- f["A"]
    }
    say(sprintf("  rho_T = %.2f %-3s N = %s (published %s); A = %s | %s",
                rho, vn, paste(Ns, collapse = ", "), paste(pub, collapse = ", "),
                paste(fmt(As, 6), collapse = ", "), v$desc))
    vrows[[length(vrows) + 1L]] <- data.frame(
      rho_T = rho, variant = vn, description = v$desc,
      prev_prior = if (isTRUE(v$point)) "point (Binomial)" else sprintf("Beta(%.6f,%.6f)", ppv[1], ppv[2]),
      N_Se_Beta_9_3 = Ns[1], N_Se_Beta_10_2 = Ns[2], N_Se_Beta_11_1 = Ns[3],
      published_9_3 = pub[1], published_10_2 = pub[2], published_11_1 = pub[3],
      diff_9_3 = Ns[1] - pub[1], diff_10_2 = Ns[2] - pub[2], diff_11_1 = Ns[3] - pub[3],
      reproduces_all_three = all(Ns == pub))
  }
}
f32bv <- do.call(rbind, vrows)
write.csv(f32bv, file.path(outdir, "f32b_figure3_variants.csv"), row.names = FALSE)
say("  variants elapsed: ", round(as.numeric(difftime(Sys.time(), t2v, units = "secs")), 1), " s")
say("")

# =============================================================================
# PART 3  (package side): configurations captured in the Shiny application
# =============================================================================
say("=== PART 3 (package side): application configurations [APP] ===")
say("Package convention for all rows: joint, one-sided 95% lower half-width, design = analysis prior,")
say("exact. The application does not state its criterion (no field for credible level or sidedness in the")
say("Calculator, 2026-11-14_2353_calculator_default.png); the comparison below tests the one-sided")
say("half-width reading, which is the convention of [M].")
t3 <- Sys.time()
cfg <- list(
  A = list(se = c(25, 7),  sp = c(30, 2), prev = c(13.56, 122.06), dse = 0.10, dsp = 0.05),
  B = list(se = c(17, 3),  sp = c(18, 2), prev = c(4, 16),         dse = 0.07, dsp = 0.05),
  C = list(se = c(17, 3),  sp = c(2, 2),  prev = c(4, 16),         dse = 0.07, dsp = 0.05),
  D = list(se = c(17, 3),  sp = c(1, 1),  prev = c(4, 16),         dse = 0.07, dsp = 0.05),
  E = list(se = c(17, 3),  sp = c(1, 1),  prev = c(4, 16),         dse = 0.07, dsp = 0.01))
pkg_cfg <- function(cn, maxn) {
  k <- cfg[[cn]]
  p <- pkg_call(prior_se = k$se, prior_sp = k$sp, prior_prev = k$prev,
                delta_se = k$dse, delta_sp = k$dsp, assurance_for = "joint",
                N_max = as.integer(maxn))
  list(N = p$res$N_total, A = p$res$joint_assurance, reached = p$res$target_reached,
       maxA = p$res$max_assurance_evaluated, NatMax = p$res$N_at_max_assurance,
       nwarn = length(p$warnings))
}
curve_cfg <- function(cn, Ns) {
  k <- cfg[[cn]]
  cvv <- bam_assurance_curve(prior_se = k$se, prior_sp = k$sp, prior_prev = k$prev,
                             delta_se = k$dse, delta_sp = k$dsp, width_type = "lower_half_width",
                             alpha_ci = 0.95, assurance_for = "joint", N_range = Ns)
  setNames(cvv$curve$assurance, cvv$curve$N)
}
# application observations: config, run label, max.n, N, A (NA = not found), capture
obs <- data.frame(
  config = c("A", "B", "B", "B", "C", "C", "C", "C", "C", "D", "E"),
  run = c("A-1", "B-2", "B-3", "B-4", "C-1", "C-2", "C-3", "C-4", "C-5", "D-1 (Sp Beta(1,1), Both)",
          "E-1 (Sp Beta(1,1), Sp half-width 0.01)"),
  max_n = c(500, 1000, 1000, 1000, 1000, 1500, 1500, 1000, 1000, 1000, 1000),
  app_N = c(321, 552, 551, 551, 551, 552, 552, 552, 552, 551, NA),
  app_A = c(0.8008, 0.8007, 0.8005, 0.8005, 0.8005, 0.8007, 0.8007, 0.8007, 0.8007, 0.8005, NA),
  capture = c("2026-11-14_2354_A_Both_calculator_output.png",
              "2026-11-15_0003_B_Both_maxn1000_rep2_output.png",
              "2026-11-15_0004_B_Both_maxn1000_rep3_output.png",
              "no screenshot (read as text)",
              "2026-11-14_2358_C_Both_maxn1000_output.png",
              "2026-11-14_2359_C_Both_maxn1500_output.png",
              "no screenshot",
              "2026-11-15_0001_C_Both_maxn1000_rep2_output.png",
              "no screenshot",
              "2026-11-15_0006_diagnostic_Sp_Beta11_Both_output.png",
              "2026-11-15_0007_diagnostic_Sp_Beta11_SpTarget001_output.png"),
  stringsAsFactors = FALSE)
# config in the app run (the one whose inputs were typed): B-2 was typed as Sp Beta(18,2)
obs$config_typed <- c("A", "B", "B", "B", "C", "C", "C", "C", "C", "D", "E")
obs$Se_prior <- vapply(obs$config_typed, function(c) sprintf("Beta(%d,%d)", cfg[[c]]$se[1], cfg[[c]]$se[2]), "")
obs$Sp_prior <- vapply(obs$config_typed, function(c) sprintf("Beta(%g,%g)", cfg[[c]]$sp[1], cfg[[c]]$sp[2]), "")
obs$prev_prior <- vapply(obs$config_typed, function(c) sprintf("Beta(%g,%g)", cfg[[c]]$prev[1], cfg[[c]]$prev[2]), "")
obs$halfwidth_Se <- vapply(obs$config_typed, function(c) cfg[[c]]$dse, 0)
obs$halfwidth_Sp <- vapply(obs$config_typed, function(c) cfg[[c]]$dsp, 0)
obs$assurance_target <- 0.8

pk_res <- list(); curv <- list()
for (cn in c("A", "B", "C", "D", "E")) {
  mx <- if (cn == "A") 500 else 1000
  t0 <- Sys.time()
  pk_res[[cn]] <- pkg_cfg(cn, mx)
  say(sprintf("  config %s, package (N_max = %d): N_total = %s; A = %s; target_reached = %s; max A evaluated = %s at N = %s; %.0f s",
              cn, mx, pk_res[[cn]]$N, fmt(pk_res[[cn]]$A, 7), pk_res[[cn]]$reached,
              fmt(pk_res[[cn]]$maxA, 7), pk_res[[cn]]$NatMax,
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
# package with the larger cap used for C-2 and C-3 (max.n = 1500)
pk_C1500 <- pkg_cfg("C", 1500)
say(sprintf("  config C, package (N_max = 1500): N_total = %s; A = %s", pk_C1500$N, fmt(pk_C1500$A, 7)))
# A at N = 551, 552 for B, C, D; 321 for A
for (cn in c("B", "C", "D")) curv[[cn]] <- curve_cfg(cn, 550:553)
curv[["A"]] <- curve_cfg("A", 320:322)
for (cn in c("A", "B", "C", "D")) {
  say("  A(N) package, config ", cn, ": ",
      paste(sprintf("N=%s %s", names(curv[[cn]]), fmt(curv[[cn]], 7)), collapse = "; "))
}
# alternative convention: two-sided equal-tailed full width = 2 x half-width (manuscript convention)
say("  alternative (manuscript) convention, two-sided full width = 2 x the half-width typed in the app:")
alt <- list()
for (cn in c("B", "C", "D")) {
  k <- cfg[[cn]]
  r <- bam_sample_size(prior_se = k$se, prior_sp = k$sp, prior_prev = k$prev,
                       delta_se = 2 * k$dse, delta_sp = 2 * k$dsp, width_type = "full_width",
                       alpha_ci = 0.95, assurance_for = "joint", method = "exact", N_max = 1000L)
  alt[[cn]] <- c(N = r$N_total, A = r$joint_assurance)
  say(sprintf("    config %s: N_total = %s; A = %s", cn, r$N_total, fmt(r$joint_assurance, 7)))
}

pkgN <- vapply(obs$config_typed, function(c) as.numeric(pk_res[[c]]$N), 0)
pkgA <- vapply(obs$config_typed, function(c) as.numeric(pk_res[[c]]$A), 0)
# C-2 and C-3 were run with max.n = 1500: same package N (checked above)
obs$package_N <- pkgN
obs$package_A_at_N <- pkgA
obs$package_A_at_app_N <- NA_real_
for (i in seq_len(nrow(obs))) {
  c <- obs$config_typed[i]
  if (!is.na(obs$app_N[i]) && !is.null(curv[[c]]) && as.character(obs$app_N[i]) %in% names(curv[[c]]))
    obs$package_A_at_app_N[i] <- curv[[c]][[as.character(obs$app_N[i])]]
}
obs$diff_N_app_minus_package <- obs$app_N - obs$package_N
obs$A_app_vs_package_round4 <- ifelse(is.na(obs$package_A_at_app_N), NA,
                                      round(obs$package_A_at_app_N, 4))
obs$A_agrees_4dp <- ifelse(is.na(obs$app_A) | is.na(obs$A_app_vs_package_round4), NA,
                           obs$app_A == obs$A_app_vs_package_round4)
obs$N_agrees <- ifelse(is.na(obs$app_N), is.na(obs$package_N), obs$app_N == obs$package_N)
obs$alt_package_N_two_sided_full_width <- vapply(obs$config_typed, function(c)
  if (c %in% names(alt)) unname(alt[[c]]["N"]) else NA_real_, 0)
obs$comment <- c(
  "app = package (N and A to 4 decimals)",
  "app N differs from package for the typed priors (552 vs 551); equals package N for Sp Beta(2,2) (config C); typed prior was edited just before",
  "app = package", "app = package",
  "app N differs from package for the typed priors (551 vs 552); equals package N for Sp Beta(18,2) (config B), the previous Sp prior",
  "app = package", "app = package", "app = package", "app = package",
  "N equal; app A 0.8005 differs from package A(551) = 0.8001 at 4 decimals; 0.8005 equals package A(551) for Sp Beta(18,2) (config B)",
  "app returns NA and an error; package reaches no N <= 1000 (target_reached = FALSE)")
write.csv(obs, file.path(outdir, "f33_package_vs_app.csv"), row.names = FALSE)
say("  table written: f33_package_vs_app.csv (", nrow(obs), " application observations)")
for (i in seq_len(nrow(obs)))
  say(sprintf("    %-40s app N = %s A = %s | package N = %s A(app N) = %s | N agree: %s; A agree (4 dp): %s",
              obs$run[i], obs$app_N[i], obs$app_A[i], obs$package_N[i],
              ifelse(is.na(obs$package_A_at_app_N[i]), "NA", fmt(obs$package_A_at_app_N[i], 4)),
              obs$N_agrees[i], obs$A_agrees_4dp[i]))
# hypothesis check, logged as such: the application output for Sp-prior edits matches the package
# evaluated with the Sp prior typed BEFORE the edit on the first run after an edit
say("  Observation (not an explanation): the 3 application runs whose N or A differ from the package at the")
say("  typed priors (B-2, C-1, D-1) all agree with the package at the Sp prior that was entered immediately")
say("  before; the repeat runs agree with the typed priors. Not established: see the report (capture requested).")

# restrictions observed in captures
restr <- data.frame(
  restriction = c(
    "Only 'Both' operative", "Maximum N control 1 to 1000, default 500", "Server accepted max.n = 1500",
    "No credible-level field in the Calculator", "No one-/two-sided choice in the Calculator",
    "No separate design and analysis priors in the Calculator", "No batch or multi-scenario input",
    "Result when no N is found within max.n", "Repeatability",
    "Assurance printed with 4 decimals", "Prior input methods (Se, Sp)",
    "Prevalence prior input", "Assurance curve plot", "Version string", "Package-side behaviour of the same inputs"),
  observation = c(
    "Text 'Currently only both sensitivity and specificity targets can be calculated using the input option'; Sensitivity-only and Specificity-only runs returned 'An error has occurred'",
    "HTML range 1 to 1000, default 500", "Calculated with max.n = 1500 (HTML control declares maximum 1000)",
    "Calculator has Set Target, two half-widths, Required Assurance, max n, Find n only",
    "Not in the Calculator; the Analysis tab has alpha (0.05, range 0 to 0.2) and One/Two tailed (default One)",
    "Not observed in the Calculator; Analysis tab is a separate tab with its own Beta(A,B) prior",
    "Not observed; one scenario per run; a 'Generate report' button exists in the Sample size tab (not tried)",
    "'is NA, which provides an assurance of NA' followed by 'An error has occurred'",
    "Same inputs gave 551 or 552 (B, C), and an error in the first B run",
    "All outputs", "Previous Results (2x2), Elicitation, Enter Parameters",
    "Menu with Previous Data, Interval, Enter Parameters; only the alpha-beta table was visible in the options tried",
    "Curve of assurance against n with dashed lines at 0.8 and at N", "'MIC Calculator v0.03'; Info cites Wilson et al. (2022), not Williamson et al. (2023)",
    "Package exact search: see f33_package_vs_app.csv; package accepts any N_max, any prior, sidedness and level through arguments"),
  capture = c(
    "2026-11-14_2353_calculator_default.png; 2026-11-14_2355_A_SensitivityOnly_output.png; 2026-11-14_2355_A_SpecificityOnly_output.png",
    "2026-11-14_2353_calculator_default.png", "2026-11-14_2359_C_Both_maxn1500_output.png",
    "2026-11-14_2353_calculator_default.png", "2026-11-14_2353_calculator_default.png; 2026-11-15_0008_analysis_tab.png",
    "2026-11-14_2353_calculator_default.png; 2026-11-15_0008_analysis_tab.png",
    "2026-11-14_2353_calculator_default.png",
    "2026-11-15_0007_diagnostic_Sp_Beta11_SpTarget001_output.png",
    "2026-11-14_2357_B_Both_maxn1000_output.png; 2026-11-15_0003_B_Both_maxn1000_rep2_output.png; 2026-11-15_0004_B_Both_maxn1000_rep3_output.png; 2026-11-14_2358_C_Both_maxn1000_output.png; 2026-11-15_0001_C_Both_maxn1000_rep2_output.png",
    "all output captures", "2026-11-14_2350_assurance_prior_method_dropdown.png; 2026-11-14_2351_assurance_enter_parameters_default.png; 2026-11-14_2351_assurance_elicitation_default.png",
    "2026-11-14_2352_assurance_prevalence_prior_default.png; 2026-11-15_0008_assurance_prevalence_previous_data.png",
    "2026-11-14_2354_A_Both_calculator_output.png", "2026-11-14_2350_landing_home.png; 2026-11-14_2350_info_tab.png",
    "this script"),
  status = c("observed", "observed", "observed", "observed", "observed (absent in Calculator)",
             "observed (absent in Calculator)", "observed (absent; 'Generate report' not tried: [not verified])",
             "observed", "observed", "observed", "observed", "observed (partly: [not verified] which of the three options uses the hidden inputs)",
             "observed", "observed", "computed"),
  stringsAsFactors = FALSE)
write.csv(restr, file.path(outdir, "f33_restrictions.csv"), row.names = FALSE)
say("  restrictions table written: f33_restrictions.csv (", nrow(restr), " rows)")
say("  part 3 elapsed: ", round(as.numeric(difftime(Sys.time(), t3, units = "secs")), 1), " s")
say("")
say("total elapsed: ", round(as.numeric(difftime(Sys.time(), tstart, units = "secs")), 1), " s")
close(logcon)
