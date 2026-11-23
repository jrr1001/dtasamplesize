#!/usr/bin/env Rscript
# =============================================================================
# success_set_shortcut.R
#
# Purpose
#   Compare three rules for the SUCCESS SET of one arm (the set of observed
#   counts x = 0..n for which the posterior credible interval is narrow
#   enough) in the scenarios of the article (P1, P2, C01-C12) and in the
#   scenarios of Wilson et al. (2022, Stat Med 41:2908; VAP example, section
#   6.1 with rho = 0.10 [and, as an extra, rho = 0.05], Supplement B.4):
#
#   (1) PRINTED rule: success iff n11 <= c1 or n11 >= c2, with
#       c1 = argmin{n11 : W(n11) >= w*} - 1 and
#       c2 = argmax{n11 : W(n11) >= w*} + 1
#       (Wilson et al., sections 3.1.1 and 3.1.3), conditional assurance 0
#       or 1 outside n1 < n_T,1 < n2. Equivalent to: success set = the
#       complement of the "hull" [min F, max F] of the failure set F.
#       It coincides with the enumeration if and only if F is contiguous.
#   (2) CODE rule: the success set as programmed in Data S2 of Wilson et al.
#       (Paper Functions.R, cutoffs() and assurance()), EMULATED here from the
#       written description of that rule kept by the authors. NO Data S2
#       code was opened, copied, adapted or executed for this script.
#   (3) ENUMERATION: every count is evaluated, success iff W(x) <= w* + TOL
#       (the package's exact method, R/bam_sample_size.R).
#
# Usage
#   Rscript success_set_shortcut.R <output_dir>
#   (writes ONLY in <output_dir>, created if missing: log, CSV tables.)
#   Libraries: if the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB
#   exist they are prepended to .libPaths() (in that order).
#   Single core; exact enumeration (no Monte Carlo, no random numbers).
#
# Reusable function (callable on its own, with explicit arguments)
#   After  options(dta.success_set.no_main = TRUE); source("success_set_shortcut.R")
#   (this defines the functions and does NOT run the analysis):
#     emulate_data_s2_assurance(N, w, design_measure, analysis_measure,
#                               design_prev, level = 0.95, tol = 1e-12,
#                               rule = c("code", "printed", "enumeration"))
#   returns A(N) (vector over N) for the ONE-ARM (marginal) objective of the
#   sensitivity ("case = 1" of Data S2), with
#     design_measure   = c(a, b)  design prior of the measure (Beta),
#     analysis_measure = c(a, b)  analysis prior of the measure (Beta),
#     design_prev      = c(a, b)  design prior of the prevalence (Beta).
#   Also: emulate_data_s2_arm_prob(n, ...) = P(success | arm of size n).
#   The output CSV wilson_A_by_N.csv holds A(N) with 15 significant digits.
#
# Conventions (named explicitly)
#   * w* is the FULL width (upper minus lower quantile) of the two-sided
#     equal-tailed Beta posterior interval at level `level` (0.95 unless the
#     scenario says otherwise: C07 and C10 use 0.90). Not a half-width and
#     not one-sided. Objective: JOINT (Se and Sp, same replication) for P1,
#     P2 and C01-C12 (each arm with its own rule, n_d ~ BetaBinomial(N,
#     prior_prev), n_nd = N - n_d); MARGINAL (Se only) for the Wilson cases.
#   * Design prior generates the counts (Beta-Binomial predictive); ANALYSIS
#     prior fixes the interval widths W(x). They coincide everywhere except
#     Supplement B.4 "Combination/Hybrid II" (flat analysis Beta(1,1) for Se).
#   * Regla 11: widths are compared with TOL = 1e-12 in the same way in the
#     three rules: failure iff W > w* + TOL. The printed rule defines failure
#     as W >= w* (success W < w*); the code and the package use success W <=
#     w*. To isolate the effect of the SET CONSTRUCTION, all three rules use the
#     same failure predicate (W > w* + TOL); the number of exact ties
#     (|W - w*| <= TOL) is counted and reported (a tie is where the >= / >
#     incoherence would matter).
#   * Empty arm (n = 0): failure in the three rules (package convention
#     P_arm(0) = 0; the code starts its sums at k = 1; for the printed rule it
#     is the declared convention, since A(n | n_T,1 = 0) is not defined there).
#   * N* = first integer N >= 2 with A(N) >= target (as the package; no
#     tolerance on A). Data S2 uses strict ">" (assurance.opt2): ties on A are
#     not expected. Non-monotone A(N) after the first crossing is counted.
#   * Code rule, as emulated (from that written description):
#       c1 = (first i in 1..n with W(i) > w) - 1, NA if none (W(0) is never
#            evaluated);
#       c2 = (last i of the R sequence 1:(n-1) with W(i) > w) + 1, n if none
#            (for n = 1, 1:(n-1) is c(1, 0) in R and the LAST element in
#            iteration order wins; this edge case is an assumption, [not
#            verified against the code]);
#       success: c1 NA -> 0..n; c1 > 1 -> 0..c1, plus c2..n only if c2 < n;
#                c1 in {0, 1} -> c2..n only if c2 < n; (else none).
#     Deviations from the enumeration are classified by cause per count x:
#       lost successes: b (lower tail dropped, c1 <= 1), c (x = n dropped,
#       c2 = n), e (interior successes inside [c1, c2], i.e. F not
#       contiguous); spurious successes: d (x = 0 not evaluated), f (x = n
#       not evaluated, c2 < n), g (other). Single-cause variants are the
#       enumeration with ONLY that cause applied; the printed rule is
#       enumeration minus the interior successes of the hull of F.
#   * Case = 2 of Data S2 (joint mode with the Sp arm of size k) is NOT
#     emulated here; the joint objective uses the paper's formulation.
#
# Published figures compared (source)
#   VAP n* = 106, A(150) = 0.88: Wilson et al. 2022, section 5.1 [W];
#   B.4 Table 1: sizes 106 (Standard, Hybrid I) and 290 (Combination, Hybrid
#   II): Data S1, B.4, L268-292 [WS1]; section 6.1 rho = 0.10: 681, rho = 0.05:
#   1643 (rho = 0.01: 2770, not computed here): section 6.1 [W]. Preliminary
#   values (earlier exploratory calculation, not evidence): code rule 106, 290, 681,
#   1643; exact enumeration 104, 290, 661, 1408.
#
# Outputs (in <output_dir>)
#   success_set_shortcut.log            text log
#   scenario_summary.csv                one row per scenario x rule/variant
#   A_by_N_all.csv                      A(N) per scenario and rule (15 sig. digits)
#   wilson_A_by_N.csv                   same, only the Wilson et al. scenarios
#   arm_size_table.csv                  per scenario x arm x n: contiguity, P, masses by cause
#   noncontiguous_arm_sizes.csv         (scenario, arm, n) with non-contiguous failure set
#   package_crosscheck.csv              enumeration vs bam_sample_size(method = "exact")
#   unimodality_stress.csv              extra: is W(x) quasi-concave (failure set contiguous
#                                       for every w) over a grid of Beta analysis priors?
# =============================================================================

TOL <- 1e-12

# ------------------------------------------------------------------ libs ----
.setup_libs <- function() {
  pk <- Sys.getenv("DTA_PKG_LIB", unset = "")
  ex <- Sys.getenv("DTA_EXTRA_LIB", unset = "")
  paths <- c(pk, ex)
  paths <- paths[nzchar(paths)]
  if (length(paths)) .libPaths(c(paths, .libPaths()))
}

# ------------------------------------------------- arm-level computations ----
# Beta-Binomial pmf of the count x = 0..n under the DESIGN prior (same formula
# as .bam_exact_width_prob() in the package).
.bb_pmf <- function(n, a_d, b_d) {
  x <- 0:n
  exp(lchoose(n, x) + lbeta(a_d + x, b_d + n - x) - lbeta(a_d, b_d))
}

# Full width of the equal-tailed ANALYSIS posterior for x = 0..n.
.width_vec <- function(n, a_a, b_a, level) {
  x <- 0:n
  lq <- (1 - level) / 2
  stats::qbeta(1 - lq, a_a + x, b_a + n - x) - stats::qbeta(lq, a_a + x, b_a + n - x)
}

# TRUE iff W(0..n) is quasi-concave (no strict decrease followed by a strict
# increase; steps with |diff| <= tol are ignored). Then the failure set
# {W > w} is contiguous for EVERY threshold w, i.e. the printed rule equals the
# enumeration whatever the target width.
.is_quasi_concave <- function(W, tol = TOL) {
  d <- diff(W)
  sg <- sign(d[abs(d) > tol])
  !any(diff(sg) > 0)
}

# Success set of the CODE rule (emulation of Data S2 cutoffs()/assurance()).
# `fail` is a logical vector over x = 0..n (W(x) > w + TOL). n >= 1.
.code_rule_sets <- function(fail, n) {
  idx <- seq_len(n)
  fi <- idx[fail[idx + 1L]]
  c1 <- if (length(fi)) fi[1L] - 1L else NA_integer_
  sq <- 1:(n - 1)                       # R semantics: n = 1 gives c(1, 0)
  sq <- sq[sq >= 0]                     # (index fail[sq + 1] must exist)
  fi2 <- sq[fail[sq + 1L]]
  c2 <- if (length(fi2)) fi2[length(fi2)] + 1L else n
  s <- rep(FALSE, n + 1L)
  if (is.na(c1)) {
    s[] <- TRUE
  } else {
    if (c1 > 1L) s[1:(c1 + 1L)] <- TRUE
    if (c2 < n) s[(min(c2, n) + 1L):(n + 1L)] <- TRUE
  }
  list(s = s, c1 = c1, c2 = c2)
}

# Everything about one arm size n: probabilities under the three rules, the
# contiguity of the failure set and the mass by deviation cause.
arm_size_stats <- function(n, a_d, b_d, a_a, b_a, w, level, tol = TOL) {
  if (n == 0L) {
    return(c(n = 0, nfail = NA, F_lo = NA, F_hi = NA, n_gap = NA, contiguous = NA, unimodal_W = NA,
             ties = 0, c1 = NA, c2 = NA,
             P_enum = 0, P_print = 0, P_code = 0,
             m_b = 0, m_c = 0, m_e = 0, m_d = 0, m_f = 0, m_g = 0,
             n_b = 0, n_c = 0, n_e = 0, n_d = 0, n_f = 0, n_g = 0))
  }
  W <- .width_vec(n, a_a, b_a, level)
  pm <- .bb_pmf(n, a_d, b_d)
  x <- 0:n
  fail <- W > w + tol
  enum_s <- !fail
  # printed rule: complement of the hull of F
  if (any(fail)) {
    wf <- which(fail) - 1L
    F_lo <- min(wf); F_hi <- max(wf)
    print_s <- x < F_lo | x > F_hi
    n_gap <- (F_hi - F_lo + 1L) - length(wf)
    nfail <- length(wf)
  } else {
    F_lo <- NA_integer_; F_hi <- NA_integer_
    print_s <- rep(TRUE, n + 1L)
    n_gap <- 0L; nfail <- 0L
  }
  cr <- .code_rule_sets(fail, n)
  code_s <- cr$s
  lost <- enum_s & !code_s
  spur <- !enum_s & code_s
  lb <- lost & !is.na(cr$c1) & cr$c1 <= 1L & x <= cr$c1
  lc <- lost & !lb & x == n & cr$c2 >= n
  le <- lost & !lb & !lc
  sd <- spur & x == 0L
  sf <- spur & !sd & x == n
  sg <- spur & !sd & !sf
  c(n = n, nfail = nfail, F_lo = F_lo, F_hi = F_hi, n_gap = n_gap,
    contiguous = as.numeric(n_gap == 0L), unimodal_W = as.numeric(.is_quasi_concave(W, tol)),
    ties = sum(abs(W - w) <= tol),
    c1 = cr$c1, c2 = cr$c2,
    P_enum = sum(pm[enum_s]), P_print = sum(pm[print_s]), P_code = sum(pm[code_s]),
    m_b = sum(pm[lb]), m_c = sum(pm[lc]), m_e = sum(pm[le]),
    m_d = sum(pm[sd]), m_f = sum(pm[sf]), m_g = sum(pm[sg]),
    n_b = sum(lb), n_c = sum(lc), n_e = sum(le), n_d = sum(sd), n_f = sum(sf), n_g = sum(sg))
}

# Table for arm sizes n_from..n_to.
arm_table <- function(n_from, n_to, a_d, b_d, a_a, b_a, w, level, tol = TOL) {
  m <- t(vapply(n_from:n_to, arm_size_stats, numeric(25),
                a_d = a_d, b_d = b_d, a_a = a_a, b_a = b_a, w = w, level = level, tol = tol))
  as.data.frame(m)
}

# ---------------------------------------------- assurance from arm tables ----
# Beta-Binomial pmf of the number of diseased n_d in N.
.prev_pmf <- function(N, a_p, b_p) {
  k <- 0:N
  exp(lchoose(N, k) + lbeta(a_p + k, b_p + N - k) - lbeta(a_p, b_p))
}

# A(N) for N in Ns. P1/P2 are vectors of P(success | arm size n) indexed n + 1
# (one per arm); P2 = NULL for the marginal (one-arm) objective.
assurance_vec <- function(Ns, a_p, b_p, Pse, Psp = NULL) {
  vapply(Ns, function(N) {
    pk <- .prev_pmf(N, a_p, b_p)
    k <- 0:N
    if (is.null(Psp)) sum(pk * Pse[k + 1L]) else sum(pk * Pse[k + 1L] * Psp[N - k + 1L])
  }, numeric(1))
}

first_cross <- function(Ns, A, target) {
  i <- which(A >= target)
  if (length(i)) Ns[i[1L]] else NA_integer_
}

# --------------------------------------------- reusable stand-alone API -----
emulate_data_s2_arm_prob <- function(n, design_measure, analysis_measure, w,
                                     level = 0.95, tol = TOL,
                                     rule = c("code", "printed", "enumeration")) {
  rule <- match.arg(rule)
  col <- switch(rule, code = "P_code", printed = "P_print", enumeration = "P_enum")
  vapply(n, function(nn) {
    arm_size_stats(as.integer(nn), design_measure[1], design_measure[2],
                   analysis_measure[1], analysis_measure[2], w, level, tol)[[col]]
  }, numeric(1))
}

emulate_data_s2_assurance <- function(N, w, design_measure, analysis_measure,
                                      design_prev, level = 0.95, tol = TOL,
                                      rule = c("code", "printed", "enumeration")) {
  rule <- match.arg(rule)
  N <- as.integer(N)
  P <- emulate_data_s2_arm_prob(0:max(N), design_measure, analysis_measure, w,
                                level, tol, rule)
  assurance_vec(N, design_prev[1], design_prev[2], P)
}

# ------------------------------------------------------------- scenarios ----
make_scenarios <- function(grid_path) {
  g <- jsonlite::fromJSON(grid_path, simplifyVector = TRUE)$scenarios
  sc <- list()
  add <- function(...) { s <- list(...); sc[[s$id]] <<- s }
  add(id = "P1", family = "joint", des_se = c(17, 3), des_sp = c(2, 2), ana_se = c(17, 3),
      ana_sp = c(2, 2), prev = c(4, 16), w_se = 0.14, w_sp = 0.10, target = 0.80,
      level = 0.95, cap0 = 500, N_hi = 1200, group = "article",
      published = "N = 678 (0.6.9 article, abstract)")
  add(id = "P2", family = "joint", des_se = c(17, 3), des_sp = c(18, 2), ana_se = c(17, 3),
      ana_sp = c(18, 2), prev = c(4, 16), w_se = 0.14, w_sp = 0.10, target = 0.80,
      level = 0.95, cap0 = 500, N_hi = 1200, group = "article",
      published = "N = 672 (0.6.9 article, Table 2)")
  for (i in seq_len(nrow(g))) {
    s <- g[i, ]
    add(id = s$id, family = "joint", des_se = as.numeric(s$prior_se[[1]]),
        des_sp = as.numeric(s$prior_sp[[1]]), ana_se = as.numeric(s$prior_se[[1]]),
        ana_sp = as.numeric(s$prior_sp[[1]]), prev = as.numeric(s$prior_prev[[1]]),
        w_se = s$delta_se, w_sp = s$delta_sp, target = s$target, level = s$level,
        cap0 = min(400L, s$N_hi), N_hi = s$N_hi, group = "grid",
        published = "LOCKED_confirmatory_grid.json (no published N)")
  }
  add(id = "WILSON_VAP", family = "se", des_se = c(25.9, 2.1), ana_se = c(25.9, 2.1),
      prev = c(29, 98), w_se = 0.16, target = 0.80, level = 0.95, cap0 = 320, N_hi = 320,
      group = "wilson",
      published = "n* = 106, A(150) = 0.88 (Wilson 2022, section 5.1)")
  add(id = "WILSON_B4_flat", family = "se", des_se = c(25.9, 2.1), ana_se = c(1, 1),
      prev = c(29, 98), w_se = 0.16, target = 0.80, level = 0.95, cap0 = 400, N_hi = 400,
      group = "wilson",
      published = "290 (Combination / Hybrid II; Data S1 B.4 Table 1, L268-292); Standard/Hybrid I = WILSON_VAP = 106")
  add(id = "WILSON_61_rho010", family = "se", des_se = c(4.5, 0.5), ana_se = c(4.5, 0.5),
      prev = c(5, 45), w_se = 0.18, target = 0.80, level = 0.95, cap0 = 800, N_hi = 800,
      group = "wilson",
      published = "681 (Wilson 2022, section 6.1, lambda = 0.9, n0 = 50, rho = 0.10)")
  add(id = "WILSON_61_rho005", family = "se", des_se = c(2.25, 0.25), ana_se = c(2.25, 0.25),
      prev = c(2.5, 47.5), w_se = 0.18, target = 0.80, level = 0.95, cap0 = 1900, N_hi = 1900,
      group = "wilson",
      published = "1643 (Wilson 2022, section 6.1, rho = 0.05); additional case")
  sc
}

VARIANTS <- c("enumeration", "printed", "code", "only_b", "only_c", "only_e", "only_d", "only_f")

# Per-n probability vectors (indexed n + 1) for each variant from the arm table.
variant_P <- function(tb) {
  P <- list(
    enumeration = tb$P_enum,
    printed = tb$P_print,
    code = tb$P_code,
    only_b = tb$P_enum - tb$m_b,
    only_c = tb$P_enum - tb$m_c,
    only_e = tb$P_enum - tb$m_e,
    only_d = tb$P_enum + tb$m_d,
    only_f = tb$P_enum + tb$m_f
  )
  P
}

# Run one scenario: build arm tables in blocks until the three main rules
# have crossed and the table covers 1.1 * N0 + 10 (or N_hi).
run_scenario <- function(s) {
  cap <- min(s$cap0, s$N_hi)
  tb_se <- NULL; tb_sp <- NULL; built <- -1L
  repeat {
    new <- (built + 1L):cap
    t1 <- arm_table(min(new), max(new), s$des_se[1], s$des_se[2], s$ana_se[1], s$ana_se[2],
                    s$w_se, s$level)
    tb_se <- rbind(tb_se, t1)
    if (s$family == "joint") {
      t2 <- arm_table(min(new), max(new), s$des_sp[1], s$des_sp[2], s$ana_sp[1], s$ana_sp[2],
                      s$w_sp, s$level)
      tb_sp <- rbind(tb_sp, t2)
    }
    built <- cap
    Ns <- 2:cap
    Pse <- variant_P(tb_se)
    Psp <- if (s$family == "joint") variant_P(tb_sp) else NULL
    A <- lapply(VARIANTS, function(v)
      assurance_vec(Ns, s$prev[1], s$prev[2], Pse[[v]], if (is.null(Psp)) NULL else Psp[[v]]))
    names(A) <- VARIANTS
    Nst <- vapply(A, function(a) first_cross(Ns, a, s$target), numeric(1))
    main <- Nst[c("enumeration", "printed", "code")]
    if (all(!is.na(main))) {
      need <- ceiling(1.1 * max(main)) + 10L
      if (cap >= min(need, s$N_hi)) break
      cap <- min(need, s$N_hi)
    } else {
      if (cap >= s$N_hi) break
      cap <- min(cap + 400L, s$N_hi)
    }
  }
  list(s = s, Ns = Ns, A = A, Nstar = Nst, tb_se = tb_se, tb_sp = tb_sp, cap = cap)
}

# ------------------------------------------------------------------ main ----
main <- function() {
  .setup_libs()
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1)
  out <- args[1]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(out, "success_set_shortcut.log")
  zz <- file(log_file, open = "wt", encoding = "UTF-8")
  sink(zz, split = TRUE)
  on.exit({ sink(); close(zz) }, add = TRUE)
  t_start <- Sys.time()
  options(width = 170)
  cat("success_set_shortcut.R\n")
  cat("Date:", format(t_start, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("R:", R.version.string, "\n")
  suppressMessages(library(dtasamplesize))
  cat("packageVersion(dtasamplesize):", as.character(utils::packageVersion("dtasamplesize")), "\n")
  cat("find.package(dtasamplesize):", find.package("dtasamplesize"), "\n")
  cat("TOL =", format(TOL), "\n\n")

  script_dir <- {
    a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
    if (length(f)) dirname(normalizePath(f[1], winslash = "/")) else getwd()
  }
  grid_path <- file.path(script_dir, "python", "LOCKED_confirmatory_grid.json")
  scen <- make_scenarios(grid_path)

  res <- list()
  summ <- list(); allA <- list(); armtab <- list(); noncont <- list()
  for (id in names(scen)) {
    t0 <- Sys.time()
    r <- run_scenario(scen[[id]])
    res[[id]] <- r
    s <- r$s
    cat(sprintf("---- %s [%s] ; objective %s ; target %.2f ; level %.2f ; table up to n = %d (%.1f s)\n",
                id, s$group, if (s$family == "joint") "JOINT" else "marginal Se", s$target,
                s$level, r$cap, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    cat("  published:", s$published, "\n")
    Nst <- r$Nstar
    Ns <- r$Ns
    for (v in VARIANTS) {
      a <- r$A[[v]]
      nst <- Nst[[v]]
      a_at <- if (is.na(nst)) NA_real_ else a[Ns == nst]
      post <- if (is.na(nst)) NA_integer_ else sum(a[Ns > nst] < s$target)
      a_at_enum <- if (is.na(Nst[["enumeration"]])) NA_real_ else a[Ns == Nst[["enumeration"]]]
      a_at_code <- if (is.na(Nst[["code"]])) NA_real_ else a[Ns == Nst[["code"]]]
      summ[[length(summ) + 1L]] <- data.frame(
        scenario = id, group = s$group, objective = if (s$family == "joint") "joint" else "se_marginal",
        rule = v, target = s$target, level = s$level, N_star = nst, A_at_N_star = a_at,
        A_at_N_star_enumeration = a_at_enum, A_at_N_star_code = a_at_code,
        n_N_above_star_with_A_below_target = post, max_A_evaluated = max(a),
        N_at_max_A = Ns[which.max(a)], N_grid_max = max(Ns), published = s$published,
        stringsAsFactors = FALSE)
      allA[[length(allA) + 1L]] <- NULL
    }
    # A(N) table for this scenario
    d <- data.frame(scenario = id, group = s$group, N = Ns,
                    A_code = r$A$code, A_printed = r$A$printed, A_enumeration = r$A$enumeration,
                    stringsAsFactors = FALSE)
    allA[[id]] <- d
    cat(sprintf("  N*: enumeration %s | printed %s | code %s | only_b %s | only_c %s | only_e %s | only_d %s | only_f %s\n",
                Nst[["enumeration"]], Nst[["printed"]], Nst[["code"]], Nst[["only_b"]],
                Nst[["only_c"]], Nst[["only_e"]], Nst[["only_d"]], Nst[["only_f"]]))
    for (v in c("enumeration", "printed", "code")) {
      nst <- Nst[[v]]
      if (!is.na(nst)) cat(sprintf("  A(N*) %-11s = %.10f at N* = %d ; A(N*-1) = %.10f\n", v,
                                   r$A[[v]][Ns == nst], nst, if (nst > 2) r$A[[v]][Ns == nst - 1] else NA))
      else cat(sprintf("  %-11s: no crossing up to N = %d (max A = %.6f at N = %d)\n", v, max(Ns),
                       max(r$A[[v]]), Ns[which.max(r$A[[v]])]))
    }
    dif_code <- r$A$code - r$A$enumeration
    dif_print <- r$A$printed - r$A$enumeration
    cat(sprintf("  max |A_code - A_enum| = %.3e (N = %d) ; max |A_printed - A_enum| = %.3e (N = %d)\n",
                max(abs(dif_code)), Ns[which.max(abs(dif_code))],
                max(abs(dif_print)), Ns[which.max(abs(dif_print))]))
    # contiguity, per arm
    arms <- list(Se = r$tb_se)
    if (s$family == "joint") arms$Sp <- r$tb_sp
    nref <- Nst[["enumeration"]]
    for (an in names(arms)) {
      tb <- arms[[an]]
      tb$scenario <- id; tb$arm <- an
      armtab[[paste(id, an)]] <- tb
      nc <- tb[!is.na(tb$contiguous) & tb$contiguous == 0, ]
      if (nrow(nc)) noncont[[paste(id, an)]] <- nc[, c("scenario", "arm", "n", "nfail", "F_lo", "F_hi", "n_gap", "ties", "P_enum", "P_print", "P_code")]
      tt <- tb[tb$n >= 1 & (is.na(nref) | tb$n <= nref), ]
      # relevant support at N*_enum: arm sizes with pmf > 1e-8
      if (!is.na(nref)) {
        pk <- .prev_pmf(nref, s$prev[1], s$prev[2]); k <- 0:nref
        rel <- if (an == "Se") k[pk > 1e-8] else (nref - k)[pk > 1e-8]
        rel <- rel[rel >= 1]
        tr <- tb[tb$n %in% rel, ]
        cat(sprintf("  arm %s: arm sizes n in 1..%d: %d non-contiguous; relevant support at N*_enum (pmf>1e-8, n = %d..%d): %d non-contiguous of %d; ties(|W-w*|<=TOL) over n<=%d: %d\n",
                    an, nref, sum(tt$contiguous == 0), min(rel), max(rel),
                    sum(tr$contiguous == 0), nrow(tr), nref, sum(tt$ties)))
        cat(sprintf("         relevant n with P_code != P_enum (>1e-15): %d ; with P_print != P_enum: %d ; max |P_code - P_enum| = %.3e ; max |P_print - P_enum| = %.3e\n",
                    sum(abs(tr$P_code - tr$P_enum) > 1e-15), sum(abs(tr$P_print - tr$P_enum) > 1e-15),
                    max(abs(tr$P_code - tr$P_enum)), max(abs(tr$P_print - tr$P_enum))))
        cat(sprintf("         cause counts over relevant n (n_b, n_c, n_e, n_d, n_f, n_g as number of counts x): b=%d c=%d e=%d d=%d f=%d g=%d ; masses: b=%.3e c=%.3e e=%.3e d=%.3e f=%.3e g=%.3e\n",
                    sum(tr$n_b), sum(tr$n_c), sum(tr$n_e), sum(tr$n_d), sum(tr$n_f), sum(tr$n_g),
                    max(tr$m_b), max(tr$m_c), max(tr$m_e), max(tr$m_d), max(tr$m_f), max(tr$m_g)))
      }
      # consistency: code probability = enumeration - lost + spurious
      comp <- tb$P_enum - tb$m_b - tb$m_c - tb$m_e + tb$m_d + tb$m_f + tb$m_g
      cat(sprintf("  arm %s: max |P_code - (P_enum - lost + spurious)| = %.2e over n = %d..%d\n",
                  an, max(abs(tb$P_code - comp)), min(tb$n), max(tb$n)))
      cat(sprintf("  arm %s: arm sizes n = 1..%d where W(x) is NOT quasi-concave in x (failure set non-contiguous for some w): %d\n",
                  an, max(tb$n), sum(tb$unimodal_W == 0, na.rm = TRUE)))
    }
    cat("\n")
  }

  # ---------------------------------------------- package cross-check -------
  cat("==== Cross-check of the enumeration against bam_sample_size(method = 'exact') ====\n")
  options(dtasamplesize.warn_small_B = FALSE)
  xc <- list()
  for (id in names(scen)) {
    s <- scen[[id]]
    if (id == "WILSON_61_rho005") next        # extra scenario, O(N^2) package cost: skipped
    r <- res[[id]]
    cl <- list(prior_se = s$des_se, prior_prev = s$prev, delta_se = s$w_se,
               target_assurance = s$target, alpha_ci = s$level, N_max = as.integer(s$N_hi),
               B = 1000, n_range = 20:60)
    if (s$family == "joint") {
      cl$prior_sp <- s$des_sp; cl$delta_sp <- s$w_sp
    } else {
      cl$assurance_for <- "se"; cl$prior_sp <- c(2, 2); cl$delta_sp <- 0.1
    }
    if (!identical(s$ana_se, s$des_se)) cl$analysis_prior_se <- s$ana_se
    pk <- suppressWarnings(do.call(bam_sample_size, cl))
    nst <- r$Nstar[["enumeration"]]
    a_mine <- if (is.na(nst)) max(r$A$enumeration) else r$A$enumeration[r$Ns == nst]
    a_pk <- if (isTRUE(pk$target_reached)) pk$joint_assurance else pk$max_assurance_evaluated
    ok <- identical(as.integer(pk$N_total), as.integer(nst)) ||
      (is.na(pk$N_total) && is.na(nst))
    cat(sprintf("  %-18s package N_total = %s ; enumeration N* = %s ; |A diff| = %.2e ; %s\n", id,
                format(pk$N_total), format(nst), abs(a_mine - a_pk), if (ok) "MATCH" else "DIFFERENT"))
    xc[[id]] <- data.frame(scenario = id, package_N_total = pk$N_total, enumeration_N_star = nst,
                           package_A = a_pk, enumeration_A = a_mine, abs_diff_A = abs(a_mine - a_pk),
                           match_N = ok, stringsAsFactors = FALSE)
  }
  xc <- do.call(rbind, xc)

  # ---------------------------------- stress: is W(x) always quasi-concave? ---
  # (extra, outside the scenarios): if W(x) is quasi-concave in x, the failure
  # set {W > w} is contiguous for every w, so the printed rule equals the
  # enumeration. Sweep over Beta analysis priors and arm sizes n = 1..300.
  cat("\n==== Stress: quasi-concavity of W(x) over a grid of analysis priors (extra) ====\n")
  hp <- c(0.1, 0.5, 1, 2, 5, 20, 50)
  st <- list()
  for (lev in c(0.90, 0.95, 0.99)) {
    for (a in hp) for (b in hp) {
      bad <- vapply(1:300, function(n) !.is_quasi_concave(.width_vec(n, a, b, lev)), logical(1))
      st[[length(st) + 1L]] <- data.frame(level = lev, a = a, b = b, n_max = 300L,
                                          n_not_quasi_concave = sum(bad),
                                          first_n = if (any(bad)) which(bad)[1L] else NA_integer_)
    }
  }
  st <- do.call(rbind, st)
  utils::write.csv(st, file.path(out, "unimodality_stress.csv"), row.names = FALSE)
  cat(sprintf("  prior grid a, b in {%s}; levels 0.90, 0.95, 0.99; n = 1..300: %d (prior, level) combinations; combinations with some non-quasi-concave n: %d; total non-quasi-concave (prior, level, n): %d\n",
              paste(hp, collapse = ", "), nrow(st), sum(st$n_not_quasi_concave > 0), sum(st$n_not_quasi_concave)))
  if (any(st$n_not_quasi_concave > 0)) print(st[st$n_not_quasi_concave > 0, ], row.names = FALSE)

  # ------------------------------------------------------------ write -------
  summ <- do.call(rbind, summ)
  allA <- do.call(rbind, allA)
  armtab <- do.call(rbind, armtab)
  noncont <- if (length(noncont)) do.call(rbind, noncont) else
    data.frame(scenario = character(), arm = character(), n = numeric(), nfail = numeric(),
               F_lo = numeric(), F_hi = numeric(), n_gap = numeric(), ties = numeric(),
               P_enum = numeric(), P_print = numeric(), P_code = numeric())
  fmt15 <- function(d) {
    for (nm in c("A_code", "A_printed", "A_enumeration")) d[[nm]] <- sprintf("%.15g", d[[nm]])
    d
  }
  utils::write.csv(summ, file.path(out, "scenario_summary.csv"), row.names = FALSE)
  utils::write.csv(fmt15(allA), file.path(out, "A_by_N_all.csv"), row.names = FALSE)
  utils::write.csv(fmt15(allA[allA$group == "wilson", ]), file.path(out, "wilson_A_by_N.csv"), row.names = FALSE)
  utils::write.csv(armtab, file.path(out, "arm_size_table.csv"), row.names = FALSE)
  utils::write.csv(noncont, file.path(out, "noncontiguous_arm_sizes.csv"), row.names = FALSE)
  utils::write.csv(xc, file.path(out, "package_crosscheck.csv"), row.names = FALSE)

  # ----------------------------------------------------- final table --------
  cat("\n==== SUMMARY: N* by rule (first N >= 2 with A >= target) ====\n")
  w <- reshape(summ[, c("scenario", "rule", "N_star")], idvar = "scenario", timevar = "rule",
               direction = "wide")
  names(w) <- sub("N_star\\.", "", names(w))
  print(w[match(names(scen), w$scenario), ], row.names = FALSE)
  cat("\nA(N*) of each main rule at its own N*, and the code/printed A at the enumeration N*:\n")
  m <- summ[summ$rule %in% c("enumeration", "printed", "code"),
            c("scenario", "rule", "N_star", "A_at_N_star", "A_at_N_star_enumeration")]
  print(m, row.names = FALSE, digits = 8)
  cat("\nPublished figures (Wilson et al.) against the rules:\n")
  g <- function(id, v) summ$N_star[summ$scenario == id & summ$rule == v]
  cat(sprintf("  VAP: published 106 | code %s | printed %s | enumeration %s\n", g("WILSON_VAP", "code"),
              g("WILSON_VAP", "printed"), g("WILSON_VAP", "enumeration")))
  va <- allA[allA$scenario == "WILSON_VAP", ]
  for (nn in c(104, 106, 150)) {
    cat(sprintf("       A(%d): code %.4f | printed %.4f | enumeration %.4f\n", nn,
                va$A_code[va$N == nn], va$A_printed[va$N == nn], va$A_enumeration[va$N == nn]))
  }
  cat(sprintf("  B.4 flat analysis prior: published 290 | code %s | printed %s | enumeration %s\n",
              g("WILSON_B4_flat", "code"), g("WILSON_B4_flat", "printed"), g("WILSON_B4_flat", "enumeration")))
  cat(sprintf("  6.1 rho = 0.10: published 681 | code %s | printed %s | enumeration %s\n",
              g("WILSON_61_rho010", "code"), g("WILSON_61_rho010", "printed"), g("WILSON_61_rho010", "enumeration")))
  cat(sprintf("  6.1 rho = 0.05 (extra): published 1643 | code %s | printed %s | enumeration %s\n",
              g("WILSON_61_rho005", "code"), g("WILSON_61_rho005", "printed"), g("WILSON_61_rho005", "enumeration")))
  cat(sprintf("\nNon-contiguous (scenario, arm, n) rows written: %d\n", nrow(noncont)))
  cat(sprintf("Total elapsed: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  cat("Output directory:", normalizePath(out, winslash = "/"), "\n")
  invisible(NULL)
}

if (!isTRUE(getOption("dta.success_set.no_main"))) main()
