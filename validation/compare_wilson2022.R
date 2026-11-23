#!/usr/bin/env Rscript
# =============================================================================
# compare_wilson2022.R
#
# Direct comparison of dtasamplesize (exact enumeration, bam_* functions) with
# the worked examples of Wilson KJ, Williamson SF, Allen AJ, Williams CJ,
# Hellyer TP, Lendrem BC. Stat Med. 2022;41(15):2908-2922. doi:10.1002/sim.9393
# [W]; its Supporting Information Data S1 [WS1].
#
# Purpose
#   For each example of [W] (all inputs stated), compute the quantity with the
#   package (exact enumeration of the Beta-Binomial outcomes) and with an
#   EMULATION of the success-set rule as it is programmed in the authors'
#   code (as described in a written description of that rule kept by the
#   authors; the emulation below was written from that description, not from
#   the code text of Data S2), and tabulate: published value (with source), package,
#   emulation, [optional] the original Data S2 code run in an isolated folder,
#   difference and explanation (when one is documented).
#   No code of Data S2 is read, copied, adapted or sourced by this script.
#
# Usage
#   Rscript compare_wilson2022.R <output_dir> [<dataS2_results.csv>]
#   The script writes ONLY into <output_dir> (created if missing). The optional
#   second argument is the CSV (columns key,value) written by the separate
#   isolated run of the original code (run_dataS2.R, outside this repository);
#   when given, its values are merged into the table by key.
#   Libraries: if the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB
#   exist they are prepended to .libPaths() (in that order).
#   Optional: WILSON_SECTIONS="a,b,c,d,e" restricts the sections (development).
#
# Conventions (named explicitly)
#   * width: TOTAL width of the equal-tailed two-sided 95% posterior interval
#     (full_width), as in [W] sections 3.1.1 and 5.1. Comparisons
#     "width <= target" use TOL = 1e-12 (width <= target + TOL), the same in
#     all validation scripts. (The original code compares width > w without
#     tolerance; they can differ only in an exact tie.)
#   * objective: MARGINAL assurance of Se (assurance_for = "se"); total N is
#     the sample size (n_T of [W]); the number of diseased is Beta-Binomial
#     from the prevalence prior.
#   * design prior vs analysis prior: the design prior generates the data and
#     weights the outcomes, the analysis prior fixes the interval widths.
#     Where the example does not state an analysis prior it equals the design
#     prior ([W] section 5.1; for section 6.1 this is an inference).
#     The prevalence has no analysis prior.
#   * empty arm: no diseased participants counts as a failure in the package
#     (P(0) = 0, documented in ?bam_sample_size); the code rule omits n_T,1 = 0
#     as well (deviation (a) of hechos 2.15).
#   * search rule: package, first N with A(N) >= 0.8 (integer scan from N = 2);
#     emulation, first N with A(N) > 0.8 (strict, as the original search), scan
#     from N = 1.
#   * emulation of the code rule, per arm of size m (analysis prior fixes the
#     widths W(i), i = 1..m, fail if W(i) > w):
#       c1 = (first i in 1..m that fails) - 1, NA if none;
#       c2 = (last i in 1..m-1 that fails) + 1, m if none (for m = 1 the index
#            sequence 1:(m-1) of R is c(1, 0));
#       success set: all 0..m if c1 is NA; 0..c1 and c2..m (only if c2 < m)
#       if c1 > 1; c2..m (only if c2 < m) if c1 in {0, 1}; empty otherwise.
#     n11 = 0 is accepted without evaluating W(0); n_T,1 = 0 is omitted.
#   * exact rule ("package"): every n11 in 0..m with W(n11) <= w is a success.
#   * distance of the sensitivity analysis: "H2" = 1 - BC (the expression
#     printed by [W] under the name Hellinger distance and used by their
#     code), epsilon = 0.00354; see ?bam_prior_contour.
#
# Sources of the published figures
#   [W]   sections 5.1 (106, 88%), 5.2 (Table 2), 5.3 (99th and 76th
#         percentiles), 6.1 (681, 1643, 2770), as listed in
#         sections 5-7 of [W].
#   [WS1] Data S1: B.4 Table 1 (106, 290, 106, 290; posterior means);
#         C Tables 2-4 (w* = 0.22, 0.18, 0.14; lines L310-442 of Data_S1.txt).
#
# Outputs (in <output_dir>)
#   compare_wilson2022.log           log of the run
#   compare_wilson2022_table.csv     one row per quantity (see columns there)
#   table2_grid_points.csv           the 63 grid points of the rule of Table 2 of [W] (file name keeps "table2")
#   suppl_C_cells.csv                Supplement C, cell by cell
#   section6_1_curves.csv            A(N) at selected N, section 6.1
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript compare_wilson2022.R <output_dir> [<dataS2_results.csv>]")
outdir <- args[1L]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ds2_csv <- if (length(args) >= 2L) args[2L] else NA_character_

libs <- c(Sys.getenv("DTA_PKG_LIB"), Sys.getenv("DTA_EXTRA_LIB"))
libs <- libs[nzchar(libs)]
if (length(libs)) .libPaths(c(libs, .libPaths()))

logfile <- file.path(outdir, "compare_wilson2022.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)

suppressPackageStartupMessages(library(dtasamplesize))
T_START <- Sys.time()
cat("compare_wilson2022.R\n")
cat("date:            ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version:       ", R.version.string, "\n")
cat("dtasamplesize:   ", as.character(packageVersion("dtasamplesize")), "\n")
cat("find.package:    ", find.package("dtasamplesize"), "\n")
cat("output dir:      ", normalizePath(outdir, winslash = "/"), "\n")
cat("Data S2 results: ", if (is.na(ds2_csv)) "(none)" else ds2_csv, "\n")

TOL     <- 1e-12
TARGET  <- 0.8
EPS     <- 0.00354
SECT    <- Sys.getenv("WILSON_SECTIONS", "a,b,c,d,e,f")
SECT    <- strsplit(SECT, ",")[[1L]]
want    <- function(s) s %in% SECT
cat("TOL =", format(TOL), "; target assurance =", TARGET, "; sections:", SECT, "\n\n")

# ---- helpers: Beta-Binomial, widths -----------------------------------------
bbpmf <- function(k, n, a, b) exp(lchoose(n, k) + lbeta(k + a, n - k + b) - lbeta(a, b))
w_full <- function(m, a, b, lev = 0.95) {
  i <- 0:m; al <- (1 - lev) / 2
  qbeta(1 - al, a + i, b + m - i) - qbeta(al, a + i, b + m - i)
}

# ---- emulation of the code rule and exact rule, per arm ---------------------
# fail flags (index i + 1, i = 0..m) for analysis prior an and width w
fail_flags <- function(m, an, w) w_full(m, an[1L], an[2L]) > w + TOL

# success probability of an arm of size m under design prior dsg, given flags
arm_P_from_flags <- function(m, fail, dsg) {
  pm <- bbpmf(0:m, m, dsg[1L], dsg[2L])
  pe <- sum(pm[!fail])
  f1 <- fail[-1L]                                   # i = 1..m
  c1 <- if (any(f1)) which(f1)[1L] - 1L else NA_integer_
  if (m >= 2L) {
    f2 <- f1[seq_len(m - 1L)]                        # i = 1..m-1
    c2 <- if (any(f2)) max(which(f2)) + 1L else m
  } else {                                           # m = 1: 1:(m-1) is c(1, 0)
    c2 <- if (fail[1L]) 1L else if (fail[2L]) 2L else 1L
  }
  idx <- if (is.na(c1)) 0:m
         else if (c1 > 1L) c(0:c1, if (c2 < m) c2:m)
         else if (c2 < m) c2:m else integer(0)
  c(code = sum(pm[idx + 1L]), exact = pe)
}

# incremental per-arm table: P_code[m], P_exact[m], m = 1..M (m = 0 -> 0)
new_arms <- function(an, dsg, w, flags_cache = NULL) {
  e <- new.env(); e$an <- an; e$dsg <- dsg; e$w <- w; e$M <- 0L
  e$Pc <- numeric(0); e$Pe <- numeric(0); e$cache <- flags_cache
  e
}
extend_arms <- function(e, M) {
  if (M <= e$M) return(invisible(e))
  for (m in (e$M + 1L):M) {
    fl <- NULL
    if (!is.null(e$cache)) {
      key <- as.character(m)
      fl <- e$cache[[key]]
      if (is.null(fl)) { fl <- fail_flags(m, e$an, e$w); e$cache[[key]] <- fl }
    } else fl <- fail_flags(m, e$an, e$w)
    p <- arm_P_from_flags(m, fl, e$dsg)
    e$Pc[m] <- p[["code"]]; e$Pe[m] <- p[["exact"]]
  }
  e$M <- M
  invisible(e)
}
# assurance at (possibly non-integer) N; rule "code" or "exact"
A_arms <- function(e, N, prev, rule = c("code", "exact")) {
  rule <- match.arg(rule)
  K <- floor(N)
  if (K < 1) return(0)
  extend_arms(e, K)
  P <- if (rule == "code") e$Pc[1:K] else e$Pe[1:K]
  sum(bbpmf(1:K, N, prev[1L], prev[2L]) * P)
}
# first N (scan from N0) with A > target (strict) or A >= target
first_cross <- function(e, prev, rule, Nmax, N0 = 1L, strict = TRUE, step = 50L) {
  last <- NA_real_
  for (N in N0:Nmax) {
    if (N > e$M) extend_arms(e, min(Nmax, N + step - 1L))
    a <- A_arms(e, N, prev, rule)
    if (if (strict) a > TARGET else a >= TARGET - 0) return(c(N = N, A = a))
  }
  c(N = NA_real_, A = NA_real_)
}

# ---- helpers for the package -----------------------------------------------
pk_nstar <- function(prior_se, prior_prev, delta, an = NULL, N_max = 3000L) {
  r <- suppressWarnings(bam_sample_size(
    prior_se = prior_se, prior_sp = c(2, 2), delta_se = delta, delta_sp = delta,
    prior_prev = prior_prev, assurance_for = "se", method = "exact",
    analysis_prior_se = an, N_max = N_max, B = 1000, n_range = 20))
  c(N = if (isTRUE(r$target_reached)) as.numeric(r$N_total) else NA_real_,
    A = if (isTRUE(r$target_reached)) r$joint_assurance else NA_real_,
    Amax = r$max_assurance_evaluated, Nmax = r$N_at_max_assurance)
}
pk_A <- function(prior_se, prior_prev, delta, N, an = NULL) {
  cv <- suppressWarnings(bam_assurance_curve(
    prior_se = prior_se, prior_sp = c(2, 2), delta_se = delta, delta_sp = delta,
    prior_prev = prior_prev, N_range = N, assurance_for = "se",
    analysis_prior_se = an))
  setNames(cv$curve$assurance, cv$curve$N)
}

# ---- results table ----------------------------------------------------------
ROWS <- list()
add_row <- function(key, section, case, inputs, quantity, published = NA_real_,
                    source = NA_character_, package = NA_real_, emulation = NA_real_,
                    explanation = "") {
  ROWS[[length(ROWS) + 1L]] <<- data.frame(
    key = key, section = section, case = case, inputs = inputs, quantity = quantity,
    published = published, published_source = source, package = package,
    emulation = emulation, dataS2_original = NA_real_,
    stringsAsFactors = FALSE, explanation = explanation)
}
sec_time <- function(label, t0) cat(sprintf("[time] %s: %.1f s\n", label, as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# =============================================================================
# (a) VAP example, [W] section 5.1
# =============================================================================
if (want("a")) {
  t0 <- Sys.time()
  cat("==== (a) VAP, [W] section 5.1 ====\n")
  cat("Inputs: Se ~ Beta(25.9, 2.1), prevalence ~ Beta(29, 98), design = analysis,\n",
      "two-sided 95%, total width 0.16, assurance 0.80, Se only (assurance_for = 'se').\n", sep = "")
  se <- c(25.9, 2.1); pv <- c(29, 98); dw <- 0.16
  inp <- "Se Beta(25.9,2.1); prev Beta(29,98); design=analysis; total width 0.16; 95% two-sided; target 0.80; Se only"
  pk <- pk_nstar(se, pv, dw)
  pkA <- pk_A(se, pv, dw, c(104L, 105L, 106L, 150L))
  e <- new_arms(se, se, dw); extend_arms(e, 200L)
  fc <- first_cross(e, pv, "code", 200L, 1L, strict = TRUE)
  fx <- first_cross(e, pv, "exact", 200L, 1L, strict = FALSE)
  ac <- sapply(c(104, 105, 106, 150), A_arms, e = e, prev = pv, rule = "code")
  ax <- sapply(c(104, 105, 106, 150), A_arms, e = e, prev = pv, rule = "exact")
  cat(sprintf("package: n* = %d, A(n*) = %.6f; A(104,105,106,150) = %s\n", pk[["N"]], pk[["A"]],
              paste(sprintf("%.6f", pkA), collapse = ", ")))
  cat(sprintf("own exact enumeration: n* = %d; A(104,105,106,150) = %s; max |own - package| = %.2e\n",
              fx[["N"]], paste(sprintf("%.6f", ax), collapse = ", "), max(abs(ax - pkA))))
  cat(sprintf("emulation (code rule): n* = %d, A(n*) = %.6f; A(104,105,106,150) = %s\n", fc[["N"]], fc[["A"]],
              paste(sprintf("%.6f", ac), collapse = ", ")))
  add_row("VAP|n_star", "a", "VAP", inp, "n* (total N, first N with A >= 0.80)", 106,
          "[W] sec 5.1", pk[["N"]], fc[["N"]],
          "Difference of 2 explained by the success-set rule as programmed in the original code, not by the exact enumeration; the printed c1/c2 rule (exact) gives 104.")
  add_row("VAP|A_N104", "a", "VAP", inp, "A(104)", NA, NA, pkA[["104"]], ac[1L], "Tracing of the 2-unit difference.")
  add_row("VAP|A_N105", "a", "VAP", inp, "A(105)", NA, NA, pkA[["105"]], ac[2L], "Tracing of the 2-unit difference.")
  add_row("VAP|A_N106", "a", "VAP", inp, "A(106)", NA, NA, pkA[["106"]], ac[3L], "Code rule gives 0.8007 (> 0.80, first crossing); exact gives 0.8049.")
  add_row("VAP|A_N150", "a", "VAP", inp, "A(150)", 0.88, "[W] sec 5.1 (88%, rounded)", pkA[["150"]], ac[4L],
          "Agreement to the published rounding (0.8826 package, 0.8825 emulation).")
  fo <- optimize(function(N) abs(A_arms(e, N, pv, "code") - TARGET), c(10, 310))
  cat(sprintf("emulation of optimize(|A - 0.8|) on [10, 310]: minimum = %.6f, ceiling = %d\n", fo$minimum, ceiling(fo$minimum)))
  add_row("VAP|optimize_minimum", "a", "VAP (P10)", inp, "optimize(|A - 0.8|, [10, 310])$minimum (non-integer N allowed)", NA, NA, NA, fo$minimum,
          "Search of the original application script: continuous minimiser of |A - 0.8|, then ceiling.")
  add_row("VAP|optimize_ceiling", "a", "VAP (P10)", inp, "ceiling(optimize minimum)", 106, "[W] sec 5.1", NA, ceiling(fo$minimum),
          "Emulates the original search procedure (P10).")
  add_row("VAP|own_exact_vs_package", "a", "VAP", inp, "max |own exact enumeration - package| over A(104,105,106,150)", NA, NA, max(abs(ax - pkA)), NA,
          "Internal check of the helper (own enumeration vs package): must be <= 1e-12.")
  sec_time("(a)", t0)
}

# =============================================================================
# (b) Section 6.1 (Figure 3 text): 681, 1643, 2770
# =============================================================================
if (want("b")) {
  t0 <- Sys.time()
  cat("\n==== (b) [W] section 6.1 ====\n")
  cat("Inputs: Se only, two-sided 95%, total width 0.18, assurance 0.80, lambda = 0.9, n0 = 50;\n",
      "design Se ~ Beta(n0*lambda*rho, n0*rho*(1-lambda)), prevalence ~ Beta(n0*rho, n0*(1-rho));\n",
      "analysis prior = design prior (inference from the legend of Fig. 3 and the label 'BAM').\n", sep = "")
  rhos <- c(0.10, 0.05, 0.01); pub <- c(681, 1643, 2770)
  s61 <- list(); curves <- list()
  for (j in 1:3) {
    rho <- rhos[j]; n0 <- 50; lam <- 0.9
    se <- c(n0 * lam * rho, n0 * rho * (1 - lam)); pv <- c(n0 * rho, n0 * (1 - rho))
    inp <- sprintf("rho=%.2f; lambda=0.9; n0=50; Se Beta(%.3g,%.3g); prev Beta(%.3g,%.3g); design=analysis; total width 0.18",
                   rho, se[1], se[2], pv[1], pv[2])
    Mx <- if (rho == 0.01) 4000L else if (rho == 0.05) 1800L else 800L
    cat(sprintf("\n-- rho = %.2f (%s), search up to N = %d --\n", rho, inp, Mx))
    tt <- Sys.time()
    pk <- pk_nstar(se, pv, 0.18, N_max = Mx)
    cat(sprintf("package: n* = %s; A(n*) = %s; max A evaluated = %.6f at N = %d\n",
                format(pk[["N"]]), format(pk[["A"]], digits = 7), pk[["Amax"]], as.integer(pk[["Nmax"]])))
    e <- new_arms(se, se, 0.18)
    fc <- first_cross(e, pv, "code", Mx, 1L, strict = TRUE)
    fe <- first_cross(e, pv, "exact", Mx, 1L, strict = FALSE)
    cat(sprintf("emulation (code rule): n* = %s; own exact enumeration: n* = %s\n", format(fc[["N"]]), format(fe[["N"]])))
    extra <- if (rho == 0.10) 500 else if (rho == 0.05) c(1000, 1020) else numeric(0)
    Ns <- sort(unique(c(extra, pub[j] + (-1:1), if (!is.na(fc[["N"]])) fc[["N"]] + (-1:0), if (!is.na(pk[["N"]])) pk[["N"]] + (-1:0))))
    pkA <- pk_A(se, pv, 0.18, Ns)
    ac <- sapply(Ns, A_arms, e = e, prev = pv, rule = "code")
    ax <- sapply(Ns, A_arms, e = e, prev = pv, rule = "exact")
    cv <- data.frame(rho = rho, N = Ns, A_package = unname(pkA), A_own_exact = ax, A_code_rule = ac)
    curves[[j]] <- cv; print(cv, digits = 6, row.names = FALSE)
    tag <- sprintf("S61|rho%.2f", rho)
    add_row(paste0(tag, "|n_star"), "b", sprintf("rho=%.2f", rho), inp, "n* (total N)", pub[j], "[W] sec 6.1 text",
            pk[["N"]], fc[["N"]],
            if (rho == 0.10) "Code rule reproduces 681; exact enumeration gives 661 (rule of the original code)."
            else if (rho == 0.05) "Code rule reproduces 1643; exact enumeration gives 1408 (rule of the original code)."
            else "NOT REPRODUCED: with the declared priors neither rule crosses 0.80 up to N = 4000 (see P2 below).")
    for (N in pub[j] + (-1:0)) {
      if (rho == 0.01 && N == pub[j] - 1) next
      add_row(sprintf("%s|A_N%d", tag, N), "b", sprintf("rho=%.2f", rho), inp, sprintf("A(%d)", N), NA, NA,
              pkA[[as.character(N)]], ac[match(N, Ns)], "")
    }
    for (N in extra)
      add_row(sprintf("%s|A_N%d", tag, N), "b", sprintf("rho=%.2f", rho), inp, sprintf("A(%d)", N), NA, NA,
              pkA[[as.character(N)]], ac[match(N, Ns)],
              "N below the overflow of choose(N, k) in the original code: used to compare the emulation with the original code at large N.")
    if (rho == 0.01) {
      add_row(paste0(tag, "|A_N", pub[j] + 1), "b", sprintf("rho=%.2f", rho), inp, sprintf("A(%d)", pub[j] + 1), NA, NA,
              pkA[[as.character(pub[j] + 1)]], ac[match(pub[j] + 1, Ns)], "")
      add_row(paste0(tag, "|max_A_to_4000"), "b", sprintf("rho=%.2f", rho), inp, "max A over N = 2..4000 (package) / A(4000) (emulation)", NA, NA,
              pk[["Amax"]], A_arms(e, 4000, pv, "code"), "No crossing of 0.80 for either rule.")
    }
    s61[[j]] <- list(e = e, pv = pv, se = se)
    cat(sprintf("[time] rho = %.2f: %.1f s\n", rho, as.numeric(difftime(Sys.time(), tt, units = "secs"))))
  }
  write.csv(do.call(rbind, curves), file.path(outdir, "section6_1_curves.csv"), row.names = FALSE)

  # ---- P2: explanations tried for 2770 -------------------------------------
  cat("\n-- P2: attempts to explain 2770 (rho = 0.01) --\n")
  e <- s61[[3]]$e; pv <- s61[[3]]$pv
  # (a) optimize() of |A - 0.8| on a bounded interval: if A stays below 0.8 the
  #     minimiser of |A - 0.8| is the upper end of the interval.
  Agrid <- sapply(seq(100, 4000, by = 100), A_arms, e = e, prev = pv, rule = "code")
  cat("A_code on N = 100, 200, ..., 4000 (every 5th shown):\n")
  print(round(setNames(Agrid[seq(5, 40, by = 5)], seq(500, 4000, by = 500)), 4))
  cat(sprintf("A_code is %s on that grid; max = %.4f at N = %d\n",
              if (all(diff(Agrid) > 0)) "increasing" else "not monotone", max(Agrid), seq(100, 4000, by = 100)[which.max(Agrid)]))
  for (U in c(2000, 3000, 4000)) {
    fo <- optimize(function(N) abs(A_arms(e, N, pv, "code") - TARGET), c(10, U))
    cat(sprintf("  optimize(|A_code - 0.8|, [10, %d]) -> minimum = %.3f, ceiling = %d (A there = %.4f)\n",
                U, fo$minimum, ceiling(fo$minimum), A_arms(e, ceiling(fo$minimum), pv, "code")))
  }
  a2770c <- A_arms(e, 2770, pv, "code"); a2770x <- A_arms(e, 2770, pv, "exact")
  cat(sprintf("A(2770): code rule = %.4f, exact = %.4f\n", a2770c, a2770x))
  # (b) other priors: n0 in {25, 75, 100}, flat analysis prior; first crossing up to N = 4000
  alts <- list(list("n0 = 25, design = analysis", 25, "dsg"), list("n0 = 75, design = analysis", 75, "dsg"),
               list("n0 = 100, design = analysis", 100, "dsg"), list("n0 = 50, flat analysis Beta(1,1)", 50, "flat"))
  altrows <- list()
  for (al in alts) {
    n0 <- al[[2L]]; rho <- 0.01; lam <- 0.9
    se <- c(n0 * lam * rho, n0 * rho * (1 - lam)); pva <- c(n0 * rho, n0 * (1 - rho))
    an <- if (al[[3L]] == "flat") c(1, 1) else se
    ea <- new_arms(an, se, 0.18)
    ccode <- first_cross(ea, pva, "code", 4000L, 1L, strict = TRUE)
    cexact <- first_cross(ea, pva, "exact", 4000L, 1L, strict = FALSE)
    cat(sprintf("  alternative [%s]: first crossing code rule = %s, exact = %s (cap 4000); A_code(2770) = %.4f\n",
                al[[1L]], format(ccode[["N"]]), format(cexact[["N"]]), A_arms(ea, 2770, pva, "code")))
    altrows[[length(altrows) + 1L]] <- data.frame(alternative = al[[1L]], n_first_code = ccode[["N"]],
                                                  n_first_exact = cexact[["N"]], A_code_2770 = A_arms(ea, 2770, pva, "code"))
  }
  write.csv(do.call(rbind, altrows), file.path(outdir, "section6_1_P2_alternatives.csv"), row.names = FALSE)
  add_row("S61|rho0.01|P2", "b", "rho=0.01", "alternatives of P2 (see section6_1_P2_alternatives.csv)",
          "2770 explained?", 2770, "[W] sec 6.1 text", NA, NA,
          "Not reproduced. A(2770) = 0.355 (code rule) and 0.646 (exact) with the declared priors; neither crosses 0.80 up to N = 4000; optimize() on a bounded interval would return the upper end; the alternatives tried (n0 = 25, 75, 100; flat analysis prior) are in section6_1_P2_alternatives.csv. Not explained.")
  sec_time("(b)", t0)
}

# =============================================================================
# (c) Prior-data conflict (sections 4.3 and 5.3)
# =============================================================================
if (want("c")) {
  t0 <- Sys.time()
  cat("\n==== (c) prior-data conflict, [W] sections 4.3 and 5.3 ====\n")
  cat("Data (Table 1C): 51 TP, 2 FN (53 with VAP); 55 FP, 42 TN (97 without); n = 150.\n",
      "Priors: prevalence Beta(29, 98); Se Beta(25.9, 2.1) (design = analysis).\n", sep = "")
  cf <- suppressWarnings(bam_prior_data_conflict(tp = 51, fn = 2, tn = 42, fp = 55,
                          prior_se = c(25.9, 2.1), prior_prev = c(29, 98)))
  print(as.data.frame(cf)[, c("quantity", "reference", "x_obs", "n", "a", "b", "p_lower", "p_upper")], digits = 6)
  g <- function(q) cf[grepl(q, cf$quantity, fixed = TRUE), ][1L, ]
  pv_row <- g("prevalence"); se_row <- g("sensitivity")
  own <- function(x, n, a, b) c(lower = sum(bbpmf(0:x, n, a, b)), upper = sum(bbpmf(x:n, n, a, b)))
  ov <- own(53, 150, 29, 98); os <- own(51, 53, 25.9, 2.1)
  cat(sprintf("own cumulative sums: P(X<=53) = %.4f, P(X>=53) = %.4f; P(X<=51) = %.4f, P(X>=51) = %.4f\n", ov[1], ov[2], os[1], os[2]))
  inpc <- "tp=51 fn=2 tn=42 fp=55; prev Beta(29,98) BB(150); Se Beta(25.9,2.1) BB(53)"
  add_row("CONF|prev_pctl_lower", "c", "prevalence (number with VAP = 53 of 150)", inpc, "P(X <= 53) x 100 (percentile)", 99, "[W] sec 5.3 (99th)",
          100 * pv_row$p_lower, 100 * ov[["lower"]], "Published percentile corresponds to the LOWER cumulative P(X <= x_obs) (P5).")
  add_row("CONF|prev_tail_upper", "c", "prevalence", inpc, "P(X >= 53) (tail of [W] sec 4.3)", NA, NA, pv_row$p_upper, ov[["upper"]],
          "The upper tail defined in section 4.3 would be 0.012, not 0.99.")
  add_row("CONF|se_pctl_lower", "c", "Se (51 TP of 53)", inpc, "P(X <= 51) x 100 (percentile)", 76, "[W] sec 5.3 (76th)",
          100 * se_row$p_lower, 100 * os[["lower"]], "Published percentile corresponds to P(X <= x_obs) (0.759).")
  add_row("CONF|se_tail_upper", "c", "Se", inpc, "P(X >= 51) (tail of [W] sec 4.3)", NA, NA, se_row$p_upper, os[["upper"]],
          "The upper tail defined in section 4.3 would be 0.390.")
  sec_time("(c)", t0)
}

# =============================================================================
# (d) Table 2 of [W] (sensitivity to the initial priors, section 5.2)
# =============================================================================
grid_points <- function(par0, phi, r, eps = EPS) {
  out <- matrix(NA_real_, length(phi), 4L, dimnames = list(NULL, c("a", "b", "d", "r")))
  for (j in seq_along(phi)) {
    a <- par0[1L] + r * cos(phi[j]); b <- par0[2L] + r * sin(phi[j])
    ok <- a >= 0 & b >= 0
    a <- a[ok]; b <- b[ok]; rr <- r[ok]
    d <- 1 - exp(lbeta((a + par0[1L]) / 2, (b + par0[2L]) / 2) - 0.5 * (lbeta(a, b) + lbeta(par0[1L], par0[2L])))
    h <- which.min(abs(d - eps))
    out[j, ] <- c(a[h], b[h], d[h], rr[h])
  }
  out
}
if (want("d")) {
  t0 <- Sys.time()
  cat("\n==== (d) Table 2 of Wilson et al.: sensitivity to the initial prior (epsilon = 0.00354, distance H2 = 1 - BC) ====\n")
  se0 <- c(9.9, 1.1); pv0 <- c(12, 43); pilot <- list(tp = 16, fn = 1, tn = 20, fp = 35)
  add_se <- c(16, 1); add_pv <- c(17, 55); dw <- 0.16; Nref <- 106L
  se1 <- se0 + add_se; pv1 <- pv0 + add_pv
  cat("Initial priors: Se Beta(9.9,1.1), prevalence Beta(12,43); pilot 16 TP, 1 FN, 35 FP, 20 TN (Table 1B);\n",
      "perturbed initial prior + pilot = design = analysis prior of the perturbed parameter;\n",
      "A evaluated at N = 106 (the published n*); total width 0.16, Se only.\n", sep = "")
  phi63 <- seq(-pi, pi, 0.1); phi629 <- seq(-pi, pi, 0.01)
  cat(sprintf("number of angles: seq(-pi, pi, 0.1) = %d; seq(-pi, pi, 0.01) = %d\n", length(phi63), length(phi629)))
  ptab <- list(); summ <- list()
  for (which in c("Se", "prev")) {
    cat(sprintf("\n-- perturbation of the %s prior --\n", which))
    par0 <- if (which == "Se") se0 else pv0
    rg <- if (which == "Se") seq(0, 4, 0.1) else seq(0, 15, 0.1)
    G <- grid_points(par0, phi63, rg)
    add <- if (which == "Se") add_se else add_pv
    P1 <- sweep(G[, 1:2], 2, add, "+")
    cat(sprintf("grid rule (63 angles, r step 0.1): max |d - eps| = %.5f; a range [%.2f, %.2f], b range [%.2f, %.2f]\n",
                max(abs(G[, "d"] - EPS)), min(G[, "a"]), max(G[, "a"]), min(G[, "b"]), max(G[, "b"])))
    # A at N = 106, package and emulation (code rule) and own exact
    Apk <- Aco <- Aex <- numeric(nrow(P1))
    for (j in seq_len(nrow(P1))) {
      if (which == "Se") { sj <- P1[j, ]; pj <- pv1 } else { sj <- se1; pj <- P1[j, ] }
      Apk[j] <- pk_A(sj, pj, dw, Nref)[[1L]]
      e <- new_arms(sj, sj, dw); extend_arms(e, Nref)
      Aco[j] <- A_arms(e, Nref, pj, "code"); Aex[j] <- A_arms(e, Nref, pj, "exact")
    }
    cat(sprintf("max |own exact - package| over the 63 points = %.2e\n", max(abs(Aex - Apk))))
    ptab[[which]] <- data.frame(parameter = which, j = seq_len(nrow(G)), phi = phi63, r = G[, "r"], a_init = G[, "a"], b_init = G[, "b"],
                                d = G[, "d"], a_post = P1[, 1], b_post = P1[, 2], A106_package = Apk, A106_code_rule = Aco, A106_own_exact = Aex)
    jc <- which.min(Aco); jC <- which.max(Aco); jp <- which.min(Apk); jP <- which.max(Apk)
    # n* in the extreme priors: package (first crossing) and emulation (optimize on [10,310] and first crossing)
    prior_pair <- function(j) if (which == "Se") list(se = P1[j, ], pv = pv1) else list(se = se1, pv = P1[j, ])
    nstar_pk <- function(j) { pp <- prior_pair(j); pk_nstar(pp$se, pp$pv, dw)[["N"]] }
    nstar_em <- function(j) {
      pp <- prior_pair(j); e <- new_arms(pp$se, pp$se, dw); extend_arms(e, 310L)
      fo <- optimize(function(N) abs(A_arms(e, N, pp$pv, "code") - TARGET), c(10, 310))
      c(opt = ceiling(fo$minimum), first = first_cross(e, pp$pv, "code", 310L, 2L, strict = TRUE)[["N"]])
    }
    nA_min_pk <- nstar_pk(jp); nA_max_pk <- nstar_pk(jP)
    nA_min_em <- nstar_em(jc); nA_max_em <- nstar_em(jC)
    cat(sprintf("package on the grid: A in [%.4f, %.4f]; n* at A-min prior = %s, at A-max prior = %s\n",
                min(Apk), max(Apk), format(nA_min_pk), format(nA_max_pk)))
    cat(sprintf("emulation on the grid: A in [%.4f, %.4f]; n* (optimize, ceiling) at A-min = %s, at A-max = %s; (first crossing) %s and %s\n",
                min(Aco), max(Aco), format(nA_min_em[["opt"]]), format(nA_max_em[["opt"]]), format(nA_min_em[["first"]]), format(nA_max_em[["first"]])))
    cat(sprintf("A-min prior (code rule): initial (%.3f, %.3f); A-max prior: initial (%.3f, %.3f)\n",
                G[jc, "a"], G[jc, "b"], G[jC, "a"], G[jC, "b"]))
    pubA <- if (which == "Se") c(0.73, 0.86) else c(0.80, 0.81)
    pubn <- if (which == "Se") c(82, 130) else c(104, 108)
    src <- "[W] sec 5.2, Table 2"
    inp <- sprintf("%s prior perturbed; contour/grid epsilon = 0.00354 (H2); A at N = 106", which)
    # table rows: grid emulation vs published and package-on-grid
    add_row(sprintf("T2|%s|A_min", which), "d", paste0(which, ": Wilson grid (63 angles, r step 0.1)"), inp, "A(106) minimum", pubA[1], src,
            min(Apk), min(Aco), "Package column: exact rule at the 63 grid points of the published procedure; emulation: code rule at the same points.")
    add_row(sprintf("T2|%s|A_max", which), "d", paste0(which, ": Wilson grid (63 angles, r step 0.1)"), inp, "A(106) maximum", pubA[2], src,
            max(Apk), max(Aco), "Same.")
    add_row(sprintf("T2|%s|nstar_at_Amin", which), "d", paste0(which, ": Wilson grid"), inp, "n* at the prior with minimum A", pubn[2], src,
            nA_min_pk, nA_min_em[["opt"]], "Published range is computed only at the two priors of extreme A (not over the whole contour); emulation uses optimize() on [10,310] and ceiling as the original procedure.")
    add_row(sprintf("T2|%s|nstar_at_Amax", which), "d", paste0(which, ": Wilson grid"), inp, "n* at the prior with maximum A", pubn[1], src,
            nA_max_pk, nA_max_em[["opt"]], "Same.")
    add_row(sprintf("T2|%s|nstar_first_cross_at_Amin", which), "d", paste0(which, ": Wilson grid"), inp, "n* (first crossing, scan from 2) at the prior with minimum A", NA, NA,
            nA_min_pk, nA_min_em[["first"]], "Information: first crossing instead of optimize().")
    add_row(sprintf("T2|%s|nstar_first_cross_at_Amax", which), "d", paste0(which, ": Wilson grid"), inp, "n* (first crossing, scan from 2) at the prior with maximum A", NA, NA,
            nA_max_pk, nA_max_em[["first"]], "Information: first crossing instead of optimize().")
    # per-point keys for merging with Data S2
    for (j in seq_len(nrow(G))) add_row(sprintf("T2|%s|A_pt%02d", which, j), "d.points", paste0(which, " grid point ", j), inp,
                                        "A(106) at the grid point", NA, NA, Apk[j], Aco[j], "")
    # package: exact contour with 63 and 629 angles (bam_sensitivity_local)
    phl <- list(list("63 angles", phi63))
    if (!nzchar(Sys.getenv("WILSON_SKIP629"))) phl[[2L]] <- list("629 angles", phi629)
    for (ph in phl) {
      tt <- Sys.time()
      rs <- suppressWarnings(bam_sensitivity_local(
        base_prior_se = se0, base_prior_sp = c(2, 2), base_prior_prev = pv0,
        which = if (which == "Se") "se" else "prev", epsilon = EPS, phi = ph[[2L]], distance = "H2",
        pilot = pilot, delta_se = dw, assurance_for = "se", n_ref = Nref))
      s <- rs$summary
      tb <- rs$table
      n_at_Amin <- tb$n_star[which.min(tb$A_nref)]; n_at_Amax <- tb$n_star[which.max(tb$A_nref)]
      cat(sprintf("package, exact contour, %s: A in [%.4f, %.4f]; n* over the contour in [%s, %s]; n* at A-min = %s, at A-max = %s; directions without root: %s, imprecise: %s [%.0f s]\n",
                  ph[[1L]], s$A_min, s$A_max, format(s$n_star_min), format(s$n_star_max), format(n_at_Amin), format(n_at_Amax),
                  format(s$n_dir_no_root), format(s$n_dir_imprecise), as.numeric(difftime(Sys.time(), tt, units = "secs"))))
      cs <- paste0(which, ": exact contour, package, ", ph[[1L]])
      add_row(sprintf("T2|%s|contour%d|A_min", which, length(ph[[2L]])), "d", cs, inp, "A(106) minimum", pubA[1], src, s$A_min, NA, "Exact contour (root in r), not the nearest grid point.")
      # Explanations (2026-11-15): the published values use the nearest grid point and
      # the success-set rule as programmed in the original code; these rows use the exact contour and the
      # exact rule (the "Wilson grid" rows above reproduce the published values).
      ex_A <- "Exact contour (root in r) with the exact success-set rule; the published value uses the nearest grid point (63 angles, r step 0.1) and the rule as programmed, reproduced in the Wilson grid rows."
      ex_n <- "n* with the exact success-set rule (here over the whole contour or at the contour point of extreme A); the published n* are computed at the two grid priors of extreme A with the rule as programmed, reproduced in the Wilson grid rows."
      add_row(sprintf("T2|%s|contour%d|A_max", which, length(ph[[2L]])), "d", cs, inp, "A(106) maximum", pubA[2], src, s$A_max, NA, ex_A)
      add_row(sprintf("T2|%s|contour%d|nstar_min_over_contour", which, length(ph[[2L]])), "d", cs, inp, "n* minimum over the whole contour", pubn[1], src, s$n_star_min, NA, ex_n)
      add_row(sprintf("T2|%s|contour%d|nstar_max_over_contour", which, length(ph[[2L]])), "d", cs, inp, "n* maximum over the whole contour", pubn[2], src, s$n_star_max, NA, ex_n)
      add_row(sprintf("T2|%s|contour%d|nstar_at_Amin", which, length(ph[[2L]])), "d", cs, inp, "n* at the contour point of minimum A", pubn[2], src, n_at_Amin, NA, ex_n)
      add_row(sprintf("T2|%s|contour%d|nstar_at_Amax", which, length(ph[[2L]])), "d", cs, inp, "n* at the contour point of maximum A", pubn[1], src, n_at_Amax, NA, ex_n)
      summ[[length(summ) + 1L]] <- data.frame(parameter = which, angles = length(ph[[2L]]), A_min = s$A_min, A_max = s$A_max,
                                              n_star_min = s$n_star_min, n_star_max = s$n_star_max, n_dir_no_root = s$n_dir_no_root,
                                              n_dir_imprecise = s$n_dir_imprecise)
    }
  }
  write.csv(do.call(rbind, ptab), file.path(outdir, "table2_grid_points.csv"), row.names = FALSE)
  write.csv(do.call(rbind, summ), file.path(outdir, "table2_package_contours.csv"), row.names = FALSE)
  cat("\nSensitivity of the package result to the angular grid:\n"); print(do.call(rbind, summ), digits = 5, row.names = FALSE)
  sec_time("(d)", t0)
}

# =============================================================================
# (e) Distinct design and analysis priors: Supplement B.4 and Supplement C
# =============================================================================
if (want("e")) {
  t0 <- Sys.time()
  cat("\n==== (e) design and analysis priors distinct: [WS1] B.4 Table 1 and Supplement C ====\n")
  se <- c(25.9, 2.1); pv <- c(29, 98); dw <- 0.16
  inp <- "design Se Beta(25.9,2.1), prev Beta(29,98); total width 0.16; Se only"
  e1 <- new_arms(se, se, dw)         # analysis = design
  e2 <- new_arms(c(1, 1), se, dw)    # analysis flat
  f1 <- first_cross(e1, pv, "code", 400L, 1L, TRUE); f2 <- first_cross(e2, pv, "code", 400L, 1L, TRUE)
  p1 <- pk_nstar(se, pv, dw); p2 <- pk_nstar(se, pv, dw, an = c(1, 1))
  cat(sprintf("Standard (analysis = design): package %d, emulation %d; Combination / Hybrid II (flat Se analysis prior): package %d (A = %.4f), emulation %d\n",
              p1[["N"]], f1[["N"]], p2[["N"]], p2[["A"]], f2[["N"]]))
  cat("Hybrid I (informative Se analysis prior, flat prevalence analysis prior): the prevalence has no analysis prior in the package nor in the criterion, so the size equals the Standard one by construction.\n")
  add_row("B4|Standard|n_star", "e", "B.4 Standard (analysis = design)", inp, "n*", 106, "[WS1] B.4 Table 1", p1[["N"]], f1[["N"]], "Same 2-unit difference as the VAP example (code rule vs exact).")
  add_row("B4|Combination|n_star", "e", "B.4 Combination (Se Beta(1,1), prev Beta(1,1) analysis)", inp, "n*", 290, "[WS1] B.4 Table 1", p2[["N"]], f2[["N"]], "Agreement; A_exact(290) = 0.8004 (rule difference does not move the crossing).")
  add_row("B4|HybridI|n_star", "e", "B.4 Hybrid I (Se informative, prev Beta(1,1) analysis)", inp, "n*", 106, "[WS1] B.4 Table 1", p1[["N"]], f1[["N"]], "Prevalence has no role in the width criterion: equals Standard (published reading, [WS1] L299-300).")
  add_row("B4|HybridII|n_star", "e", "B.4 Hybrid II (Se Beta(1,1), prev informative analysis)", inp, "n*", 290, "[WS1] B.4 Table 1", p2[["N"]], f2[["N"]], "Equals Combination.")
  cat(sprintf("A(290) with flat analysis: package %.4f; emulation (code) %.4f\n", pk_A(se, pv, dw, 290L, an = c(1, 1))[[1L]], A_arms(e2, 290, pv, "code")))
  a290p <- pk_A(se, pv, dw, c(289L, 290L), an = c(1, 1))
  a290c <- sapply(c(289, 290), A_arms, e = e2, prev = pv, rule = "code")
  add_row("B4|A_N289_nonTRUE", "e", "B.4 flat analysis prior", inp, "A(289)", NA, NA, a290p[["289"]], a290c[1L], "")
  add_row("B4|A_N290_nonTRUE", "e", "B.4 flat analysis prior", inp, "A(290)", NA, NA, a290p[["290"]], a290c[2L], "Both rules cross 0.80 at 290.")
  # posterior means (arithmetic) of the B.4 table
  est <- c(Standard_Se = 76.9 / 81, Combination_Se = 52 / 55, Standard_prev = 82 / 277, Combination_prev = 54 / 152)
  pubest <- c(0.949, 0.945, 0.296, 0.355)
  for (i in seq_along(est)) add_row(paste0("B4|est|", names(est)[i]), "e", "B.4 posterior mean", "51/53 TP; 53 with VAP of 150", names(est)[i],
                                    pubest[i], "[WS1] B.4 Table 1", est[[i]], NA, "Arithmetic of the Beta-posterior means (not a sample-size result).")
  # ---- Supplement C ---------------------------------------------------------
  sc_txt <- "
w prev se bam25 bam50 bam75 noninf
0.22 0.9 0.70 56 28 1 79
0.22 0.9 0.75 50 20 1 73
0.22 0.9 0.80 41 9 1 66
0.22 0.9 0.85 28 1 1 56
0.22 0.7 0.70 82 53 24 103
0.22 0.7 0.75 74 43 12 96
0.22 0.7 0.80 62 30 2 86
0.22 0.7 0.85 48 12 2 74
0.22 0.5 0.70 132 99 68 149
0.22 0.5 0.75 120 86 55 138
0.22 0.5 0.80 104 69 36 124
0.22 0.5 0.85 83 46 13 106
0.22 0.3 0.70 256 214 179 259
0.22 0.3 0.75 235 191 157 240
0.22 0.3 0.80 207 164 126 215
0.22 0.3 0.85 171 125 88 185
0.18 0.9 0.70 97 68 40 119
0.18 0.9 0.75 88 57 28 110
0.18 0.9 0.80 75 43 13 98
0.18 0.9 0.85 57 24 1 83
0.18 0.9 0.90 35 1 1 65
0.18 0.7 0.70 136 105 77 156
0.18 0.7 0.75 125 92 62 144
0.18 0.7 0.80 108 74 43 128
0.18 0.7 0.85 86 50 17 109
0.18 0.7 0.90 57 19 2 85
0.18 0.5 0.70 211 175 144 224
0.18 0.5 0.75 194 157 124 207
0.18 0.5 0.80 171 131 98 184
0.18 0.5 0.85 139 99 64 157
0.18 0.5 0.90 97 56 23 123
0.18 0.3 0.70 399 348 310 390
0.18 0.3 0.75 367 316 277 360
0.18 0.3 0.80 326 275 234 322
0.18 0.3 0.85 272 220 178 274
0.18 0.3 0.90 199 149 108 215
0.14 0.9 0.70 178 147 119 197
0.14 0.9 0.75 164 131 101 182
0.14 0.9 0.80 143 108 77 161
0.14 0.9 0.85 115 79 47 136
0.14 0.9 0.90 79 42 1 104
0.14 0.7 0.70 243 210 180 258
0.14 0.7 0.75 224 188 157 238
0.14 0.7 0.80 198 159 127 211
0.14 0.7 0.85 163 122 88 178
0.14 0.7 0.90 115 75 41 137
0.14 0.5 0.70 368 326 292 372
0.14 0.5 0.75 340 296 260 343
0.14 0.5 0.80 302 256 218 305
0.14 0.5 0.85 251 203 166 257
0.14 0.5 0.90 183 136 99 198
0.14 0.3 0.70 680 614 569 648
0.14 0.3 0.75 629 562 515 598
0.14 0.3 0.80 561 495 445 534
0.14 0.3 0.85 473 406 357 452
0.14 0.3 0.90 352 291 244 349
"
  sc <- read.table(text = sc_txt, header = TRUE)
  if (nzchar(Sys.getenv("WILSON_SC_ROWS"))) sc <- sc[as.integer(strsplit(Sys.getenv("WILSON_SC_ROWS"), ",")[[1L]]), ]
  cat(sprintf("\nSupplement C ([WS1] Tables 2-4, L310-442): %d cells x 4 columns transcribed (BAM 25, 50, 75, Non-inf).\n", nrow(sc)))
  cat("Design: Se ~ Beta(n0*se*prev, n0*prev*(1-se)), prev ~ Beta(n0*prev, n0*(1-prev)), n0 = 25, 50, 75; BAM: analysis = design;\n",
      "Non-inf: analysis Beta(1,1), design prior with n0 = 75 (P3) or n0 = 25 (as described in [W] sec 6.1).\n", sep = "")
  flags_by_w <- list()
  res <- vector("list", nrow(sc) * 5L); ri <- 0L
  Mcap <- 1000L
  for (i in seq_len(nrow(sc))) {
    wi <- sc$w[i]; ri_prev <- sc$prev[i]; ls <- sc$se[i]
    wkey <- as.character(wi)
    if (is.null(flags_by_w[[wkey]])) flags_by_w[[wkey]] <- new.env()
    pvi <- c(75 * ri_prev, 75 * (1 - ri_prev))
    cols <- list(bam25 = 25, bam50 = 50, bam75 = 75, noninf = 75, noninf25 = 25)
    for (cn in names(cols)) {
      n0 <- cols[[cn]]
      sei <- c(n0 * ri_prev * ls, n0 * ri_prev * (1 - ls)); pvn <- c(n0 * ri_prev, n0 * (1 - ri_prev))
      flat <- cn %in% c("noninf", "noninf25")
      e <- new_arms(if (flat) c(1, 1) else sei, sei, wi, flags_cache = if (flat) flags_by_w[[wkey]] else NULL)
      fc <- first_cross(e, pvn, "code", Mcap, 1L, TRUE)
      fx <- first_cross(e, pvn, "exact", Mcap, 1L, FALSE)
      ri <- ri + 1L
      res[[ri]] <- data.frame(w = wi, prev = ri_prev, se = ls, column = cn, n0_design = n0,
                              published = if (cn == "noninf25") NA else sc[i, cn], emulation_code = fc[["N"]], own_exact_from1 = fx[["N"]],
                              package = NA_real_, stringsAsFactors = FALSE)
    }
    if (i %% 8L == 0L) cat(sprintf("  emulation: %d of %d cells done [%.0f s]\n", i, nrow(sc), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  scr <- do.call(rbind, res[seq_len(ri)])
  # package (exact enumeration) for the four published columns
  cat("package (exact) for the four published columns of every cell ...\n")
  for (r in seq_len(nrow(scr))) {
    if (scr$column[r] == "noninf25") next
    n0 <- scr$n0_design[r]; flat <- scr$column[r] == "noninf"
    sei <- c(n0 * scr$prev[r] * scr$se[r], n0 * scr$prev[r] * (1 - scr$se[r])); pvn <- c(n0 * scr$prev[r], n0 * (1 - scr$prev[r]))
    scr$package[r] <- pk_nstar(sei, pvn, scr$w[r], an = if (flat) c(1, 1) else NULL, N_max = 1500L)[["N"]]
    if (r %% 40L == 0L) cat(sprintf("  package: %d of %d rows [%.0f s]\n", r, nrow(scr), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  scr$dataS2_original <- NA_real_
  if (!is.na(ds2_csv) && file.exists(ds2_csv)) {
    d2c <- read.csv(ds2_csv, stringsAsFactors = FALSE)
    kk <- sprintf("SC|w%.2f|prev%.1f|se%.2f|%s|n_star", scr$w, scr$prev, scr$se, scr$column)
    scr$dataS2_original <- d2c$value[match(kk, d2c$key)]
  }
  write.csv(scr, file.path(outdir, "suppl_C_cells.csv"), row.names = FALSE)
  summ_c <- list()
  for (wv in c(0.22, 0.18, 0.14)) for (cn in c("bam25", "bam50", "bam75", "noninf", "noninf25")) {
    d <- scr[scr$w == wv & scr$column == cn, ]
    n <- nrow(d)
    m_code <- if (cn == "noninf25") NA else sum(d$emulation_code == d$published)
    m_pk <- if (cn == "noninf25") NA else sum(d$package == d$published | (d$published == 1 & d$package == 2 & d$emulation_code == 1))
    summ_c[[length(summ_c) + 1L]] <- data.frame(w = wv, column = cn, cells = n, match_code_rule = m_code, match_package = m_pk)
  }
  summ_c <- do.call(rbind, summ_c)
  cat("\nSupplement C: number of cells that coincide with the published value (published N = 1 cells counted as a match for the package when it returns 2, the scan starting at N = 2):\n")
  print(summ_c, row.names = FALSE)
  write.csv(summ_c, file.path(outdir, "suppl_C_summary.csv"), row.names = FALSE)
  # cells that differ
  dif <- scr[scr$column != "noninf25" & (scr$emulation_code != scr$published | (scr$package != scr$published & !(scr$published == 1 & scr$package == 2))), ]
  cat("\nCells where the emulation or the package differ from the published value:\n"); print(dif, row.names = FALSE)
  for (i in seq_len(nrow(summ_c))) {
    s <- summ_c[i, ]
    if (s$column == "noninf25") next
    add_row(sprintf("SC|w%.2f|%s|match", s$w, s$column), "e", sprintf("Suppl. C w* = %.2f, column %s", s$w, s$column),
            if (s$column == "noninf") "analysis Beta(1,1); design prior n0 = 75" else sprintf("design = analysis; n0 = %s", sub("bam", "", s$column)),
            "cells that coincide with the published value (of the table)", s$cells, "[WS1] C Tables 2-4", s$match_package, s$match_code_rule,
            if (s$column == "noninf") "Column coincides with a design prior of n0 = 75 (and not 25: see the noninf25 check in suppl_C_cells.csv), contrary to the n0 = 25 described in [W] sec 6.1 (P3)." else
              "The code-rule emulation matches every cell; the package (exact enumeration) differs by one unit in the other cells, explained by the success-set rule as programmed (cells listed in suppl_C_cells.csv).")
  }
  for (wv in c(0.22, 0.18, 0.14)) {
    d <- scr[scr$w == wv & scr$column == "noninf25", ]
    d$pub <- sc$noninf[sc$w == wv]
    add_row(sprintf("SC|w%.2f|noninf25|match", wv), "e", sprintf("Suppl. C w* = %.2f, Non-inf with design n0 = 25", wv), "analysis Beta(1,1); design prior n0 = 25 (as described in [W] sec 6.1)",
            "cells whose emulated size coincides with the published Non-inf", nrow(d), "[WS1] C Tables 2-4", NA, sum(d$emulation_code == d$pub), "Does not coincide (P3).")
  }
  # a few Data S2 keys: cells with published value differing between rules
  sec_time("(e)", t0)
}

# =============================================================================
# (f) Supplement F (joint Se and Sp)
# =============================================================================
if (want("f")) {
  cat("\n==== (f) Supplement F (joint Se and Sp) ====\n")
  cat("Not reproduced: [WS1] F gives only curves and box plots (no numerical sample sizes), and does not state the Sp prior,\n",
      "the Se grid step, the search range or the analysis prior. Figures were not digitised.\n",
      "Declared as a limitation: the joint two-sided mode of the package has no numerical external reproduction.\n", sep = "")
  add_row("SF|limitation", "f", "Supplement F", "Se, Sp joint; widths 0.18 and 0.22; prevalence 0.3; n0 = 50", "numerical reproduction", NA,
          "[WS1] F (figures only)", NA, NA, "Not reproducible: no numerical results and missing inputs. Declared as a limitation.")
}

# =============================================================================
# Merge with the original code run (optional) and write the table
# =============================================================================
tab <- do.call(rbind, ROWS)
if (!is.na(ds2_csv) && file.exists(ds2_csv)) {
  d2 <- read.csv(ds2_csv, stringsAsFactors = FALSE)
  m <- match(tab$key, d2$key)
  tab$dataS2_original <- d2$value[m]
  # derived rows from per-point values
  for (which in c("Se", "prev")) {
    ks <- sprintf("T2|%s|A_pt%02d", which, 1:63)
    v <- d2$value[match(ks, d2$key)]
    if (all(!is.na(v))) {
      tab$dataS2_original[tab$key == sprintf("T2|%s|A_min", which)] <- min(v)
      tab$dataS2_original[tab$key == sprintf("T2|%s|A_max", which)] <- max(v)
    }
  }
  cat(sprintf("\nData S2 original values merged for %d keys of %d in the file.\n", sum(!is.na(tab$dataS2_original)), nrow(d2)))
  if (exists("ptab")) for (pn in names(ptab)) {
    pt <- ptab[[pn]]
    ga <- d2$value[match(sprintf("T2|%s|grid_a_pt%02d", pn, pt$j), d2$key)]
    gb <- d2$value[match(sprintf("T2|%s|grid_b_pt%02d", pn, pt$j), d2$key)]
    ap <- d2$value[match(sprintf("T2|%s|A_pt%02d", pn, pt$j), d2$key)]
    cat(sprintf("Table 2 grid of Wilson et al., %s: max |a_emulated - a_DataS2| = %.2e, max |b| = %.2e; max |A_emulation - A_DataS2| = %.2e (63 points)\n",
                pn, max(abs(pt$a_init - ga)), max(abs(pt$b_init - gb)), max(abs(pt$A106_code_rule - ap))))
  }
  if (exists("scr")) {
    dd <- scr[!is.na(scr$dataS2_original), ]
    cat(sprintf("Supplement C cells with the original code: %d; emulation coincides in %d.\n", nrow(dd), sum(dd$emulation_code == dd$dataS2_original)))
  }
}
tab$diff_package_vs_published  <- tab$package   - tab$published
tab$diff_emulation_vs_published <- tab$emulation - tab$published
tab$diff_emulation_vs_dataS2   <- tab$emulation - tab$dataS2_original
cols <- c("key", "section", "case", "inputs", "quantity", "published", "published_source", "package", "emulation", "dataS2_original",
          "diff_package_vs_published", "diff_emulation_vs_published", "diff_emulation_vs_dataS2", "explanation")
tab2 <- tab[, cols]
write.csv(tab2, file.path(outdir, "compare_wilson2022_table.csv"), row.names = FALSE)
cat(sprintf("\nTable written: %d rows -> compare_wilson2022_table.csv\n", nrow(tab2)))
short <- tab2[tab2$section != "d.points", c("key", "published", "package", "emulation", "dataS2_original")]
print(short, digits = 6, row.names = FALSE)
cat(sprintf("\nTotal elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), T_START, units = "mins"))))
cat("sessionInfo (short):", R.version.string, "; dtasamplesize", as.character(packageVersion("dtasamplesize")), "\n")
