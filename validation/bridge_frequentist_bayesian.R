#!/usr/bin/env Rscript
# =============================================================================
# bridge_frequentist_bayesian.R
#
# Purpose
#   Approximate "bridge" demonstration between the frequentist joint-width
#   sample size and the Bayesian assurance sample size of the worked example
#   (N = 672, Se Beta(17,3), Sp Beta(18,2), prevalence Beta(4,16), full-width
#   targets 0.14 (Se) and 0.10 (Sp), target 0.80, 95% equal-tailed intervals).
#   It is a 2 x 2 factorial plus extra decomposition cells:
#     DESIGN prior   : point (Se = 0.85, Sp = 0.90, prevalence = 0.20) or
#                      Beta (the example priors);
#                      extras: only the prevalence uncertain, only Se and Sp
#                      uncertain.
#     ANALYSIS prior : vague (Jeffreys Beta(0.5,0.5) or uniform Beta(1,1), same
#                      for Se and Sp) or informative (the example priors, Se
#                      Beta(17,3) and Sp Beta(18,2)); reference row: frequentist
#                      Wilson score interval.
#   Frequentist exact reference: the enumeration of the Wilson-interval event
#   (manuscript v0.6.9, L277-278): N = 578.
#
# Usage
#   Rscript bridge_frequentist_bayesian.R <output_dir>
#   (writes ONLY in <output_dir>: log and CSV). Libraries: if the environment
#   variables DTA_PKG_LIB and DTA_EXTRA_LIB exist they are prepended to
#   .libPaths() (in that order). Single core, exact enumeration, no random
#   numbers.
#
# Conventions (named explicitly)
#   * Targets are FULL widths (upper minus lower limit) of the two-sided
#     interval: 0.14 (Se) and 0.10 (Sp); NOT half-widths, NOT one-sided.
#   * Bayesian interval: equal-tailed (quantiles 0.025 and 0.975) of the Beta
#     posterior Beta(a + x, b + n - x) under the ANALYSIS prior.
#   * Success rule: width <= target + TOL, TOL = 1e-12 (rule 11), for each arm
#     in the same replication (joint objective).
#   * POINT design prior (how it is represented): the design is NOT a
#     degenerate Beta; the true values are fixed (Se = 0.85, Sp = 0.90,
#     prevalence = 0.20, the means of the example priors 17/20, 18/20 and 4/20),
#     the number of diseased n_d is Binomial(N, 0.20) and each arm count is
#     Binomial(n_arm, true value). This is the design of joint_sample_size()
#     with design = "cohort" and of an earlier exploratory script (not distributed).
#   * BETA design prior: n_d ~ BetaBinomial(N, 4, 16); arm count ~
#     BetaBinomial(n_arm, a, b) (the package's exact method).
#   * Empty arm (n = 0): failure (package convention P_arm(0) = 0).
#   * N* = first integer N >= 2 with A(N) >= 0.80 (as the package; no
#     tolerance on A). Search up to N = N_CAP = 1200 (declared).
#   * Frequentist Wilson interval: score interval with z = qnorm(0.975), limits
#     clamped to [0, 1], as in joint_sample_size() and in the preliminary
#     script.
#
# Preliminary values (earlier exploratory calculation, not distributed;
# NOT evidence): point design with Wilson 578 (A = 0.8012447), Jeffreys 576
# (0.8014446), uniform 579 (0.8025248); point design with informative analysis
# priors 466. Article: Beta/informative = 672 (A = 0.8002692084, Table 2).
#
# Outputs (in <output_dir>)
#   bridge_frequentist_bayesian.log   text log (ends with the table and the sentence)
#   bridge_table.csv                  one row per design x analysis
#   bridge_A_by_N.csv                 A(N), N = 2..N_CAP, 15 significant digits
#   package_crosscheck.csv            Beta-design cells against bam_sample_size()
# =============================================================================

TOL <- 1e-12
N_CAP <- 1200L

.setup_libs <- function() {
  paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
  paths <- paths[nzchar(paths)]
  if (length(paths)) .libPaths(c(paths, .libPaths()))
}

main <- function() {
  .setup_libs()
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1)
  out <- args[1]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  zz <- file(file.path(out, "bridge_frequentist_bayesian.log"), open = "wt", encoding = "UTF-8")
  sink(zz, split = TRUE)
  on.exit({ sink(); close(zz) }, add = TRUE)
  t_start <- Sys.time()
  options(width = 170)
  suppressMessages(library(dtasamplesize))
  cat("bridge_frequentist_bayesian.R\n")
  cat("Date:", format(t_start, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("R:", R.version.string, "\n")
  cat("packageVersion(dtasamplesize):", as.character(utils::packageVersion("dtasamplesize")), "\n")
  cat("find.package(dtasamplesize):", find.package("dtasamplesize"), "\n")
  cat("TOL =", format(TOL), "; N_CAP =", N_CAP, "\n\n")

  z <- stats::qnorm(0.975)
  delta <- c(se = 0.14, sp = 0.10)
  target <- 0.80
  prior_se <- c(17, 3); prior_sp <- c(18, 2); prior_prev <- c(4, 16)
  pt <- c(se = 17 / 20, sp = 18 / 20, prev = 4 / 20)
  cat(sprintf("Point design values: Se = %.2f, Sp = %.2f, prevalence = %.2f (means of the example priors)\n\n",
              pt["se"], pt["sp"], pt["prev"]))

  # ---- widths (n = 0..N_CAP), success indicator per count: width <= delta + TOL
  wil_width <- function(n) {
    x <- 0:n
    if (n == 0L) return(NA_real_)
    p <- x / n; den <- 1 + z^2 / n
    c0 <- (p + z^2 / (2 * n)) / den
    m <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / den
    pmin(1, c0 + m) - pmax(0, c0 - m)
  }
  beta_width <- function(a, b) function(n) {
    x <- 0:n
    stats::qbeta(0.975, a + x, b + n - x) - stats::qbeta(0.025, a + x, b + n - x)
  }
  analyses <- list(
    wilson = list(label = "Frequentist Wilson score (reference)", se = wil_width, sp = wil_width),
    jeffreys = list(label = "Vague, Jeffreys Beta(0.5,0.5)", se = beta_width(.5, .5), sp = beta_width(.5, .5)),
    uniform = list(label = "Vague, uniform Beta(1,1)", se = beta_width(1, 1), sp = beta_width(1, 1)),
    informative = list(label = "Informative (example) Se Beta(17,3), Sp Beta(18,2)",
                       se = beta_width(17, 3), sp = beta_width(18, 2))
  )
  # shared widths where the two arms use the same analysis prior
  cat("Building width tables (qbeta enumeration up to n =", N_CAP, ") ...\n")
  t_w <- Sys.time()
  ok_list <- function(wfun, d) {
    lapply(0:N_CAP, function(n) if (n == 0L) logical(0) else wfun(n) <= d + TOL)
  }
  ok <- list()
  for (an in names(analyses)) {
    A <- analyses[[an]]
    if (an %in% c("wilson", "jeffreys", "uniform")) {
      wl <- lapply(0:N_CAP, function(n) if (n == 0L) numeric(0) else A$se(n))
      ok[[an]] <- list(se = lapply(wl, function(w) w <= delta["se"] + TOL),
                       sp = lapply(wl, function(w) w <= delta["sp"] + TOL))
    } else {
      ok[[an]] <- list(se = ok_list(A$se, delta["se"]), sp = ok_list(A$sp, delta["sp"]))
    }
  }
  cat(sprintf("  done in %.1f s\n\n", as.numeric(difftime(Sys.time(), t_w, units = "secs"))))

  # ---- design pmfs and per-arm success probabilities
  bb_pmf <- function(n, a, b) { x <- 0:n; exp(lchoose(n, x) + lbeta(a + x, b + n - x) - lbeta(a, b)) }
  pmf_fun <- function(kind, a, b, p) {
    if (kind == "point") function(n) stats::dbinom(0:n, n, p) else function(n) bb_pmf(n, a, b)
  }
  arm_P <- function(okl, pmf) {
    vapply(0:N_CAP, function(n) if (n == 0L) 0 else sum(pmf(n)[okl[[n + 1L]]]), numeric(1))
  }
  designs <- list(
    point = list(label = "Point (Se .85, Sp .90, prev .20)", se = "point", sp = "point", prev = "point"),
    beta = list(label = "Beta (example: Se B(17,3), Sp B(18,2), prev B(4,16))", se = "beta", sp = "beta", prev = "beta"),
    beta_prev_only = list(label = "Extra: only prevalence uncertain (Beta(4,16)); Se, Sp point", se = "point", sp = "point", prev = "beta"),
    beta_meas_only = list(label = "Extra: only Se and Sp uncertain (Beta); prevalence point .20", se = "beta", sp = "beta", prev = "point")
  )
  prev_pmf <- function(kind, N) {
    if (kind == "point") stats::dbinom(0:N, N, pt["prev"]) else bb_pmf(N, prior_prev[1], prior_prev[2])
  }

  rows <- list(); curves <- list(); Aall <- list()
  for (dn in names(designs)) {
    D <- designs[[dn]]
    pm_se <- pmf_fun(D$se, prior_se[1], prior_se[2], pt["se"])
    pm_sp <- pmf_fun(D$sp, prior_sp[1], prior_sp[2], pt["sp"])
    for (an in names(analyses)) {
      Pse <- arm_P(ok[[an]]$se, pm_se)
      Psp <- arm_P(ok[[an]]$sp, pm_sp)
      Ns <- 2:N_CAP
      A <- vapply(Ns, function(N) {
        k <- 0:N
        sum(prev_pmf(D$prev, N) * Pse[k + 1L] * Psp[N - k + 1L])
      }, numeric(1))
      i <- which(A >= target)
      nst <- if (length(i)) Ns[i[1L]] else NA_integer_
      rows[[length(rows) + 1L]] <- data.frame(
        design = dn, design_label = D$label, analysis = an, analysis_label = analyses[[an]]$label,
        N_star = nst, A_at_N_star = if (is.na(nst)) NA_real_ else A[Ns == nst],
        A_at_N_star_minus_1 = if (is.na(nst)) NA_real_ else A[Ns == nst - 1L],
        n_N_above_star_below_target = if (is.na(nst)) NA_integer_ else sum(A[Ns > nst] < target),
        max_A = max(A), stringsAsFactors = FALSE)
      Aall[[paste(dn, an)]] <- data.frame(design = dn, analysis = an, N = Ns,
                                          A = sprintf("%.15g", A), stringsAsFactors = FALSE)
    }
  }
  tab <- do.call(rbind, rows)
  utils::write.csv(tab, file.path(out, "bridge_table.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, Aall), file.path(out, "bridge_A_by_N.csv"), row.names = FALSE)

  # ---- cross-check of the Beta-design cells with the package (exact method)
  cat("==== Cross-check with bam_sample_size(method = 'exact') (Beta design) ====\n")
  options(dtasamplesize.warn_small_B = FALSE)
  xc <- list()
  chk <- list(informative = list(), jeffreys = list(analysis_prior_se = c(.5, .5), analysis_prior_sp = c(.5, .5)),
              uniform = list(analysis_prior_se = c(1, 1), analysis_prior_sp = c(1, 1)))
  for (an in names(chk)) {
    cl <- c(list(prior_se = prior_se, prior_sp = prior_sp, prior_prev = prior_prev,
                 delta_se = delta["se"], delta_sp = delta["sp"], target_assurance = target,
                 N_max = as.integer(N_CAP), B = 1000, n_range = 20:60), chk[[an]])
    cl$delta_se <- unname(cl$delta_se); cl$delta_sp <- unname(cl$delta_sp)
    t0 <- Sys.time()
    pk <- suppressWarnings(do.call(bam_sample_size, cl))
    mine <- tab[tab$design == "beta" & tab$analysis == an, ]
    ok_n <- identical(as.integer(pk$N_total), as.integer(mine$N_star))
    cat(sprintf("  Beta design / %-11s: package N_total = %s (A = %.10f) ; script N* = %s (A = %.10f) ; %s (%.0f s)\n",
                an, format(pk$N_total), pk$joint_assurance, format(mine$N_star), mine$A_at_N_star,
                if (ok_n) "MATCH" else "DIFFERENT", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    xc[[an]] <- data.frame(analysis = an, package_N_total = pk$N_total, package_A = pk$joint_assurance,
                           script_N_star = mine$N_star, script_A = mine$A_at_N_star,
                           abs_diff_A = abs(pk$joint_assurance - mine$A_at_N_star), match = ok_n)
  }
  utils::write.csv(do.call(rbind, xc), file.path(out, "package_crosscheck.csv"), row.names = FALSE)

  # ---- tables
  cat("\n==== Results: N* (first N with joint assurance >= 0.80) ====\n")
  show <- tab[, c("design", "analysis", "N_star", "A_at_N_star", "A_at_N_star_minus_1",
                  "n_N_above_star_below_target")]
  print(show, row.names = FALSE, digits = 10)

  g <- function(d, a) tab$N_star[tab$design == d & tab$analysis == a]
  prelim <- data.frame(
    cell = c("point / Wilson (frequentist reference)", "point / Jeffreys", "point / uniform",
             "point / informative", "Beta / informative"),
    preliminary = c(578, 576, 579, 466, 672),
    obtained = c(g("point", "wilson"), g("point", "jeffreys"), g("point", "uniform"),
                 g("point", "informative"), g("beta", "informative")),
    source = c("earlier exploratory calculation (not distributed)", "earlier exploratory calculation (not distributed)",
               "earlier exploratory calculation (not distributed)", "earlier exploratory calculation (not distributed)",
               "article, Table 2 (N = 672)"))
  prelim$difference <- prelim$obtained - prelim$preliminary
  cat("\n==== Preliminary values against obtained ====\n")
  print(prelim, row.names = FALSE)

  cat("\n==== 2 x 2 factorial (N*), rows = design, columns = analysis ====\n")
  f <- data.frame(design = c("point", "Beta"),
                  vague_Jeffreys = c(g("point", "jeffreys"), g("beta", "jeffreys")),
                  vague_uniform = c(g("point", "uniform"), g("beta", "uniform")),
                  informative = c(g("point", "informative"), g("beta", "informative")),
                  Wilson_frequentist = c(g("point", "wilson"), g("beta", "wilson")))
  print(f, row.names = FALSE)
  cat("(The 'Beta / Wilson' cell is the Wilson-interval event under the Beta design; it is an extra, not part of the 2 x 2.)\n")
  cat("\nExtra decomposition cells (design uncertainty by source):\n")
  ex <- tab[tab$design %in% c("beta_prev_only", "beta_meas_only"), c("design", "analysis", "N_star", "A_at_N_star")]
  print(ex, row.names = FALSE, digits = 10)

  # ---- sentence for the text
  f0 <- g("point", "wilson"); fj <- g("point", "jeffreys"); fu <- g("point", "uniform")
  bv_j <- g("beta", "jeffreys"); bv_u <- g("beta", "uniform")
  pi_ <- g("point", "informative"); bi <- g("beta", "informative")
  cat("\n==== Sentence for the text (numbers as obtained) ====\n")
  s1 <- sprintf(paste0(
    "With a point design prior, the frequentist Wilson-interval event gives N = %d, and Bayesian credible intervals ",
    "with a vague analysis prior give N = %d (Jeffreys) and %d (uniform), so the interval type alone moves N by at most %d. ",
    "Replacing the point design by the Beta design priors of the example, with the same vague analysis priors, changes N to %d (Jeffreys) and %d (uniform), ",
    "a change of %+d and %+d relative to the point-design values (effect of design uncertainty). ",
    "Replacing the vague analysis priors by the informative priors of the example, with the point design, changes N from %d-%d to %d ",
    "(effect of the informative analysis prior, %+d to %+d). ",
    "The worked example combines both changes (Beta design, informative analysis): N = %d, ",
    "that is %+d relative to the frequentist reference %d; the design uncertainty and the informative analysis prior act in opposite directions. ",
    "Along the path point design with Jeffreys analysis, Beta design with Jeffreys analysis, Beta design with informative analysis: %d, %d, %d (%+d, then %+d)."),
    f0, fj, fu, max(abs(c(fj, fu) - f0)),
    bv_j, bv_u, bv_j - fj, bv_u - fu,
    fj, fu, pi_, pi_ - fj, pi_ - fu,
    bi, bi - f0, f0,
    fj, bv_j, bi, bv_j - fj, bi - bv_j)
  cat(strwrap(s1, width = 110), sep = "\n")
  cat("\n(The sentence is generated from the table; check the sign of each effect before reusing the wording: ",
      sprintf("design effect on vague analysis = %+d / %+d, analysis effect under point design = %+d / %+d, ",
              bv_j - fj, bv_u - fu, pi_ - fj, pi_ - fu),
      sprintf("analysis effect under Beta design = %+d / %+d.)\n", bi - bv_j, bi - bv_u), sep = "")
  cat(sprintf("\nTotal elapsed: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  cat("Output directory:", normalizePath(out, winslash = "/"), "\n")
}

main()
