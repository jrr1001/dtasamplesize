#!/usr/bin/env Rscript
# =============================================================================
# exact_vs_mc_resolution.R
#
# Purpose
#   Value of the EXACT evaluation of the joint assurance against the package's
#   own Monte Carlo estimate (method = "monte_carlo"), in the two worked
#   examples of the article: P1 (N* = 678, vague Sp Beta(2,2)) and P2
#   (N* = 672, informative Sp Beta(18,2)). It is NOT a comparison with Wilson et
#   al. or Williamson et al., which also evaluate in closed form.
#   Reported:
#     (1) distance of A(N*) and A(N*-1) to the target 0.80;
#     (2) replicates B needed so that 1.96 MCSE is smaller than that distance,
#         i.e. to resolve the integer N;
#     (3) range of N compatible with 0.80 under typical B (2000, 10000, 50000):
#         N whose exact A lies within 1.96 MCSE(B) of the target;
#     (4) effect of the package's anti-noise acceptance rule (Monte Carlo mode)
#         on the N returned: run with declared seeds and several B, against
#         the exact N*;
#     (5) computing times (exact vs Monte Carlo with each B).
#
# Usage
#   Rscript exact_vs_mc_resolution.R <output_dir>
#   (writes ONLY in <output_dir>: log and CSV). Libraries: if the environment
#   variables DTA_PKG_LIB and DTA_EXTRA_LIB exist they are prepended to
#   .libPaths() (in that order). Single core.
#
# Conventions (named explicitly)
#   * Targets are FULL widths of the equal-tailed 95% credible interval
#     (width_type = "full_width"): 0.14 (Se) and 0.10 (Sp); two-sided; joint
#     objective (Se and Sp widths met in the same replication); design prior =
#     analysis prior; target assurance 0.80.
#   * Priors: P1 Se Beta(17,3), Sp Beta(2,2), prevalence Beta(4,16); P2 as P1
#     with Sp Beta(18,2) (as validation/compat_v069_check.R).
#   * Width comparison with TOL = 1e-12 (rule 11): width <= delta + TOL in the
#     script's own enumeration; the package compares width <= delta exactly.
#   * Empty arm: failure (P_arm(0) = 0), as the package.
#   * N* = first N >= 2 with A(N) >= 0.80 (no tolerance on A), as the package.
#   * Anti-noise rule of the package (R/bam_sample_size.R, Details, "Anti-noise
#     acceptance rule", method = "monte_carlo" only): an N is accepted only if
#     the lower bound of the ONE-SIDED 95% confidence interval of its Monte
#     Carlo estimate reaches the target:
#         A_hat - qnorm(0.95) * MCSE >= target,  MCSE = sqrt(A_hat (1 - A_hat) / B).
#     (qnorm(0.95) = 1.645.) The exact method does not apply it.
#   * Resolution criterion (2): the preliminary convention of the plan: 1.96
#     MCSE with p(1 - p) = 0.16 (p = target 0.80). Three figures: B_min_margin
#     (1.96 MCSE < the nearer of the two distances, resolves N* against N*-1 AND
#     against the target on both sides), B_max_margin (the farther one) and
#     B_half_diff (MCSE < half of A(N*) - A(N*-1)).
#   * Compatible range (3): N with |A(N) - 0.80| <= 1.96 sqrt(A(N)(1 - A(N))/B).
#   * Anti-noise proxy (4): A_req(B) solves A - 1.645 sqrt(A(1-A)/B) = 0.80;
#     N_rule(B) = first N with exact A(N) >= A_req(B) (deterministic proxy of
#     the N the rule would return if the Monte Carlo estimate were exactly the
#     exact probability). The returned N differs from it by Monte Carlo noise.
#   * Monte Carlo seeds (declared, fresh; not the registered seeds of the
#     LOCKED grid): B = 2000: 11, 22, 33, 44, 55; B = 10000: 11, 22, 33;
#     B = 50000: 11, 22, 33. Additionally seed 1 with the automatic grid and with
#     the explicit window, B = 2000 and 10000 (consistency check of the window).
#   * The Monte Carlo calls use an explicit N_range (a window around N*,
#     declared per B below) to save time. The package resets the seed at every N
#     (R/bam_sample_size.R, search loop), so the estimate at each N does not
#     depend on the grid; the window returns the same N_total as the automatic
#     grid whenever the crossing lies inside it (checked with seed 1).
#   * Timings depend on the machine and on other jobs running at the same time;
#     they are only meaningful as orders of magnitude. The exact timing
#     includes the legacy per-arm fields (B = 5000 default) in the "default"
#     row; the "minimal" row sets B = 1000 and n_range = 20:60.
#
# Outputs (in <output_dir>)
#   exact_vs_mc_resolution.log   text log
#   exact_A_near_target.csv      A(N) for N = 2..N_cap, 15 significant digits
#   distances_and_replicates.csv (1) and (2)
#   compatible_range.csv         (3)
#   mc_runs.csv                  (4): every Monte Carlo call
#   mc_summary.csv               (4): by scenario and B
#   timings.csv                  (5)
# =============================================================================

TOL <- 1e-12
N_CAP <- 900L
set.seed(20261115)   # no random numbers are drawn by this script itself; the package seeds its own runs

.setup_libs <- function() {
  paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
  paths <- paths[nzchar(paths)]
  if (length(paths)) .libPaths(c(paths, .libPaths()))
}

bb_pmf <- function(n, a, b) { x <- 0:n; exp(lchoose(n, x) + lbeta(a + x, b + n - x) - lbeta(a, b)) }

# P(width <= delta | arm size n), n = 0..N_CAP; P(0) = 0.
arm_prob <- function(a, b, delta, level = 0.95) {
  lq <- (1 - level) / 2
  vapply(0:N_CAP, function(n) {
    if (n == 0L) return(0)
    x <- 0:n
    w <- stats::qbeta(1 - lq, a + x, b + n - x) - stats::qbeta(lq, a + x, b + n - x)
    sum(bb_pmf(n, a, b)[w <= delta + TOL])
  }, numeric(1))
}

exact_A <- function(Pse, Psp, prev) {
  Ns <- 2:N_CAP
  A <- vapply(Ns, function(N) {
    k <- 0:N
    sum(bb_pmf(N, prev[1], prev[2]) * Pse[k + 1L] * Psp[N - k + 1L])
  }, numeric(1))
  data.frame(N = Ns, A = A)
}

main <- function() {
  .setup_libs()
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1)
  out <- args[1]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  zz <- file(file.path(out, "exact_vs_mc_resolution.log"), open = "wt", encoding = "UTF-8")
  sink(zz, split = TRUE)
  on.exit({ sink(); close(zz) }, add = TRUE)
  t_start <- Sys.time()
  suppressMessages(library(dtasamplesize))
  options(dtasamplesize.warn_small_B = FALSE, width = 170)
  cat("exact_vs_mc_resolution.R\n")
  cat("Date:", format(t_start, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("R:", R.version.string, "\n")
  cat("packageVersion(dtasamplesize):", as.character(utils::packageVersion("dtasamplesize")), "\n")
  cat("find.package(dtasamplesize):", find.package("dtasamplesize"), "\n")
  cat("TOL =", format(TOL), "; N_CAP =", N_CAP, "\n\n")

  target <- 0.80
  scen <- list(
    P1 = list(prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16), prelim = c(N = 678, d1 = 3.49e-4, d0 = 3.15e-4)),
    P2 = list(prior_se = c(17, 3), prior_sp = c(18, 2), prior_prev = c(4, 16), prelim = c(N = 672, d1 = NA, d0 = NA))
  )
  base_args <- function(s) list(prior_se = s$prior_se, prior_sp = s$prior_sp, prior_prev = s$prior_prev,
                                delta_se = 0.14, delta_sp = 0.10, target_assurance = target)
  Bs_typ <- c(2000, 10000, 50000)
  quick <- identical(Sys.getenv("DTA_QUICK"), "1")
  if (quick) scen <- scen["P1"]

  # ---------------------------------------------------------------- exact ---
  cat("==== (1) Exact enumeration (own implementation), cross-checked with the package ====\n")
  ex <- list(); Acurves <- list(); dist_rows <- list(); tim <- list()
  for (id in names(scen)) {
    s <- scen[[id]]
    t0 <- Sys.time()
    Pse <- arm_prob(s$prior_se[1], s$prior_se[2], 0.14)
    Psp <- arm_prob(s$prior_sp[1], s$prior_sp[2], 0.10)
    cur <- exact_A(Pse, Psp, s$prior_prev)
    t_own <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    i <- which(cur$A >= target)[1]
    Nst <- cur$N[i]
    A1 <- cur$A[i]; A0 <- cur$A[i - 1L]
    ex[[id]] <- list(cur = cur, Nst = Nst, A1 = A1, A0 = A0, Pse = Pse, Psp = Psp)
    Acurves[[id]] <- data.frame(scenario = id, N = cur$N, A = sprintf("%.15g", cur$A), stringsAsFactors = FALSE)
    # package, default call (timed)
    ba <- base_args(s)
    tm_def <- system.time(pk <- suppressWarnings(do.call(bam_sample_size, ba)))[["elapsed"]]
    tm_min <- system.time(pk_min <- suppressWarnings(do.call(bam_sample_size,
                            c(ba, list(B = 1000, n_range = 20:60)))))[["elapsed"]]
    pk0 <- suppressWarnings(do.call(bam_sample_size, c(ba, list(N_range = Nst - 1L, B = 1000, n_range = 20:60))))
    cat(sprintf("%s: own enumeration N* = %d, A(N*) = %.10f, A(N*-1) = %.10f (%.1f s)\n", id, Nst, A1, A0, t_own))
    cat(sprintf("    package exact: N_total = %s, A(N*) = %.10f (|diff| = %.2e); A(N*-1) = %.10f (|diff| = %.2e)\n",
                format(pk$N_total), pk$joint_assurance, abs(pk$joint_assurance - A1),
                pk0$max_assurance_evaluated, abs(pk0$max_assurance_evaluated - A0)))
    d1 <- A1 - target; d0 <- target - A0
    cat(sprintf("    distance to target: A(N*) - 0.80 = %.4e ; 0.80 - A(N*-1) = %.4e ; A(N*) - A(N*-1) = %.4e\n", d1, d0, A1 - A0))
    if (!is.na(s$prelim["d1"]))
      cat(sprintf("    preliminary (plan): %.2e and %.2e (published in the plan as 3.49e-4 and 3.15e-4)\n", s$prelim["d1"], s$prelim["d0"]))
    dmin <- min(d1, d0); dmax <- max(d1, d0)
    Bmin <- 0.16 * (1.96 / dmin)^2
    Bmax <- 0.16 * (1.96 / dmax)^2
    Bhalf <- 0.16 / ((A1 - A0) / 2)^2
    cat(sprintf("(2) replicates for 1.96 MCSE < nearer distance: B_min_margin = %.0f ; < farther distance: B_max_margin = %.0f ; MCSE < half(A(N*)-A(N*-1)): B_half_diff = %.0f\n",
                Bmin, Bmax, Bhalf))
    cat(sprintf("    1.96 MCSE at typical B (p = 0.80): B = 2000: %.4f ; 10000: %.4f ; 50000: %.4f ; ratio to the nearer distance: %.0f, %.0f, %.0f\n",
                1.96 * sqrt(0.16 / 2000), 1.96 * sqrt(0.16 / 10000), 1.96 * sqrt(0.16 / 50000),
                1.96 * sqrt(0.16 / 2000) / dmin, 1.96 * sqrt(0.16 / 10000) / dmin, 1.96 * sqrt(0.16 / 50000) / dmin))
    dist_rows[[id]] <- data.frame(scenario = id, N_star = Nst, A_N_star = A1, A_N_star_minus_1 = A0,
      dist_N_star = d1, dist_N_star_minus_1 = d0, diff_A = A1 - A0,
      package_N_total = pk$N_total, package_A_N_star = pk$joint_assurance,
      package_A_N_star_minus_1 = pk0$max_assurance_evaluated,
      B_min_margin = Bmin, B_max_margin = Bmax, B_half_diff = Bhalf,
      plan_dist_N_star = unname(s$prelim["d1"]), plan_dist_N_star_minus_1 = unname(s$prelim["d0"]), row.names = NULL)
    tim[[length(tim) + 1L]] <- data.frame(scenario = id, method = "exact (own enumeration up to N = 900)", B = NA, grid = "N = 2..900",
                                           elapsed_s = t_own, note = "script's own qbeta enumeration, no legacy fields")
    tim[[length(tim) + 1L]] <- data.frame(scenario = id, method = "exact (bam_sample_size defaults)", B = 5000, grid = "automatic integer scan",
                                           elapsed_s = tm_def, note = "includes legacy per-arm Monte Carlo fields (B = 5000, n_range 20:500)")
    tim[[length(tim) + 1L]] <- data.frame(scenario = id, method = "exact (bam_sample_size, B = 1000, n_range = 20:60)", B = 1000,
                                           grid = "automatic integer scan", elapsed_s = tm_min,
                                           note = "minimal legacy fields; the headline exact result does not depend on B")
    cat(sprintf("    timing: own enumeration %.1f s; bam_sample_size exact defaults %.1f s; minimal-legacy %.1f s\n\n", t_own, tm_def, tm_min))
  }
  utils::write.csv(do.call(rbind, Acurves), file.path(out, "exact_A_near_target.csv"), row.names = FALSE)
  utils::write.csv(do.call(rbind, dist_rows), file.path(out, "distances_and_replicates.csv"), row.names = FALSE)

  # ------------------------------------------------------- compatible range --
  cat("==== (3) Range of N compatible with 0.80 at typical B (|A(N) - 0.80| <= 1.96 MCSE(B), MCSE at A(N)) ====\n")
  comp <- list()
  for (id in names(scen)) {
    cur <- ex[[id]]$cur
    for (B in Bs_typ) {
      ok <- abs(cur$A - target) <= 1.96 * sqrt(cur$A * (1 - cur$A) / B)
      rng <- range(cur$N[ok])
      # contiguous block around N*
      comp[[length(comp) + 1L]] <- data.frame(scenario = id, B = B, N_min = rng[1], N_max = rng[2],
        n_values = sum(ok), N_star = ex[[id]]$Nst,
        below = ex[[id]]$Nst - rng[1], above = rng[2] - ex[[id]]$Nst, contiguous = (rng[2] - rng[1] + 1) == sum(ok))
      cat(sprintf("  %s, B = %5d: N in [%d, %d] (%d values; N* - %d, N* + %d)\n", id, B, rng[1], rng[2],
                  sum(ok), ex[[id]]$Nst - rng[1], rng[2] - ex[[id]]$Nst))
    }
  }
  comp <- do.call(rbind, comp)
  utils::write.csv(comp, file.path(out, "compatible_range.csv"), row.names = FALSE)

  # --------------------------------------------------------- anti-noise rule -
  cat("\n==== (4) Anti-noise acceptance rule: A_hat - qnorm(0.95) MCSE >= 0.80 (package, Monte Carlo mode) ====\n")
  cat("Deterministic proxy: A_req(B) solves A - 1.645 sqrt(A(1-A)/B) = 0.80; N_rule(B) = first N with exact A(N) >= A_req(B).\n")
  Areq <- function(B) {
    f <- function(A) A - stats::qnorm(0.95) * sqrt(A * (1 - A) / B) - target
    stats::uniroot(f, c(target, 0.999), tol = 1e-13)$root
  }
  proxy <- list()
  for (id in names(scen)) {
    cur <- ex[[id]]$cur
    for (B in if (quick) c(Bs_typ, 500) else Bs_typ) {
      ar <- Areq(B)
      nr <- cur$N[which(cur$A >= ar)[1]]
      proxy[[paste(id, B)]] <- data.frame(scenario = id, B = B, A_required = ar, N_rule_proxy = nr,
                                          N_star_exact = ex[[id]]$Nst, shift_vs_exact = nr - ex[[id]]$Nst)
      cat(sprintf("  %s, B = %5d: A_req = %.5f ; N_rule(B) = %d ; exact N* = %d ; shift = %+d\n", id, B, ar, nr,
                  ex[[id]]$Nst, nr - ex[[id]]$Nst))
    }
  }
  proxy <- do.call(rbind, proxy)

  windows <- function(Nst, B) if (B >= 50000) (Nst - 40L):(Nst + 80L) else (Nst - 60L):(Nst + 160L)
  runs <- list()
  plan <- list(list(B = 2000, seeds = c(11, 22, 33, 44, 55)),
               list(B = 10000, seeds = c(11, 22, 33)),
               list(B = 50000, seeds = c(11, 22, 33)))
  check_Bs <- c(2000, 10000)
  if (quick) {   # smoke test only (DTA_QUICK=1); the reported run never sets it
    plan <- list(list(B = 500, seeds = 11)); check_Bs <- 500
  }
  for (id in names(scen)) {
    s <- scen[[id]]; ba <- base_args(s); Nst <- ex[[id]]$Nst; cur <- ex[[id]]$cur
    # consistency check: automatic grid vs window, seed 1 (also the automatic-grid timings)
    for (B in check_Bs) {
      tm_auto <- system.time(ma <- suppressWarnings(do.call(bam_sample_size,
                  c(ba, list(method = "monte_carlo", B = B, seed = 1)))))[["elapsed"]]
      tm_win <- system.time(mw <- suppressWarnings(do.call(bam_sample_size,
                  c(ba, list(method = "monte_carlo", B = B, seed = 1, N_range = windows(Nst, B))))))[["elapsed"]]
      cat(sprintf("  %s check seed 1, B = %5d: automatic grid N = %s (%.0f s) ; window N = %s (%.0f s) ; %s\n",
                  id, B, format(ma$N_total), tm_auto, format(mw$N_total), tm_win,
                  if (identical(ma$N_total, mw$N_total)) "SAME" else "DIFFERENT"))
      tim[[length(tim) + 1L]] <- data.frame(scenario = id, method = "monte_carlo automatic grid, seed 1", B = B,
        grid = paste0(min(ma$N_range_used), ":", max(ma$N_range_used)), elapsed_s = tm_auto,
        note = "bam_sample_size defaults otherwise (includes legacy per-arm fields at this B)")
      tim[[length(tim) + 1L]] <- data.frame(scenario = id, method = "monte_carlo explicit window, seed 1", B = B,
        grid = paste0(min(windows(Nst, B)), ":", max(windows(Nst, B))), elapsed_s = tm_win,
        note = "window around N*; includes legacy per-arm fields at this B")
      runs[[length(runs) + 1L]] <- data.frame(scenario = id, B = B, seed = 1, mode = "automatic grid",
        N_returned = ma$N_total, A_hat = ma$joint_assurance, mcse = ma$assurance_mcse,
        A_exact_at_returned = if (is.na(ma$N_total)) NA_real_ else cur$A[cur$N == ma$N_total],
        N_exact = Nst, elapsed_s = tm_auto)
    }
    for (pl in plan) {
      for (sd in pl$seeds) {
        wdw <- windows(Nst, pl$B)
        tm <- system.time(m <- suppressWarnings(do.call(bam_sample_size,
                c(ba, list(method = "monte_carlo", B = pl$B, seed = sd, N_range = wdw)))))[["elapsed"]]
        nr <- m$N_total
        runs[[length(runs) + 1L]] <- data.frame(scenario = id, B = pl$B, seed = sd, mode = "window",
          N_returned = nr, A_hat = m$joint_assurance, mcse = m$assurance_mcse,
          A_exact_at_returned = if (is.na(nr)) NA_real_ else cur$A[cur$N == nr],
          N_exact = Nst, elapsed_s = tm)
        cat(sprintf("  %s, B = %5d, seed %2d: N returned = %s ; A_hat = %.4f (MCSE %.4f) ; exact A there = %s ; N - N* = %s ; %.0f s\n",
                    id, pl$B, sd, format(nr), m$joint_assurance, m$assurance_mcse,
                    if (is.na(nr)) "NA" else sprintf("%.4f", cur$A[cur$N == nr]),
                    if (is.na(nr)) "NA" else sprintf("%+d", nr - Nst), tm))
      }
    }
  }
  runs <- do.call(rbind, runs)
  utils::write.csv(runs, file.path(out, "mc_runs.csv"), row.names = FALSE)
  summ <- do.call(rbind, lapply(split(runs[runs$mode == "window", ], list(runs$scenario[runs$mode == "window"], runs$B[runs$mode == "window"])), function(d) {
    pr <- proxy[proxy$scenario == d$scenario[1] & proxy$B == d$B[1], ]
    data.frame(scenario = d$scenario[1], B = d$B[1], n_seeds = nrow(d), N_exact = d$N_exact[1],
               N_returned_min = min(d$N_returned, na.rm = TRUE), N_returned_median = stats::median(d$N_returned, na.rm = TRUE),
               N_returned_max = max(d$N_returned, na.rm = TRUE), n_NA = sum(is.na(d$N_returned)),
               shift_median_vs_exact = stats::median(d$N_returned, na.rm = TRUE) - d$N_exact[1],
               N_rule_proxy = pr$N_rule_proxy, A_hat_median = stats::median(d$A_hat, na.rm = TRUE),
               A_exact_at_returned_median = stats::median(d$A_exact_at_returned, na.rm = TRUE))
  }))
  rownames(summ) <- NULL
  utils::write.csv(summ, file.path(out, "mc_summary.csv"), row.names = FALSE)
  cat("\nSummary of the Monte Carlo runs (explicit window, declared seeds):\n")
  print(summ, row.names = FALSE, digits = 6)

  tim <- do.call(rbind, tim)
  # window timings of the declared-seed runs
  tw <- aggregate(elapsed_s ~ scenario + B, data = runs[runs$mode == "window", ], FUN = mean)
  tw$method <- "monte_carlo explicit window, mean over declared seeds"
  tw$grid <- NA; tw$note <- "mean elapsed per call"; names(tw)[names(tw) == "elapsed_s"] <- "elapsed_s"
  tim <- rbind(tim, tw[, names(tim)])
  utils::write.csv(tim, file.path(out, "timings.csv"), row.names = FALSE)
  cat("\n==== (5) Timings (seconds, this machine, one core; orders of magnitude only) ====\n")
  print(tim[, c("scenario", "method", "B", "grid", "elapsed_s")], row.names = FALSE, digits = 4)

  cat(sprintf("\nTotal elapsed: %.1f s\n", as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
  cat("Output directory:", normalizePath(out, winslash = "/"), "\n")
}

main()
