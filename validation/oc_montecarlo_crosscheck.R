#!/usr/bin/env Rscript
# =============================================================================
# oc_montecarlo_crosscheck.R
#
# Monte Carlo cross-check of a subset of the exact operating characteristics
# computed by validation/operating_characteristics.R (success probability and
# coverage at fixed true values) and validation/misspecified_priors.R
# (detection probability of the prior-data conflict check).
#
# The script shares no code with those two scripts. It simulates the study
# (K ~ Binomial(N, prevalence); x ~ Binomial(K, Se); y ~ Binomial(N - K, Sp)),
# computes each interval with the analysis method of its design and, for the
# conflict check, calls bam_prior_data_conflict() of the package directly.
#
# Conventions (the same as in operating_characteristics.R)
#   * Success = full width of the Se interval <= 0.14 and full width of the
#     Sp interval <= 0.10, both with tolerance TOL = 1e-12.
#   * An empty arm (no diseased or no non-diseased participants) is a failure.
#   * Two-sided 95% intervals. Wald: truncated to [0, 1] (x = 0 or x = n gives
#     width 0). Wilson: score interval. Bayesian: equal-tailed, analysis prior
#     equal to the design prior.
#   * Coverage of Se is computed over the replicates with at least one
#     diseased participant (coverage of Sp likewise, with at least one
#     non-diseased participant).
#   * Detection: the two-sided flag of bam_prior_data_conflict() (default
#     alpha = 0.05), per quantity and for "any quantity".
#   * Criterion: |z| <= 3 in each cell. For success and coverage z uses the
#     binomial Monte Carlo standard error (MCSE) of the simulated proportion
#     (floored at 1e-4 for success); for detection it uses sqrt(p (1 - p) / R)
#     with p the exact value. The criterion is evaluated by
#     make_manuscript_tables_R1.R (anchor S3.9f-z).
#
# Cells
#   Success and coverage: R = 20000 replicates per cell, seeds 11 to 18.
#   Detection: R = 5000 replicates per scenario, seeds 31 to 33.
#   The exact values to compare with are those printed by
#   operating_characteristics.R and misspecified_priors.R (success, coverage
#   of Se and detection), at four decimals.
#
# Usage
#   Rscript validation/oc_montecarlo_crosscheck.R <output_dir>
#   Writes ONLY into <output_dir>: oc_montecarlo_crosscheck.log,
#   oc_crosscheck_success_coverage.csv and oc_crosscheck_conflict_detection.csv.
#   Libraries: if the environment variables DTA_PKG_LIB and DTA_EXTRA_LIB
#   exist they are prepended to .libPaths() (in that order). Single core; the
#   detection part takes a few minutes.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) stop("Usage: Rscript oc_montecarlo_crosscheck.R <output_dir>")
out <- args[1L]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

paths <- c(Sys.getenv("DTA_PKG_LIB", unset = ""), Sys.getenv("DTA_EXTRA_LIB", unset = ""))
paths <- paths[nzchar(paths)]
if (length(paths)) .libPaths(c(paths, .libPaths()))
suppressPackageStartupMessages(library(dtasamplesize))

logf <- file.path(out, "oc_montecarlo_crosscheck.log")
sink(logf, split = TRUE)
cat("oc_montecarlo_crosscheck.R\n")
cat("Date: ", format(Sys.time()), "; dtasamplesize ", as.character(packageVersion("dtasamplesize")),
    " (", find.package("dtasamplesize"), "); ", R.version.string, "\n", sep = "")

TOL <- 1e-12
z <- qnorm(0.975)
wald_ci <- function(x, n) {
  p <- x / n; h <- z * sqrt(p * (1 - p) / n)
  cbind(pmax(0, p - h), pmin(1, p + h))
}
wilson_ci <- function(x, n) {
  p <- x / n; d <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / d
  h <- z / d * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  cbind(pmax(0, c0 - h), pmin(1, c0 + h))
}
bayes_ci <- function(x, n, ab) {
  cbind(qbeta(0.025, ab[1] + x, ab[2] + n - x), qbeta(0.975, ab[1] + x, ab[2] + n - x))
}
sim <- function(N, se, sp, prev, R, seed) {
  set.seed(seed)
  K <- rbinom(R, N, prev); x <- rbinom(R, K, se); y <- rbinom(R, N - K, sp)
  data.frame(K = K, x = x, y = y, M = N - K)
}
arm <- function(method, x, n, ab = NULL) {
  ci <- switch(method, wald = wald_ci(x, n), wilson = wilson_ci(x, n), bayes = bayes_ci(x, n, ab))
  ci[n == 0, ] <- NA
  ci
}
cell <- function(label, design, N, se, sp, prev, method, ab_se = NULL, ab_sp = NULL, R = 20000, seed = 1) {
  s <- sim(N, se, sp, prev, R, seed)
  cs <- arm(method, s$x, s$K, ab_se); cp <- arm(method, s$y, s$M, ab_sp)
  ws <- cs[, 2] - cs[, 1]; wp <- cp[, 2] - cp[, 1]
  ok <- !is.na(ws) & !is.na(wp) & ws <= 0.14 + TOL & wp <= 0.10 + TOL
  covse <- (cs[, 1] <= se & se <= cs[, 2])[s$K > 0]
  covsp <- (cp[, 1] <= sp & sp <= cp[, 2])[s$M > 0]
  data.frame(label = label, design = design, N = N, Se = se, Sp = sp, prev = prev, R = R,
             success = mean(ok), success_mcse = sqrt(mean(ok) * (1 - mean(ok)) / R),
             cov_se = mean(covse), cov_se_mcse = sqrt(mean(covse) * (1 - mean(covse)) / length(covse)),
             cov_sp = mean(covsp), cov_sp_mcse = sqrt(mean(covsp) * (1 - mean(covsp)) / length(covsp)))
}

# Priors of the worked example: P1 (vague Sp prior, N = 678) and P2 (Sp prior
# Beta(18, 2), N = 672).
P1 <- list(se = c(17, 3), sp = c(2, 2), prev = c(4, 16))
P2 <- list(se = c(17, 3), sp = c(18, 2), prev = c(4, 16))

res <- rbind(
  cell("planning point", "Buderer Wald", 500, 0.85, 0.90, 0.20, "wald", seed = 11),
  cell("planning point", "Joint Wilson", 580, 0.85, 0.90, 0.20, "wilson", seed = 12),
  cell("planning point", "BAM P1", 678, 0.85, 0.90, 0.20, "bayes", P1$se, P1$sp, seed = 13),
  cell("planning point", "BAM P2", 672, 0.85, 0.90, 0.20, "bayes", P2$se, P2$sp, seed = 14),
  cell("Se 0.75", "BAM P1", 678, 0.75, 0.90, 0.20, "bayes", P1$se, P1$sp, seed = 15),
  cell("Se 0.80", "Joint Wilson", 580, 0.80, 0.90, 0.20, "wilson", seed = 16),
  cell("prevalence 0.10", "Buderer Wald", 500, 0.85, 0.90, 0.10, "wald", seed = 17),
  cell("Se 0.65 (S02)", "BAM P1", 678, 0.65, 0.90, 0.20, "bayes", P1$se, P1$sp, seed = 18))

# Exact values (four decimals) of operating_characteristics.R (rows 1 to 7)
# and misspecified_priors.R (row 8, scenario S02 of P1).
exact <- c(0.5265, 0.8070, 0.9994, 0.9992, 0.8349, 0.3572, 0.0561, 0.1649)
exact_cov_se <- c(0.9364, 0.9511, 0.9621, 0.9622, NA, NA, NA, NA)
res$exact_success <- exact
res$z_success <- (res$success - exact) / pmax(res$success_mcse, 1e-4)
res$exact_cov_se <- exact_cov_se
res$z_cov_se <- (res$cov_se - exact_cov_se) / res$cov_se_mcse
cat("\n==== Success probability and coverage (R = 20000 per cell) ====\n")
print(res, digits = 4, row.names = FALSE)

# Detection of the prior-data conflict check: two-sided rule of the package
# (flag_two_sided), per quantity and for any quantity.
det <- function(label, pr, N, se, sp, prev, R = 5000, seed = 21) {
  s <- sim(N, se, sp, prev, R, seed)
  f <- t(vapply(seq_len(R), function(i) {
    r <- bam_prior_data_conflict(tp = s$x[i], fn = s$K[i] - s$x[i], tn = s$y[i], fp = s$M[i] - s$y[i],
                                 prior_se = pr$se, prior_sp = pr$sp, prior_prev = pr$prev)
    g <- function(q) { v <- r$flag_two_sided[grepl(q, r$quantity)]; if (length(v)) any(v) else FALSE }
    c(se = g("sensitivity"), sp = g("specificity"), prev = g("prevalence"))
  }, logical(3)))
  m <- colMeans(f); anyf <- mean(apply(f, 1, any))
  data.frame(label = label, R = R, det_se = m["se"], det_sp = m["sp"], det_prev = m["prev"], det_any = anyf,
             mcse_max = sqrt(0.25 / R), row.names = NULL)
}
t0 <- Sys.time()
dres <- rbind(det("P1 S00 (0.85, 0.90, 0.20)", P1, 678, 0.85, 0.90, 0.20, seed = 31),
              det("P1 S02 (0.65, 0.90, 0.20)", P1, 678, 0.65, 0.90, 0.20, seed = 32),
              det("P2 S04 (0.85, 0.70, 0.20)", P2, 672, 0.85, 0.70, 0.20, seed = 33))
# Exact detection probabilities (four decimals) of misspecified_priors.R
# (two-sided flag with the analysis prior of each design).
dres$exact_se <- c(0.0000, 0.5210, 0.0000)
dres$exact_sp <- c(0.2679, 0.2679, 0.9663)
dres$exact_any <- c(0.2679, 0.6494, 0.9663)
# z = (Monte Carlo - exact) / sqrt(exact (1 - exact) / R); NA when the exact value is 0 or 1
zf <- function(m, e, R) ifelse(e > 0.001 & e < 0.999, (m - e) / sqrt(e * (1 - e) / R), NA_real_)
dres$z_se <- zf(dres$det_se, dres$exact_se, dres$R)
dres$z_sp <- zf(dres$det_sp, dres$exact_sp, dres$R)
dres$z_any <- zf(dres$det_any, dres$exact_any, dres$R)
cat("\n==== Detection of the conflict check (R = 5000 per scenario) ====\n")
print(dres, digits = 4, row.names = FALSE)
cat("Detection time: ", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), " s\n", sep = "")
write.csv(res, file.path(out, "oc_crosscheck_success_coverage.csv"), row.names = FALSE)
write.csv(dres, file.path(out, "oc_crosscheck_conflict_detection.csv"), row.names = FALSE)
cat("Criterion: |z| <= 3 in each cell (binomial MCSE as declared).\n")
sink()
