#' Maximum Allowed Evaluation Budget for the Net-Benefit Ceiling
#'
#' Hard upper bound on \code{ss_unified}'s \code{nb_B_ceiling} argument (and
#' on \code{nb_assurance_ceiling}'s own \code{B_ceiling}), enforced with a
#' clear \code{stop()} rather than silently attempting the computation. This
#' guards against a pathologically large request (a typo, or a value copied
#' from an unrelated context) consuming excessive compute time; see
#' \code{\link{nb_assurance_ceiling}} for why the default sits far below
#' this cap and for the memory-exhaustion incident this cap and the lower
#' default were both introduced to prevent.
#' @keywords internal
#' @noRd
NB_B_CEILING_MAX <- 20000000L  # 2e7

#' Deterministic Gauss-Legendre Quadrature Nodes and Weights on [0, 1]
#'
#' Nodes and weights for \code{n}-point Gauss-Legendre quadrature on the
#' unit interval, computed by the Golub-Welsch algorithm: the nodes are the
#' eigenvalues of the symmetric tridiagonal Jacobi matrix for the Legendre
#' three-term recurrence, and the weights are \code{2 * (first component of
#' each normalized eigenvector)^2} (Golub & Welsch 1969); both are then
#' mapped from the canonical interval \eqn{[-1, 1]} to \eqn{[0, 1]}. Used by
#' \code{\link{nb_assurance_ceiling}} to integrate over the Se and Sp priors
#' deterministically, in place of Monte Carlo sampling.
#'
#' @param n Number of nodes (positive integer).
#' @return List with elements \code{x} (nodes in \eqn{(0, 1)}, strictly
#'   interior and ascending -- so a Beta density that diverges at 0 or 1,
#'   e.g. under a shape parameter below 1, is never evaluated exactly at the
#'   singularity) and \code{w} (weights, summing to 1).
#' @references
#' Golub GH, Welsch JH (1969). Calculation of Gauss quadrature rules.
#' \emph{Math Comp} 23:221-230. \doi{10.1090/S0025-5718-69-99647-1}
#' @keywords internal
#' @noRd
gauss_legendre_unit <- function(n) {
  if (n < 1L) {
    stop("gauss_legendre_unit(): n must be a positive integer.", call. = FALSE)
  }
  if (n == 1L) return(list(x = 0.5, w = 1))
  i <- seq_len(n - 1L)
  off_diag <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- off_diag
  J[cbind(i + 1L, i)] <- off_diag
  eig <- eigen(J, symmetric = TRUE)
  ord <- order(eig$values)
  x <- eig$values[ord]
  w <- 2 * (eig$vectors[1, ord])^2
  list(x = (x + 1) / 2, w = w / 2)
}

#' Wilson Score Confidence Interval for a Proportion
#'
#' Computes the Wilson score confidence interval for a binomial proportion.
#'
#' @param x Number of successes.
#' @param n Number of trials.
#' @param alpha Significance level (default 0.05).
#' @return Named numeric vector with elements \code{lower}, \code{upper},
#'   and \code{width}.
#' @references
#' Wilson EB (1927). Probable inference, the law of succession, and
#' statistical inference. \emph{J Am Stat Assoc} 22:209-212.
#' @keywords internal
#' @export
wilson_ci <- function(x, n, alpha = 0.05) {
  p_hat <- x / n
  z <- stats::qnorm(1 - alpha / 2)
  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
  lower <- pmax(0, center - margin)
  upper <- pmin(1, center + margin)
  c(lower = lower, upper = upper, width = upper - lower)
}

#' Wald Confidence Interval for a Proportion
#'
#' Computes the standard Wald confidence interval for a binomial proportion.
#'
#' @param x Number of successes.
#' @param n Number of trials.
#' @param alpha Significance level (default 0.05).
#' @return Named numeric vector with elements \code{lower}, \code{upper},
#'   and \code{width}.
#' @keywords internal
#' @export
wald_ci <- function(x, n, alpha = 0.05) {
  p_hat <- x / n
  se <- sqrt(p_hat * (1 - p_hat) / n)
  z <- stats::qnorm(1 - alpha / 2)
  lower <- pmax(0, p_hat - z * se)
  upper <- pmin(1, p_hat + z * se)
  c(lower = lower, upper = upper, width = upper - lower)
}

#' Hanley-McNeil Variance of AUC
#'
#' Computes the variance of the area under the ROC curve using the
#' Hanley-McNeil approximation.
#'
#' @param auc Estimated AUC.
#' @param n_cases Number of diseased individuals.
#' @param n_controls Number of non-diseased individuals.
#' @return Numeric variance estimate.
#' @references
#' Hanley JA, McNeil BJ (1982). The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. \emph{Radiology}
#' 143:29-36.
#' @keywords internal
#' @export
hanley_mcneil_var <- function(auc, n_cases, n_controls) {
  Q1 <- auc / (2 - auc)
  Q2 <- 2 * auc^2 / (1 + auc)
  var_auc <- (auc * (1 - auc) +
    (n_cases - 1) * (Q1 - auc^2) +
    (n_controls - 1) * (Q2 - auc^2)) / (n_cases * n_controls)
  var_auc
}

#' Buderer Sample Size for a Single Proportion
#'
#' Computes the required number of diseased (or non-diseased) individuals
#' to estimate sensitivity (or specificity) with a given precision using
#' the Buderer formula.
#'
#' @param p Expected proportion (sensitivity or specificity).
#' @param d Desired precision (half-width of the confidence interval).
#' @param alpha Significance level (default 0.05).
#' @return Integer sample size (ceiling).
#' @references
#' Buderer NMF (1996). Statistical methodology: I. Incorporating the
#' prevalence of disease into the sample size calculation for sensitivity
#' and specificity. \emph{Acad Emerg Med} 3:895-900.
#' \doi{10.1111/j.1553-2712.1996.tb03538.x}
#' @keywords internal
#' @export
buderer_n <- function(p, d, alpha = 0.05) {
  z <- stats::qnorm(1 - alpha / 2)
  n <- ceiling(z^2 * p * (1 - p) / d^2)
  n
}

#' Buderer Total Sample Size Across Both Arms
#'
#' Total sample size that simultaneously supplies enough diseased subjects
#' for the sensitivity target and enough non-diseased subjects for the
#' specificity target, using the Buderer per-proportion formula on each arm.
#' Equal to
#' \code{ceiling(max(buderer_n(Se, d_se, alpha) / prev,
#' buderer_n(Sp, d_sp, alpha) / (1 - prev)))}.
#'
#' @param Se Expected sensitivity.
#' @param Sp Expected specificity.
#' @param d_se Desired precision (half-width) for Se.
#' @param d_sp Desired precision (half-width) for Sp.
#' @param prev Disease prevalence.
#' @param alpha Significance level (default 0.05).
#' @return Integer total sample size (ceiling).
#' @keywords internal
#' @noRd
buderer_total_N <- function(Se, Sp, d_se, d_sp, prev, alpha = 0.05) {
  ceiling(max(buderer_n(Se, d_se, alpha) / prev,
              buderer_n(Sp, d_sp, alpha) / (1 - prev)))
}

#' Warn When B Is Small
#'
#' Issues a warning when the number of Monte Carlo replications \code{B} is
#' small enough that the Monte Carlo error of the reported assurance /
#' probability may be substantial. The warning is controlled by the
#' \code{dtasamplesize.warn_small_B} option (default \code{TRUE}), so that
#' scripts, tests, and examples that deliberately use a small \code{B} for
#' speed can silence it with
#' \code{options(dtasamplesize.warn_small_B = FALSE)}.
#'
#' @param B Number of Monte Carlo replications.
#' @return \code{invisible(NULL)}. Called for its warning side effect.
#' @keywords internal
#' @noRd
warn_small_B <- function(B) {
  if (B > 0 && B < 1000 &&
        isTRUE(getOption("dtasamplesize.warn_small_B", TRUE))) {
    warning(
      "B = ", B, " is small; Monte Carlo error in the assurance estimate ",
      "may be substantial. B >= 1000 (ideally 5000) is recommended for ",
      "final results. Set options(dtasamplesize.warn_small_B = FALSE) to ",
      "silence this.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Save the Caller's RNG State (Kind and Seed)
#'
#' Captures both the active RNG kind (\code{RNGkind()}: the uniform, normal,
#' and discrete-sampling generators) and the current \code{.Random.seed}
#' value (or its absence), so that a function which temporarily reseeds the
#' generator can restore the caller's \strong{exact} state afterwards --
#' including the kind, not just the seed vector. Restoring only
#' \code{.Random.seed} is not enough: \code{set.seed()} called later by
#' unrelated code with no explicit \code{kind} argument reuses whichever
#' kind is \strong{currently active}, not whatever kind \code{.Random.seed}
#' happens to encode, so a caller left on a non-default generator after a
#' function returns would silently draw from that generator instead of
#' their own. This is the same mechanism \code{\link{ss_unified}} already
#' uses (it needs L'Ecuyer-CMRG internally); \code{save_rng_state()} and
#' \code{\link{restore_rng_state}} generalize it for use by any function
#' in this package that seeds the generator.
#'
#' @return A list with elements \code{kind} (the \code{RNGkind()} vector),
#'   \code{seed_present} (logical), and \code{seed} (the
#'   \code{.Random.seed} value, or \code{NULL} when absent), for use with
#'   \code{\link{restore_rng_state}}.
#' @keywords internal
#' @noRd
save_rng_state <- function() {
  kind <- RNGkind()
  seed_present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  seed <- if (seed_present) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  list(kind = kind, seed_present = seed_present, seed = seed)
}

#' Restore an RNG State Saved by \code{save_rng_state()}
#'
#' Meant to be called from \code{on.exit(restore_rng_state(state), add =
#' TRUE)}, immediately after \code{state <- save_rng_state()}, so that the
#' caller's generator kind and \code{.Random.seed} are put back exactly as
#' found, including when the calling function exits via an error.
#'
#' @param state Object returned by \code{\link{save_rng_state}}.
#' @return \code{invisible(NULL)}. Called for its side effect.
#' @keywords internal
#' @noRd
restore_rng_state <- function(state) {
  # Setting the kind (even back to what it already was) can itself create
  # .Random.seed with an arbitrary value when none existed before, so the
  # exact seed vector -- or its absence -- is restored AFTER the kind,
  # overwriting whatever RNGkind() just did.
  suppressWarnings(RNGkind(
    kind = state$kind[1],
    normal.kind = state$kind[2],
    sample.kind = state$kind[3]
  ))
  if (state$seed_present) {
    assign(".Random.seed", state$seed, envir = .GlobalEnv)
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }
  invisible(NULL)
}

#' One-Sided Wilson Score Lower Confidence Bound for a Proportion
#'
#' Computes the lower limit of a one-sided Wilson score confidence interval
#' for a proportion \code{phat} estimated from an effective sample size
#' \code{n} (\code{n} need not be a literal count of Bernoulli trials --
#' \code{\link{ss_unified}}'s \code{decision = "isotonic"} rule uses it with
#' the pooled effective sample size of a smoothed block of grid points,
#' since the quantity being bounded there is a block-pooled proportion, not
#' a raw per-N count).
#'
#' Unlike the ordinary normal-approximation (Wald) bound
#' \code{phat - z * sqrt(phat * (1 - phat) / n)}, the Wilson bound does not
#' collapse to a vacuous value whenever \code{phat} sits at 0 or 1: at
#' \code{phat = 1} the Wald standard error is \strong{exactly} 0 regardless
#' of \code{n}, so the Wald lower bound is exactly \code{1} even from a
#' single successful replicate (\code{n = 1}) -- reporting total certainty
#' from the smallest possible amount of evidence, an artifact of the normal
#' approximation to the estimate rather than a real statistical conclusion.
#' The Wilson interval instead inverts the normal approximation to the
#' \emph{score} (the standardized difference between \code{phat} and the
#' hypothesized proportion), which keeps its bound strictly inside
#' \eqn{(0, 1)} whenever \code{n} is finite, even at the extremes. It
#' converges to the same Wald bound as \code{n -> Inf} -- the two differ by
#' \eqn{O(z^2/n)}, on the order of \eqn{10^{-4}} already at the \code{B =
#' 20000} used for this package's own published results and therefore
#' immaterial to any of them -- so it can replace the Wald formula
#' everywhere a lower confidence bound on a Monte Carlo proportion is
#' needed without moving results computed at realistic \code{B}.
#'
#' @param phat Estimated proportion (e.g. \code{pass_count / B}, or a
#'   pooled/fitted proportion), in \eqn{[0, 1]}.
#' @param n Effective sample size behind \code{phat} (may be a vector,
#'   recycled against \code{phat}).
#' @param z Critical value for the desired one-sided confidence level
#'   (e.g. \code{qnorm(0.95)} for a one-sided 95\% lower bound).
#' @return Numeric (vector), the Wilson lower confidence limit(s), each in
#'   \eqn{[0, 1]}.
#' @references
#' Wilson EB (1927). Probable inference, the law of succession, and
#' statistical inference. \emph{J Am Stat Assoc} 22:209-212.
#' @keywords internal
#' @noRd
wilson_lower <- function(phat, n, z) {
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half_width <- z * sqrt(pmax(phat * (1 - phat) / n + z^2 / (4 * n^2), 0)) / denom
  pmax(0, pmin(1, center - half_width))
}

#' Width of the Wilson Score Interval
#'
#' Convenience wrapper returning only the width of the Wilson score
#' interval for a proportion. Vectorized in \code{x}. Byte-identical to the
#' \code{width} element of \code{\link{wilson_ci}} (same clamping to
#' \eqn{[0, 1]}, same order of operations).
#'
#' @param x Number of successes (scalar or vector).
#' @param n Number of trials.
#' @param alpha Significance level (default 0.05).
#' @return Numeric Wilson interval width(s).
#' @keywords internal
#' @noRd
wilson_width <- function(x, n, alpha = 0.05) {
  p_hat <- x / n
  z <- stats::qnorm(1 - alpha / 2)
  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n)) / n) / denom
  pmin(1, center + margin) - pmax(0, center - margin)
}
