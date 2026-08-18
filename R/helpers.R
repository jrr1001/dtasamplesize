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
