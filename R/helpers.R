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
