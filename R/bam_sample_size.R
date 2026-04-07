#' Bayesian Assurance Method for DTA Sample Size
#'
#' Finds the minimum number of diseased (and non-diseased) individuals
#' such that the probability of the posterior credible interval achieving
#' the target width is at least \code{target_assurance}.
#'
#' @details The parameters \code{delta_se} and \code{delta_sp} are the
#'   \strong{full width} of the posterior credible interval (not half-width).
#'   For example, \code{delta_se = 0.14} corresponds to a half-width of 0.07.
#'   This differs from other functions in this package that use half-width.
#'
#' @param prior_se Numeric vector \code{c(alpha, beta)} for Beta prior on
#'   sensitivity. Default \code{c(17, 3)} (E[Se]=0.85, ~20 pseudo-observations).
#' @param prior_sp Numeric vector \code{c(alpha, beta)} for Beta prior on
#'   specificity. Default \code{c(2, 2)} (vague).
#' @param delta_se Target full width for Se credible interval. Default 0.14.
#' @param delta_sp Target full width for Sp credible interval. Default 0.10.
#' @param target_assurance Minimum assurance probability. Default 0.80.
#' @param prior_prev Numeric vector \code{c(alpha, beta)} for Beta prior on
#'   prevalence. Default \code{c(6, 14)} (E[prev]=0.30).
#' @param n_range Integer vector of candidate n values to search.
#'   Default \code{20:500}.
#' @param B Number of MC replications. Default 5000.
#' @param alpha_ci Credible interval level. Default 0.95.
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{n_diseased}{Minimum diseased sample for Se assurance.}
#'     \item{n_non_diseased}{Minimum non-diseased sample for Sp assurance.}
#'     \item{assurance_se}{Achieved assurance for Se at \code{n_diseased}.}
#'     \item{assurance_sp}{Achieved assurance for Sp at \code{n_non_diseased}.}
#'     \item{N_total_median}{Median total N accounting for prevalence uncertainty.}
#'     \item{N_total_P75}{75th percentile of total N.}
#'     \item{N_total_P90}{90th percentile of total N.}
#'     \item{buderer_n_se}{Buderer sample size for comparison.}
#'   }
#' @references
#' Wilson KJ et al. (2022). Bayesian sample size determination for
#' diagnostic accuracy studies. \emph{Stat Med} 41:2908-2922.
#' \doi{10.1002/sim.9393}
#' @examples
#' result <- bam_sample_size(B = 500, n_range = 20:500)
#' print(result)
#' @export
bam_sample_size <- function(prior_se = c(17, 3),
                            prior_sp = c(2, 2),
                            delta_se = 0.14,
                            delta_sp = 0.10,
                            target_assurance = 0.80,
                            prior_prev = c(6, 14),
                            n_range = 20:500,
                            B = 5000,
                            alpha_ci = 0.95,
                            seed = 2026) {
  # Validate inputs
  stopifnot(length(prior_se) == 2, all(prior_se > 0))
  stopifnot(length(prior_sp) == 2, all(prior_sp > 0))
  stopifnot(length(prior_prev) == 2, all(prior_prev > 0))
  stopifnot(delta_se > 0, delta_sp > 0)
  stopifnot(target_assurance > 0, target_assurance < 1)
  stopifnot(B >= 1)

  a_se <- prior_se[1]
  b_se <- prior_se[2]
  a_sp <- prior_sp[1]
  b_sp <- prior_sp[2]

  ci_lower_q <- (1 - alpha_ci) / 2
  ci_upper_q <- 1 - ci_lower_q

  # --- Search for n_se ---
  n_se <- NA_integer_
  assurance_se_achieved <- NA_real_

  for (n in n_range) {
    set.seed(seed)
    se_true <- stats::rbeta(B, a_se, b_se)
    x <- stats::rbinom(B, n, se_true)
    # Posterior: Beta(a_se + x, b_se + n - x)
    post_a <- a_se + x
    post_b <- b_se + n - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance <- mean(ci_width <= delta_se)
    if (assurance >= target_assurance) {
      n_se <- as.integer(n)
      assurance_se_achieved <- assurance
      break
    }
  }

  if (is.na(n_se)) {
    warning("No n in n_range achieved target assurance for Se. ",
            "Consider expanding n_range.")
    n_se <- max(n_range)
    # Compute assurance at max
    set.seed(seed)
    se_true <- stats::rbeta(B, a_se, b_se)
    x <- stats::rbinom(B, n_se, se_true)
    post_a <- a_se + x
    post_b <- b_se + n_se - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance_se_achieved <- mean(ci_width <= delta_se)
  }

  # --- Search for n_sp ---
  n_sp <- NA_integer_
  assurance_sp_achieved <- NA_real_

  for (n in n_range) {
    set.seed(seed)
    sp_true <- stats::rbeta(B, a_sp, b_sp)
    x <- stats::rbinom(B, n, sp_true)
    post_a <- a_sp + x
    post_b <- b_sp + n - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance <- mean(ci_width <= delta_sp)
    if (assurance >= target_assurance) {
      n_sp <- as.integer(n)
      assurance_sp_achieved <- assurance
      break
    }
  }

  if (is.na(n_sp)) {
    warning("No n in n_range achieved target assurance for Sp. ",
            "Consider expanding n_range.")
    n_sp <- max(n_range)
    set.seed(seed)
    sp_true <- stats::rbeta(B, a_sp, b_sp)
    x <- stats::rbinom(B, n_sp, sp_true)
    post_a <- a_sp + x
    post_b <- b_sp + n_sp - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance_sp_achieved <- mean(ci_width <= delta_sp)
  }

  # --- Total N distribution accounting for prevalence uncertainty ---
  set.seed(seed)
  prev_draws <- stats::rbeta(B, prior_prev[1], prior_prev[2])
  N_total_draws <- n_se / prev_draws
  N_total_median <- stats::median(N_total_draws)
  N_total_P75 <- stats::quantile(N_total_draws, 0.75, names = FALSE)
  N_total_P90 <- stats::quantile(N_total_draws, 0.90, names = FALSE)

  # --- Buderer comparison ---
  buderer_n_se <- buderer_n(a_se / (a_se + b_se), delta_se / 2)

  structure(
    list(
      method = "Bayesian Assurance Method (BAM) for DTA Sample Size",
      n_diseased = n_se,
      n_non_diseased = n_sp,
      n_total = ceiling(N_total_median),
      assurance_se = assurance_se_achieved,
      assurance_sp = assurance_sp_achieved,
      N_total_median = ceiling(N_total_median),
      N_total_P75 = ceiling(N_total_P75),
      N_total_P90 = ceiling(N_total_P90),
      buderer_n_se = buderer_n_se,
      prior_se = prior_se,
      prior_sp = prior_sp,
      prior_prev = prior_prev,
      delta_se = delta_se,
      delta_sp = delta_sp,
      target_assurance = target_assurance,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
