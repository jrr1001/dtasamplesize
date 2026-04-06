#' Joint Sample Size for Sensitivity, Specificity, and AUC
#'
#' Finds the minimum total sample size \code{N} where sensitivity,
#' specificity, and AUC simultaneously achieve target precision with
#' joint probability at least \code{target_prob}. AUC precision is
#' computed deterministically via the Hanley-McNeil variance approximation,
#' while Se and Sp are evaluated via Monte Carlo simulation.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param AUC Expected AUC. Default 0.80.
#' @param delta_se Half-width target for Se. Default 0.07.
#' @param delta_sp Half-width target for Sp. Default 0.05.
#' @param delta_auc Half-width target for AUC. Default 0.05.
#' @param prev Disease prevalence. Default 0.20.
#' @param target_prob Minimum joint probability. Default 0.80.
#' @param N_range Range of total N to search. Default \code{seq(100, 800, by = 10)}.
#' @param B MC replications. Default 5000.
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{n_total}{Minimum total N achieving joint target.}
#'     \item{n_diseased}{Number of diseased at optimal N.}
#'     \item{n_non_diseased}{Number of non-diseased at optimal N.}
#'     \item{joint_prob}{Achieved joint probability at optimal N.}
#'     \item{buderer_N}{Buderer-based total N for comparison.}
#'   }
#' @note AUC CI width is computed using the Hanley-McNeil variance
#'   approximation rather than empirical Mann-Whitney statistics,
#'   making the AUC component deterministic per N. This is an acceptable
#'   approximation for the search phase.
#' @references
#' Hanley JA, McNeil BJ (1982). The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. \emph{Radiology}
#' 143:29-36.
#'
#' Buderer NMF (1996). Statistical methodology: I. Incorporating the
#' prevalence of disease into the sample size calculation for sensitivity
#' and specificity. \emph{Acad Emerg Med} 3:895-900.
#' \doi{10.1111/j.1553-2712.1996.tb03538.x}
#' @examples
#' result <- joint_sample_size(B = 1000, N_range = seq(100, 700, by = 20))
#' print(result)
#' @export
joint_sample_size <- function(Se = 0.85,
                              Sp = 0.90,
                              AUC = 0.80,
                              delta_se = 0.07,
                              delta_sp = 0.05,
                              delta_auc = 0.05,
                              prev = 0.20,
                              target_prob = 0.80,
                              N_range = seq(100, 800, by = 10),
                              B = 5000,
                              seed = 2026) {
  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1, AUC > 0.5, AUC <= 1)
  stopifnot(delta_se > 0, delta_sp > 0, delta_auc > 0)
  stopifnot(prev > 0, prev < 1)
  stopifnot(target_prob > 0, target_prob < 1)
  stopifnot(B >= 1)

  z <- stats::qnorm(1 - 0.05 / 2)  # for Wilson CI at alpha=0.05
  target_se_width <- 2 * delta_se
  target_sp_width <- 2 * delta_sp
  target_auc_width <- 2 * delta_auc

  optimal_N <- NA_integer_
  joint_prob_achieved <- NA_real_

  for (N in N_range) {
    n_d <- floor(N * prev)
    n_nd <- N - n_d
    if (n_d < 2 || n_nd < 2) next

    # --- AUC: deterministic via Hanley-McNeil ---
    var_auc <- hanley_mcneil_var(AUC, n_d, n_nd)
    auc_width <- 2 * stats::qnorm(0.975) * sqrt(var_auc)
    auc_pass <- auc_width <= target_auc_width

    # If AUC alone fails, skip MC for Se/Sp
    if (!auc_pass) next

    # --- Se and Sp: Monte Carlo with Wilson CI ---
    set.seed(seed)
    x_se <- stats::rbinom(B, n_d, Se)
    se_hat <- x_se / n_d

    x_sp <- stats::rbinom(B, n_nd, Sp)
    sp_hat <- x_sp / n_nd

    # Wilson CI width for Se (vectorized)
    denom_se <- 1 + z^2 / n_d
    center_se <- (se_hat + z^2 / (2 * n_d)) / denom_se
    margin_se <- z * sqrt((se_hat * (1 - se_hat) + z^2 / (4 * n_d)) / n_d) / denom_se
    se_lower <- pmax(0, center_se - margin_se)
    se_upper <- pmin(1, center_se + margin_se)
    se_width <- se_upper - se_lower

    # Wilson CI width for Sp (vectorized)
    denom_sp <- 1 + z^2 / n_nd
    center_sp <- (sp_hat + z^2 / (2 * n_nd)) / denom_sp
    margin_sp <- z * sqrt((sp_hat * (1 - sp_hat) + z^2 / (4 * n_nd)) / n_nd) / denom_sp
    sp_lower <- pmax(0, center_sp - margin_sp)
    sp_upper <- pmin(1, center_sp + margin_sp)
    sp_width <- sp_upper - sp_lower

    # Joint: all three pass
    pass <- (se_width <= target_se_width) &
      (sp_width <= target_sp_width)
    # AUC pass is scalar TRUE at this point
    joint_prob <- mean(pass)

    if (joint_prob >= target_prob) {
      optimal_N <- as.integer(N)
      joint_prob_achieved <- joint_prob
      break
    }
  }

  if (is.na(optimal_N)) {
    warning("No N in N_range achieved target joint probability. ",
            "Consider expanding N_range.")
    optimal_N <- max(N_range)
    joint_prob_achieved <- joint_prob  # last computed value
  }

  n_d_final <- floor(optimal_N * prev)
  n_nd_final <- optimal_N - n_d_final

  # Buderer comparison
  buderer_N <- ceiling(max(
    buderer_n(Se, delta_se) / prev,
    buderer_n(Sp, delta_sp) / (1 - prev)
  ))

  structure(
    list(
      method = "Joint Sample Size for Se + Sp + AUC",
      n_total = optimal_N,
      n_diseased = n_d_final,
      n_non_diseased = n_nd_final,
      joint_prob = joint_prob_achieved,
      buderer_N = buderer_N,
      Se = Se,
      Sp = Sp,
      AUC = AUC,
      prev = prev,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
