#' Sample Size for Net Benefit (Decision Curve Analysis)
#'
#' Finds the minimum total sample size such that a study can
#' \emph{conclude}, with \code{(1 - alpha)} confidence, that the index
#' test has positive net benefit and outperforms the default strategies
#' (treat-all and treat-none) across a range of threshold probabilities.
#'
#' @details For each threshold probability \code{pt}, the net benefit of
#'   the test and of the treat-all strategy are
#'   \deqn{NB = prev \cdot Se - (1 - prev)(1 - Sp)\,\frac{pt}{1 - pt},}
#'   \deqn{NB_{all} = prev - (1 - prev)\,\frac{pt}{1 - pt}.}
#'   The treat-none strategy has \eqn{NB = 0} by definition. A test is
#'   useful at \code{pt} only when \eqn{NB > 0} \strong{and}
#'   \eqn{NB > NB_{all}}.
#'
#'   The sample size criterion is \strong{inference-based}: for each Monte
#'   Carlo study the lower limit of the two-sided \code{(1 - alpha)}
#'   confidence interval for \eqn{NB} (normal approximation, with
#'   \deqn{Var(NB) = (n_d/N)^2 Se(1-Se)/n_d + (n_{nd}/N)^2 (pt/(1-pt))^2 Sp(1-Sp)/n_{nd})}
#'   is required to exceed both 0 and \eqn{NB_{all}}. The required \code{N}
#'   is the smallest value for which this happens with probability at least
#'   \code{target_prob} ("assurance"). This differs from earlier versions,
#'   which used the point estimate of \eqn{NB} and therefore returned the
#'   smallest \code{N} whenever the test was clearly useful.
#'
#'   When the test is \strong{not} useful at a threshold under the assumed
#'   parameters (\eqn{NB \le 0} or \eqn{NB \le NB_{all}}), no sample size
#'   can demonstrate superiority; the threshold is flagged
#'   \code{feasible = FALSE} and \code{N_required = NA}.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param prev Disease prevalence. Default 0.20.
#' @param pt_range Threshold probabilities to evaluate.
#'   Default \code{seq(0.10, 0.50, 0.05)}.
#' @param target_prob Assurance: probability that the \code{(1 - alpha)}
#'   confidence interval for net benefit excludes both 0 and the treat-all
#'   net benefit. Default 0.80.
#' @param alpha Confidence-interval significance level. Default 0.05.
#' @param N_range Range of total N to search.
#'   Default \code{seq(50, 1500, by = 10)}.
#' @param B MC replications. Default 5000.
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{N_by_pt}{Data frame with columns \code{pt}, \code{feasible},
#'       \code{N_required}, \code{prob_achieved}, \code{NB_true},
#'       \code{NB_treat_all}.}
#'     \item{N_conservative}{Maximum required N across feasible thresholds
#'       (\code{NA} if none feasible).}
#'   }
#' @note The net benefit sampling distribution is approximated by the
#'   normal distribution; at very small \code{n_d} or \code{n_nd} this
#'   approximation degrades. Thresholds at which the test is not useful
#'   under the assumed parameters are reported as infeasible rather than
#'   assigned an (unachievable) sample size.
#' @references
#' Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method
#' for evaluating prediction models. \emph{Med Decis Making} 26:565-574.
#' \doi{10.1177/0272989X06295361}
#'
#' Vickers AJ, van Calster B, Steyerberg EW (2019). A simple,
#' step-by-step guide to interpreting decision curve analysis.
#' \emph{Diagn Progn Res} 3:18. \doi{10.1186/s41512-019-0064-7}
#' @examples
#' result <- ss_net_benefit(pt_range = c(0.10, 0.30), B = 1000)
#' print(result)
#' @export
ss_net_benefit <- function(Se = 0.85,
                           Sp = 0.90,
                           prev = 0.20,
                           pt_range = seq(0.10, 0.50, 0.05),
                           target_prob = 0.80,
                           alpha = 0.05,
                           N_range = seq(50, 1500, by = 10),
                           B = 5000,
                           seed = 2026) {
  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(prev > 0, prev < 1)
  stopifnot(all(pt_range > 0 & pt_range < 1))
  stopifnot(target_prob > 0, target_prob < 1)
  stopifnot(alpha > 0, alpha < 1)
  stopifnot(B >= 1)

  z <- stats::qnorm(1 - alpha / 2)
  results_list <- vector("list", length(pt_range))

  for (i in seq_along(pt_range)) {
    pt <- pt_range[i]
    odds_pt <- pt / (1 - pt)
    NB_true <- prev * Se - (1 - prev) * (1 - Sp) * odds_pt
    NB_treat_all <- prev - (1 - prev) * odds_pt

    # A test can only be shown superior if it is genuinely useful:
    # NB > 0 (beats treat-none) and NB > NB_treat_all (beats treat-all).
    feasible <- (NB_true > 0) && (NB_true > NB_treat_all)

    found_N <- NA_integer_
    found_prob <- NA_real_
    best_prob <- 0

    if (feasible) {
      for (N in N_range) {
        n_d <- floor(N * prev)
        n_nd <- N - n_d
        if (n_d < 2 || n_nd < 2) next

        set.seed(seed)
        TP <- stats::rbinom(B, n_d, Se)
        FP <- stats::rbinom(B, n_nd, 1 - Sp)
        Se_hat <- TP / n_d
        Sp_hat <- 1 - FP / n_nd

        w_d <- n_d / N
        w_nd <- n_nd / N
        NB_hat <- w_d * Se_hat - w_nd * (1 - Sp_hat) * odds_pt

        # Variance of NB_hat (n_d, n_nd fixed given N; TP, FP independent)
        var_nb <- w_d^2 * Se_hat * (1 - Se_hat) / n_d +
          w_nd^2 * odds_pt^2 * Sp_hat * (1 - Sp_hat) / n_nd
        se_nb <- sqrt(pmax(var_nb, 0))
        nb_lower <- NB_hat - z * se_nb

        # Study is "successful" if the CI lower limit clears both
        # competing strategies (treat-none = 0, treat-all).
        success <- (nb_lower > 0) & (nb_lower > NB_treat_all)
        prob <- mean(success)
        best_prob <- max(best_prob, prob)

        if (prob >= target_prob) {
          found_N <- as.integer(N)
          found_prob <- prob
          break
        }
      }

      if (is.na(found_N)) {
        warning("Threshold pt = ", pt, ": target assurance not reached ",
                "within N_range (best = ", round(best_prob, 3),
                "). Consider expanding N_range.")
        found_prob <- best_prob
      }
    }

    results_list[[i]] <- data.frame(
      pt = pt,
      feasible = feasible,
      N_required = found_N,
      prob_achieved = if (feasible) found_prob else NA_real_,
      NB_true = NB_true,
      NB_treat_all = NB_treat_all,
      stringsAsFactors = FALSE
    )
  }

  N_by_pt <- do.call(rbind, results_list)
  rownames(N_by_pt) <- NULL

  # N_conservative is the worst-case requirement over *feasible* thresholds.
  # If a feasible threshold could not be met within N_range, the true
  # worst case exceeds the search range and is unknown, so report NA
  # rather than the (smaller) maximum over the thresholds that did converge.
  feas <- N_by_pt$feasible
  unmet_feasible <- any(feas & is.na(N_by_pt$N_required))
  achievable <- N_by_pt$N_required[feas & !is.na(N_by_pt$N_required)]
  if (unmet_feasible) {
    warning("At least one feasible threshold did not reach target ",
            "assurance within N_range; N_conservative is reported as NA. ",
            "Expand N_range for a finite worst-case sample size.")
    N_conservative <- NA_integer_
  } else if (length(achievable) > 0) {
    N_conservative <- max(achievable)
  } else {
    N_conservative <- NA_integer_
  }

  # n_diseased / n_total derived from the conservative N (if any)
  n_d_final <- if (is.na(N_conservative)) NA_integer_ else floor(N_conservative * prev)

  structure(
    list(
      method = "Sample Size for Net Benefit (Decision Curve Analysis)",
      n_diseased = n_d_final,
      n_total = N_conservative,
      N_by_pt = N_by_pt,
      N_conservative = N_conservative,
      Se = Se,
      Sp = Sp,
      prev = prev,
      alpha = alpha,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
