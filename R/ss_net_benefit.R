#' Sample Size for Net Benefit (Decision Curve Analysis)
#'
#' Finds the minimum total sample size such that the estimated net benefit
#' exceeds zero and the treat-all net benefit with probability at least
#' \code{target_prob}, across a range of threshold probabilities.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param prev Disease prevalence. Default 0.20.
#' @param pt_range Threshold probabilities to evaluate.
#'   Default \code{seq(0.10, 0.50, 0.05)}.
#' @param target_prob Probability that NB > 0 and NB > treat-all at each
#'   threshold. Default 0.80.
#' @param N_range Range of total N to search.
#'   Default \code{seq(100, 600, by = 10)}.
#' @param B MC replications. Default 5000.
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{N_by_pt}{Data frame with columns \code{pt}, \code{N_required},
#'       \code{prob_achieved}, \code{NB_true}.}
#'     \item{N_conservative}{Maximum N across all thresholds.}
#'   }
#' @references
#' Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method
#' for evaluating prediction models. \emph{Med Decis Making} 26:565-574.
#' \doi{10.1177/0272989X06295361}
#' @examples
#' result <- ss_net_benefit(pt_range = c(0.10, 0.30), B = 1000)
#' print(result)
#' @export
ss_net_benefit <- function(Se = 0.85,
                           Sp = 0.90,
                           prev = 0.20,
                           pt_range = seq(0.10, 0.50, 0.05),
                           target_prob = 0.80,
                           N_range = seq(100, 600, by = 10),
                           B = 5000,
                           seed = 2026) {
  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(prev > 0, prev < 1)
  stopifnot(all(pt_range > 0 & pt_range < 1))
  stopifnot(target_prob > 0, target_prob < 1)
  stopifnot(B >= 1)

  results_list <- vector("list", length(pt_range))

  for (i in seq_along(pt_range)) {
    pt <- pt_range[i]
    odds_pt <- pt / (1 - pt)
    NB_true <- prev * Se - (1 - prev) * (1 - Sp) * odds_pt
    NB_treat_all <- prev - (1 - prev) * odds_pt

    found_N <- NA_integer_
    found_prob <- NA_real_

    for (N in N_range) {
      n_d <- floor(N * prev)
      n_nd <- N - n_d
      if (n_d < 2 || n_nd < 2) next

      set.seed(seed)
      TP <- stats::rbinom(B, n_d, Se)
      FP <- stats::rbinom(B, n_nd, 1 - Sp)
      Se_hat <- TP / n_d
      Sp_hat <- 1 - FP / n_nd

      NB_hat <- (n_d / N) * Se_hat - (n_nd / N) * (1 - Sp_hat) * odds_pt
      NB_all_hat <- (n_d / N) - (n_nd / N) * odds_pt

      success <- (NB_hat > 0) & (NB_hat > NB_all_hat)
      prob <- mean(success)

      if (prob >= target_prob) {
        found_N <- as.integer(N)
        found_prob <- prob
        break
      }
    }

    if (is.na(found_N)) {
      found_N <- max(N_range)
      found_prob <- prob
    }

    results_list[[i]] <- data.frame(
      pt = pt,
      N_required = found_N,
      prob_achieved = found_prob,
      NB_true = NB_true,
      NB_treat_all = NB_treat_all,
      stringsAsFactors = FALSE
    )
  }

  N_by_pt <- do.call(rbind, results_list)
  rownames(N_by_pt) <- NULL
  N_conservative <- max(N_by_pt$N_required)

  # Use the max N for the n_diseased/n_total outputs
  n_d_final <- floor(N_conservative * prev)

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
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
