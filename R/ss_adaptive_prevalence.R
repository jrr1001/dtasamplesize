#' Adaptive Sample Size with Blinded Prevalence Re-estimation
#'
#' Simulates a two-stage adaptive design where prevalence is estimated
#' at an interim stage and the sample size is adjusted upward if the
#' observed prevalence is lower than initially assumed.
#'
#' @details The parameters \code{d_se} and \code{d_sp} are
#'   \strong{half-widths} of the confidence interval. The full CI width
#'   target is \code{2 * d}.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param d_se Precision for Se (half-width). Default 0.07.
#' @param d_sp Precision for Sp (half-width). Default 0.05.
#' @param prev_initial Initially assumed prevalence. Default 0.30.
#' @param prev_true_range Numeric vector of true prevalence scenarios.
#'   Default \code{c(0.18, 0.25, 0.30, 0.35, 0.42)}.
#' @param fraction_stage1 Fraction of initial N for stage 1. Default 0.40.
#' @param loss_rate Expected loss-to-follow-up rate. Default 0.10.
#' @param B MC replications. Default 5000.
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{N_initial}{Initial planned N (before adaptation).}
#'     \item{n_stage1}{Number of subjects in stage 1.}
#'     \item{results}{Data frame with one row per true prevalence scenario,
#'       showing N_final mean/median/P75, precision achieved.}
#'   }
#' @references
#' Stark M, Zapf A (2020). Sample size calculation and re-estimation for
#' diagnostic accuracy studies based on sensitivity and specificity.
#' \emph{Stat Methods Med Res} 29:2958-2971.
#' \doi{10.1177/0962280220913588}
#' @examples
#' result <- ss_adaptive_prevalence(B = 500, prev_true_range = c(0.20, 0.30))
#' print(result)
#' @export
ss_adaptive_prevalence <- function(Se = 0.85,
                                   Sp = 0.90,
                                   d_se = 0.07,
                                   d_sp = 0.05,
                                   prev_initial = 0.30,
                                   prev_true_range = c(0.18, 0.25, 0.30, 0.35, 0.42),
                                   fraction_stage1 = 0.40,
                                   loss_rate = 0.10,
                                   B = 5000,
                                   seed = 2026) {
  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(d_se > 0, d_sp > 0)
  stopifnot(prev_initial > 0, prev_initial < 1)
  stopifnot(all(prev_true_range > 0 & prev_true_range < 1))
  stopifnot(fraction_stage1 > 0, fraction_stage1 < 1)
  stopifnot(loss_rate >= 0, loss_rate < 1)
  stopifnot(B >= 1)

  # Required sample sizes for Se and Sp (Buderer)
  n_se <- buderer_n(Se, d_se)
  n_sp <- buderer_n(Sp, d_sp)

  # Initial N
  N_initial <- ceiling(max(n_se / prev_initial, n_sp / (1 - prev_initial)))
  N_initial_adj <- ceiling(N_initial / (1 - loss_rate))

  # Stage 1 sample size
  n_stage1 <- ceiling(fraction_stage1 * N_initial)

  z <- stats::qnorm(1 - 0.05 / 2)
  target_se_width <- 2 * d_se
  target_sp_width <- 2 * d_sp

  results_list <- vector("list", length(prev_true_range))

  for (i in seq_along(prev_true_range)) {
    prev_true <- prev_true_range[i]
    set.seed(seed)

    N_finals <- numeric(B)
    precision_met <- logical(B)

    for (b in seq_len(B)) {
      # Stage 1: estimate prevalence
      D_stage1 <- stats::rbinom(1, n_stage1, prev_true)
      prev_hat <- max(D_stage1 / n_stage1, 0.05)

      # Recalculate N
      N_new <- max(n_se / prev_hat, n_sp / (1 - prev_hat))
      N_final <- max(ceiling(N_new), N_initial)
      N_final_adj <- ceiling(N_final / (1 - loss_rate))

      N_finals[b] <- N_final_adj

      # Simulate full study
      n_d <- stats::rbinom(1, N_final_adj, prev_true)
      n_nd <- N_final_adj - n_d
      n_d <- max(n_d, 1)
      n_nd <- max(n_nd, 1)

      se_obs <- stats::rbinom(1, n_d, Se) / n_d
      sp_obs <- stats::rbinom(1, n_nd, Sp) / n_nd

      # Wilson CI for Se
      denom_se <- 1 + z^2 / n_d
      center_se <- (se_obs + z^2 / (2 * n_d)) / denom_se
      margin_se <- z * sqrt((se_obs * (1 - se_obs) +
        z^2 / (4 * n_d)) / n_d) / denom_se
      width_se <- 2 * margin_se

      # Wilson CI for Sp
      denom_sp <- 1 + z^2 / n_nd
      center_sp <- (sp_obs + z^2 / (2 * n_nd)) / denom_sp
      margin_sp <- z * sqrt((sp_obs * (1 - sp_obs) +
        z^2 / (4 * n_nd)) / n_nd) / denom_sp
      width_sp <- 2 * margin_sp

      precision_met[b] <- (width_se <= target_se_width) &
        (width_sp <= target_sp_width)
    }

    results_list[[i]] <- data.frame(
      prev_true = prev_true,
      N_final_mean = mean(N_finals),
      N_final_median = stats::median(N_finals),
      N_final_P75 = stats::quantile(N_finals, 0.75, names = FALSE),
      precision_achieved = mean(precision_met),
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  structure(
    list(
      method = "Adaptive Sample Size with Prevalence Re-estimation",
      n_diseased = ceiling(N_initial_adj * prev_initial),
      n_total = N_initial_adj,
      N_initial = N_initial,
      N_initial_adj = N_initial_adj,
      n_stage1 = n_stage1,
      n_se = n_se,
      n_sp = n_sp,
      results = results,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
