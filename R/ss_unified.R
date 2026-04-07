#' Unified Monte Carlo Framework for DTA Sample Size
#'
#' Integrates parameter uncertainty, imperfect reference standards,
#' and multiple accuracy metrics into a single sample size calculation.
#' Simulates the full data-generating process including prevalence
#' uncertainty, reference standard misclassification, and optionally
#' AUC and net benefit constraints.
#'
#' @details The parameters \code{delta_se}, \code{delta_sp}, and
#'   \code{delta_auc} are \strong{half-widths} of the confidence interval.
#'   The target full CI width checked internally is \code{2 * delta}.
#'
#' @param prior_se Beta prior for Se: \code{c(alpha, beta)}.
#'   Default \code{c(17, 3)}.
#' @param prior_sp Beta prior for Sp: \code{c(alpha, beta)}.
#'   Default \code{c(2, 2)}.
#' @param prior_prev Beta prior for prevalence: \code{c(alpha, beta)}.
#'   Default \code{c(5, 15)}.
#' @param Se_ref Reference standard sensitivity. Default 0.92.
#' @param Sp_ref Reference standard specificity. Default 0.95.
#' @param loss_rate Expected losses. Default 0.10.
#' @param delta_se Half-width target for Se. Default 0.07.
#' @param delta_sp Half-width target for Sp. Default 0.05.
#' @param delta_auc Half-width target for AUC. Default 0.06.
#'   Set to 0 to skip AUC check.
#' @param check_nb Logical, check Net Benefit > 0. Default \code{FALSE}.
#' @param pt_range Threshold range for NB if \code{check_nb}.
#'   Default \code{c(0.15, 0.40)}.
#' @param target_assurance Joint assurance target. Default 0.80.
#' @param N_range Search range for total N. Default \code{seq(200, 800, by = 20)}.
#' @param B MC replications. Default 5000.
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{joint_assurance}{Achieved joint assurance at optimal N.}
#'     \item{comparison}{Data frame comparing methods.}
#'   }
#' @export
ss_unified <- function(prior_se = c(17, 3),
                       prior_sp = c(2, 2),
                       prior_prev = c(5, 15),
                       Se_ref = 0.92,
                       Sp_ref = 0.95,
                       loss_rate = 0.10,
                       delta_se = 0.07,
                       delta_sp = 0.05,
                       delta_auc = 0.06,
                       check_nb = FALSE,
                       pt_range = c(0.15, 0.40),
                       target_assurance = 0.80,
                       N_range = seq(200, 1000, by = 20),
                       B = 5000,
                       seed = 2026) {
  # Validate inputs
  stopifnot(length(prior_se) == 2, all(prior_se > 0))
  stopifnot(length(prior_sp) == 2, all(prior_sp > 0))
  stopifnot(length(prior_prev) == 2, all(prior_prev > 0))
  stopifnot(Se_ref > 0, Se_ref <= 1, Sp_ref > 0, Sp_ref <= 1)
  stopifnot(Se_ref + Sp_ref > 1)
  stopifnot(loss_rate >= 0, loss_rate < 1)
  stopifnot(delta_se > 0, delta_sp > 0, delta_auc >= 0)
  stopifnot(B >= 1)

  z <- stats::qnorm(0.975)
  target_se_width <- 2 * delta_se
  target_sp_width <- 2 * delta_sp
  target_auc_width <- 2 * delta_auc

  # Prior means for comparison table
  E_prev <- prior_prev[1] / sum(prior_prev)
  E_se <- prior_se[1] / sum(prior_se)

  optimal_N <- NA_integer_
  joint_assurance_achieved <- NA_real_

  for (N in N_range) {
    set.seed(seed)
    pass_count <- 0L
    valid_count <- 0L

    for (b in seq_len(B)) {
      # Draw parameters from priors
      prev_b <- stats::rbeta(1, prior_prev[1], prior_prev[2])
      Se_true <- stats::rbeta(1, prior_se[1], prior_se[2])
      Sp_true <- stats::rbeta(1, prior_sp[1], prior_sp[2])

      # Generate true disease status
      n_d <- stats::rbinom(1, N, prev_b)
      n_nd <- N - n_d
      if (n_d < 5 || n_nd < 5) next
      valid_count <- valid_count + 1L

      # Generate T x R using conditional independence given D
      # Among n_d truly diseased: T and R independent
      #   P(T+,R+|D+) = Se_true * Se_ref
      #   P(T+,R-|D+) = Se_true * (1-Se_ref)
      #   P(T-,R+|D+) = (1-Se_true) * Se_ref
      #   P(T-,R-|D+) = (1-Se_true) * (1-Se_ref)
      p_d <- c(Se_true * Se_ref,
               Se_true * (1 - Se_ref),
               (1 - Se_true) * Se_ref,
               (1 - Se_true) * (1 - Se_ref))
      cells_d <- stats::rmultinom(1, n_d, p_d)
      # cells_d: [T+R+, T+R-, T-R+, T-R-] among truly diseased

      # Among n_nd truly non-diseased:
      #   P(T+,R+|D-) = (1-Sp_true) * (1-Sp_ref)
      #   P(T+,R-|D-) = (1-Sp_true) * Sp_ref
      #   P(T-,R+|D-) = Sp_true * (1-Sp_ref)
      #   P(T-,R-|D-) = Sp_true * Sp_ref
      p_nd <- c((1 - Sp_true) * (1 - Sp_ref),
                (1 - Sp_true) * Sp_ref,
                Sp_true * (1 - Sp_ref),
                Sp_true * Sp_ref)
      cells_nd <- stats::rmultinom(1, n_nd, p_nd)

      # Observed 2x2 (T x R)
      # a = T+R+ total, b = T+R-, c = T-R+, d = T-R-
      a <- cells_d[1] + cells_nd[1]   # T+ & R+
      b_cell <- cells_d[2] + cells_nd[2]   # T+ & R-
      c_cell <- cells_d[3] + cells_nd[3]   # T- & R+
      d_cell <- cells_d[4] + cells_nd[4]   # T- & R-

      n_ref_pos <- a + c_cell          # R+
      n_ref_neg <- b_cell + d_cell     # R-

      if (n_ref_pos < 2 || n_ref_neg < 2) next

      # Observed Se_obs = P(T+|R+), Sp_obs = P(T-|R-)
      Se_obs <- a / n_ref_pos
      Sp_obs <- d_cell / n_ref_neg

      # Wilson CI for Se_obs
      denom_se <- 1 + z^2 / n_ref_pos
      center_se <- (Se_obs + z^2 / (2 * n_ref_pos)) / denom_se
      margin_se <- z * sqrt((Se_obs * (1 - Se_obs) +
        z^2 / (4 * n_ref_pos)) / n_ref_pos) / denom_se
      width_se <- 2 * margin_se
      se_pass <- (width_se <= target_se_width)

      # Wilson CI for Sp_obs
      denom_sp <- 1 + z^2 / n_ref_neg
      center_sp <- (Sp_obs + z^2 / (2 * n_ref_neg)) / denom_sp
      margin_sp <- z * sqrt((Sp_obs * (1 - Sp_obs) +
        z^2 / (4 * n_ref_neg)) / n_ref_neg) / denom_sp
      width_sp <- 2 * margin_sp
      sp_pass <- (width_sp <= target_sp_width)

      # AUC check (if active)
      auc_pass <- TRUE
      if (delta_auc > 0) {
        # Binormal AUC approximation from observed Se, Sp
        Se_obs_clip <- min(max(Se_obs, 0.01), 0.99)
        Sp_obs_clip <- min(max(Sp_obs, 0.01), 0.99)
        auc_approx <- stats::pnorm(
          (stats::qnorm(Se_obs_clip) + stats::qnorm(Sp_obs_clip)) / sqrt(2)
        )
        var_auc <- hanley_mcneil_var(auc_approx, n_ref_pos, n_ref_neg)
        auc_width <- 2 * z * sqrt(max(var_auc, 0))
        auc_pass <- (auc_width <= target_auc_width)
      }

      # Net Benefit check (if active)
      nb_pass <- TRUE
      if (check_nb) {
        for (pt in pt_range) {
          odds_pt <- pt / (1 - pt)
          nb_hat <- (n_ref_pos / N) * Se_obs -
            (n_ref_neg / N) * (1 - Sp_obs) * odds_pt
          nb_all <- (n_ref_pos / N) - (n_ref_neg / N) * odds_pt
          if (nb_hat <= 0 || nb_hat <= nb_all) {
            nb_pass <- FALSE
            break
          }
        }
      }

      all_pass <- se_pass & sp_pass & auc_pass & nb_pass
      if (all_pass) pass_count <- pass_count + 1L
    }

    if (valid_count > 0) {
      joint_assurance <- pass_count / valid_count
    } else {
      joint_assurance <- 0
    }

    if (joint_assurance >= target_assurance) {
      optimal_N <- as.integer(N)
      joint_assurance_achieved <- joint_assurance
      break
    }
  }

  if (is.na(optimal_N)) {
    warning("No N in N_range achieved target assurance. ",
            "Consider expanding N_range.")
    optimal_N <- max(N_range)
    joint_assurance_achieved <- joint_assurance
  }

  N_enrolled <- ceiling(optimal_N / (1 - loss_rate))

  # Comparison table
  buderer_N_se <- ceiling(buderer_n(E_se, delta_se) / E_prev)
  imperfect_res <- ss_imperfect_ref(
    Se = E_se,
    Sp = prior_sp[1] / sum(prior_sp),
    d_se = delta_se, d_sp = delta_sp,
    prev = E_prev,
    Se_ref = Se_ref, Sp_ref = Sp_ref,
    loss_rate = loss_rate, B = 0,
    sensitivity_table = FALSE
  )

  comparison <- data.frame(
    method = c("Buderer (classical)", "Imperfect ref (Staquet)",
               "Unified (this method)"),
    N = c(buderer_N_se, imperfect_res$N_adjusted_loss, N_enrolled),
    stringsAsFactors = FALSE
  )

  n_d_final <- floor(optimal_N * E_prev)

  structure(
    list(
      method = "Unified Monte Carlo Framework for DTA Sample Size",
      n_diseased = n_d_final,
      n_total = N_enrolled,
      N_effective = optimal_N,
      N_enrolled = N_enrolled,
      joint_assurance = joint_assurance_achieved,
      comparison = comparison,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
