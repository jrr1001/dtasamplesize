#' Sample Size Adjustment for Imperfect Reference Standard
#'
#' Inflates the Buderer sample size to account for misclassification by
#' an imperfect reference standard, using the variance inflation factor
#' derived from the Staquet correction.
#'
#' @param Se Expected sensitivity of index test. Default 0.85.
#' @param Sp Expected specificity of index test. Default 0.90.
#' @param d_se Desired precision for Se. Default 0.07.
#' @param d_sp Desired precision for Sp. Default 0.05.
#' @param prev Disease prevalence. Default 0.30.
#' @param Se_ref Sensitivity of reference standard. Default 0.90.
#' @param Sp_ref Specificity of reference standard. Default 0.95.
#' @param loss_rate Expected loss-to-follow-up rate. Default 0.10.
#' @param alpha Significance level. Default 0.05.
#' @param B MC replications for validation. Default 5000. Set to 0 to
#'   skip MC validation.
#' @param seed Random seed. Default 2026.
#' @param sensitivity_table Logical. If \code{TRUE} (default), compute a
#'   table of inflation factors for varying Se_ref and Sp_ref.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{inflation_factor_se}{Variance inflation factor for Se.}
#'     \item{n_diseased_unadjusted}{Buderer n for Se without correction.}
#'     \item{n_diseased_adjusted}{Corrected n for Se.}
#'     \item{n_nondiseased_unadjusted}{Buderer n for Sp without correction.}
#'     \item{n_nondiseased_adjusted}{Corrected n for Sp.}
#'     \item{N_unadjusted}{Total N without correction.}
#'     \item{N_adjusted}{Total N with correction.}
#'     \item{N_adjusted_loss}{Total N with correction and loss adjustment.}
#'     \item{sensitivity_table}{Data frame of VIF values (if requested).}
#'     \item{mc_validation}{MC validation results (if B > 0).}
#'   }
#' @references
#' Staquet M et al. (1981). Methodology for the assessment of new
#' dichotomous diagnostic tests. \emph{J Chronic Dis} 34:599-610.
#' \doi{10.1016/0021-9681(81)90059-X}
#' @note The Staquet correction assumes conditional independence between
#'   index test and reference standard given true disease status. This
#'   assumption may be violated when both measure the same construct.
#' @examples
#' result <- ss_imperfect_ref(B = 0)
#' print(result)
#' @export
ss_imperfect_ref <- function(Se = 0.85,
                             Sp = 0.90,
                             d_se = 0.07,
                             d_sp = 0.05,
                             prev = 0.30,
                             Se_ref = 0.90,
                             Sp_ref = 0.95,
                             loss_rate = 0.10,
                             alpha = 0.05,
                             B = 5000,
                             seed = 2026,
                             sensitivity_table = TRUE) {
  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(d_se > 0, d_sp > 0)
  stopifnot(prev > 0, prev < 1)
  stopifnot(Se_ref > 0, Se_ref <= 1, Sp_ref > 0, Sp_ref <= 1)
  stopifnot(Se_ref + Sp_ref > 1)  # must be informative
  stopifnot(loss_rate >= 0, loss_rate < 1)

  # Variance inflation factor
  VIF <- 1 / (Se_ref + Sp_ref - 1)^2

  # Unadjusted sample sizes (Buderer)
  n_unadj_se <- buderer_n(Se, d_se, alpha)
  n_unadj_sp <- buderer_n(Sp, d_sp, alpha)

  # Adjusted sample sizes
  n_adj_se <- ceiling(n_unadj_se * VIF)
  n_adj_sp <- ceiling(n_unadj_sp * VIF)

  # Total N
  N_unadj <- ceiling(max(n_unadj_se / prev, n_unadj_sp / (1 - prev)))
  N_adj <- ceiling(max(n_adj_se / prev, n_adj_sp / (1 - prev)))
  N_adj_loss <- ceiling(N_adj / (1 - loss_rate))

  # Sensitivity table
  sens_table <- NULL
  if (sensitivity_table) {
    se_ref_grid <- seq(0.80, 1.00, 0.05)
    sp_ref_grid <- seq(0.85, 1.00, 0.05)
    grid <- expand.grid(Se_ref = se_ref_grid, Sp_ref = sp_ref_grid)
    grid <- grid[grid$Se_ref + grid$Sp_ref > 1, , drop = FALSE]
    grid$VIF <- 1 / (grid$Se_ref + grid$Sp_ref - 1)^2
    grid$n_adj_se <- ceiling(n_unadj_se * grid$VIF)
    grid$N_adj <- ceiling(pmax(grid$n_adj_se / prev,
                               ceiling(n_unadj_sp * grid$VIF) / (1 - prev)))
    sens_table <- grid
  }

  # MC validation
  mc_validation <- NULL
  if (B > 0) {
    set.seed(seed)

    # Generate true disease status
    true_disease <- stats::rbinom(B * N_adj, 1, prev)

    # For each subject, generate reference test result
    ref_result <- ifelse(
      true_disease == 1,
      stats::rbinom(length(true_disease), 1, Se_ref),
      1L - stats::rbinom(length(true_disease), 1, Sp_ref)
    )

    # Generate index test result
    index_result <- ifelse(
      true_disease == 1,
      stats::rbinom(length(true_disease), 1, Se),
      1L - stats::rbinom(length(true_disease), 1, Sp)
    )

    # Reshape into B datasets of size N_adj
    dim(true_disease) <- c(N_adj, B)
    dim(ref_result) <- c(N_adj, B)
    dim(index_result) <- c(N_adj, B)

    z_val <- stats::qnorm(1 - alpha / 2)
    target_width_se <- 2 * d_se

    # For adjusted n: use all N_adj subjects
    # Apparent Se = P(index+ | ref+)
    ref_pos <- colSums(ref_result == 1)
    index_pos_given_ref_pos <- colSums(index_result == 1 & ref_result == 1)
    se_apparent <- index_pos_given_ref_pos / pmax(ref_pos, 1)
    se_se_adj <- sqrt(se_apparent * (1 - se_apparent) / pmax(ref_pos, 1))
    width_adj <- 2 * z_val * se_se_adj
    P_width_adj <- mean(width_adj <= target_width_se)

    # For unadjusted n: use only first N_unadj subjects
    n_use <- min(N_unadj, N_adj)
    ref_pos_unadj <- colSums(ref_result[seq_len(n_use), , drop = FALSE] == 1)
    idx_pos_unadj <- colSums(
      index_result[seq_len(n_use), , drop = FALSE] == 1 &
        ref_result[seq_len(n_use), , drop = FALSE] == 1
    )
    se_apparent_unadj <- idx_pos_unadj / pmax(ref_pos_unadj, 1)
    se_se_unadj <- sqrt(se_apparent_unadj * (1 - se_apparent_unadj) /
                          pmax(ref_pos_unadj, 1))
    width_unadj <- 2 * z_val * se_se_unadj
    P_width_unadj <- mean(width_unadj <= target_width_se)

    mc_validation <- data.frame(
      scenario = c("unadjusted", "adjusted"),
      N = c(N_unadj, N_adj),
      P_width_target = c(P_width_unadj, P_width_adj),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      method = "Sample Size for Imperfect Reference Standard (Staquet Correction)",
      n_diseased = n_adj_se,
      n_total = N_adj_loss,
      inflation_factor_se = VIF,
      n_diseased_unadjusted = n_unadj_se,
      n_diseased_adjusted = n_adj_se,
      n_nondiseased_unadjusted = n_unadj_sp,
      n_nondiseased_adjusted = n_adj_sp,
      N_unadjusted = N_unadj,
      N_adjusted = N_adj,
      N_adjusted_loss = N_adj_loss,
      Se_ref = Se_ref,
      Sp_ref = Sp_ref,
      VIF = VIF,
      sensitivity_table = sens_table,
      mc_validation = mc_validation,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
