
#' Sample Size Adjustment for Imperfect Reference Standard
#'
#' Inflates the Buderer sample size to account for misclassification by
#' an imperfect reference standard, using the Rogan-Gladen variance
#' inflation factor \eqn{1 / (Se_{ref} + Sp_{ref} - 1)^2}.
#'
#' @details The parameters \code{d_se} and \code{d_sp} are
#'   \strong{half-widths} of the confidence interval. The Buderer formula
#'   uses these directly; the full CI width is \code{2 * d}.
#'
#'   \strong{How weak a reference standard is too weak.} The inflation
#'   factor is \eqn{VIF = 1 / J_{ref}^2} where \eqn{J_{ref} = Se_{ref} +
#'   Sp_{ref} - 1} is the Youden index of the reference standard. \eqn{VIF}
#'   diverges as \eqn{J_{ref} \to 0}: a reference with
#'   \code{Se_ref = Sp_ref = 0.51} has \eqn{J_{ref} = 0.02} and
#'   \eqn{VIF = 2500}, which multiplies the required sample size by 2500 and
#'   (in versions <= 0.2.0) made the Monte Carlo validation try to allocate
#'   several gigabytes before failing. Two guards now bound this:
#'   \code{min_youden} refuses a reference standard that is barely better
#'   than a coin flip, and \code{max_mc_cells} bounds the memory of the
#'   Monte Carlo validation (which allocates \code{B * N_adjusted} cells).
#'   Both raise an informative error naming the offending quantities rather
#'   than failing obscurely.
#'
#' @param Se Expected sensitivity of index test. Default 0.85.
#' @param Sp Expected specificity of index test. Default 0.90.
#' @param d_se Desired precision for Se (half-width). Default 0.07.
#' @param d_sp Desired precision for Sp (half-width). Default 0.05.
#' @param prev Disease prevalence. Default 0.30.
#' @param Se_ref Sensitivity of reference standard. Default 0.90.
#' @param Sp_ref Specificity of reference standard. Default 0.95.
#' @param loss_rate Expected loss-to-follow-up rate. Default 0.10.
#' @param alpha Significance level. Default 0.05.
#' @param B MC replications for validation. Default 5000. Set to 0 to
#'   skip MC validation. A warning is issued when \code{0 < B < 1000},
#'   since the Monte Carlo error of the validation probabilities may then
#'   be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @param sensitivity_table Logical. If \code{TRUE} (default), compute a
#'   table of inflation factors for varying Se_ref and Sp_ref.
#' @param min_youden Minimum acceptable Youden index of the reference
#'   standard, \eqn{Se_{ref} + Sp_{ref} - 1}. Default 0.5, which caps the
#'   variance inflation factor at 4. A reference standard below this is
#'   too weak for the Rogan-Gladen inflation to give a usable sample size;
#'   the function stops with an informative error. Lower it deliberately
#'   (and expect a very large \code{N}) if you really intend to plan
#'   against such a reference. See \code{Details}.
#' @param max_mc_cells Upper bound on \code{B * N_adjusted}, the number of
#'   cells the Monte Carlo validation allocates. Default \code{2e7}
#'   (roughly a few hundred MB). Exceeding it stops with an informative
#'   error suggesting a smaller \code{B} or \code{B = 0}. This guard also
#'   catches an extreme \code{prev}, which inflates \code{N_adjusted}
#'   without touching the Youden index.
#' @param delta_se Alias for \code{d_se}, matching the \code{delta_*}
#'   naming used elsewhere in the package. If supplied (non-\code{NULL}),
#'   it takes precedence over \code{d_se}. Default \code{NULL} (use
#'   \code{d_se}).
#' @param delta_sp Alias for \code{d_sp}. If supplied (non-\code{NULL}),
#'   it takes precedence over \code{d_sp}. Default \code{NULL} (use
#'   \code{d_sp}).
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{inflation_factor_se}{Variance inflation factor for Se.}
#'     \item{n_diseased_unadjusted}{Buderer n for Se without correction.}
#'     \item{n_diseased_adjusted}{Corrected n for Se.}
#'     \item{n_nondiseased_unadjusted}{Buderer n for Sp without correction.}
#'     \item{n_nondiseased_adjusted}{Corrected n for Sp.}
#'     \item{N_unadjusted}{Total N without correction, i.e. the classical
#'       Buderer total (same value as \code{N_buderer}).}
#'     \item{N_buderer}{Alias of \code{N_unadjusted}: the classical Buderer
#'       total N before the Rogan-Gladen inflation is applied. The
#'       inflation factor actually realised on the total N is
#'       \code{N_adjusted / N_unadjusted} (not exactly \code{VIF}, because
#'       of the separate ceiling/prevalence-weighting steps on the Se and
#'       Sp arms).}
#'     \item{N_adjusted}{Total N with correction.}
#'     \item{N_adjusted_loss}{Total N with correction and loss adjustment.}
#'     \item{sensitivity_table}{Data frame of VIF values (if requested),
#'       with both a \code{VIF} column and an identical \code{inflation_factor}
#'       column (the latter is the more descriptive name; both are kept for
#'       backward compatibility).}
#'     \item{mc_validation}{Data frame (if \code{B > 0}) with, for the
#'       unadjusted and adjusted sample sizes, the probability that the
#'       \strong{apparent} sensitivity CI achieves the target width
#'       (\code{P_width_target}), the mean apparent sensitivity
#'       \eqn{P(index+\,|\,ref+)} (\code{se_apparent}), the true Se
#'       (\code{se_true}) and the verification \code{bias}
#'       (\code{se_apparent - se_true}).}
#'   }
#' @references
#' Staquet M et al. (1981). Methodology for the assessment of new
#' dichotomous diagnostic tests. \emph{J Chronic Dis} 34:599-610.
#' \doi{10.1016/0021-9681(81)90059-X}
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. \emph{Am J Epidemiol} 107:71-76.
#' \doi{10.1093/oxfordjournals.aje.a112510}
#' @note The variance inflation factor \eqn{1/(Se_{ref} + Sp_{ref} - 1)^2}
#'   is the variance multiplier of the Rogan-Gladen misclassification
#'   correction, applied here as an approximation to the loss of precision
#'   incurred when an imperfect reference standard is used. It assumes
#'   conditional independence between the index test and the reference
#'   standard given true disease status; this may be violated when both
#'   measure the same construct.
#'
#'   \strong{Important:} the factor restores the \emph{precision} of the
#'   estimator but does not remove its \emph{bias}. The apparent
#'   sensitivity \eqn{P(index+\,|\,ref+)} estimated against an imperfect
#'   reference converges to a value below the true Se (reported as
#'   \code{bias} in \code{mc_validation}); recovering an unbiased estimate
#'   additionally requires a bias-correction or latent-class analysis at
#'   the analysis stage.
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
                             sensitivity_table = TRUE,
                             min_youden = 0.5,
                             max_mc_cells = 2e7,
                             delta_se = NULL,
                             delta_sp = NULL) {
  # delta_se/delta_sp are aliases of d_se/d_sp, matching the delta_*
  # naming used elsewhere in the package; when supplied they take
  # precedence, but default behaviour (both NULL) is unchanged.
  if (!is.null(delta_se)) d_se <- delta_se
  if (!is.null(delta_sp)) d_sp <- delta_sp

  # --- preserve the caller's RNG state -------------------------------
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      suppressWarnings(rm(".Random.seed", envir = .GlobalEnv)),
      add = TRUE
    )
  }

  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(d_se > 0, d_sp > 0)
  stopifnot(prev > 0, prev < 1)
  stopifnot(Se_ref > 0, Se_ref <= 1, Sp_ref > 0, Sp_ref <= 1)
  stopifnot(Se_ref + Sp_ref > 1)  # must be informative
  stopifnot(loss_rate >= 0, loss_rate < 1)
  stopifnot(min_youden > 0, min_youden <= 1)
  stopifnot(max_mc_cells > 0)
  warn_small_B(B)

  # Youden index of the reference standard and the Rogan-Gladen VIF
  youden_ref <- Se_ref + Sp_ref - 1
  VIF <- 1 / youden_ref^2

  # --- refuse a reference standard so weak that the VIF explodes -----
  if (youden_ref < min_youden) {
    stop(
      sprintf(
        paste0(
          "Reference standard too weak: Youden index = Se_ref + Sp_ref - 1 = ",
          "%.4f (< min_youden = %.4f).\n",
          "  The Rogan-Gladen variance inflation factor 1 / Youden^2 = %.1f ",
          "would multiply the\n  required sample size by %.1f, which is not a ",
          "usable study design.\n",
          "  Use a better reference standard, or lower min_youden ",
          "deliberately if you\n  really intend to plan against this one ",
          "(and expect a very large N)."
        ),
        youden_ref, min_youden, VIF, VIF
      ),
      call. = FALSE
    )
  }

  # Unadjusted sample sizes (Buderer)
  n_unadj_se <- buderer_n(Se, d_se, alpha)
  n_unadj_sp <- buderer_n(Sp, d_sp, alpha)

  # Adjusted sample sizes
  n_adj_se <- ceiling(n_unadj_se * VIF)
  n_adj_sp <- ceiling(n_unadj_sp * VIF)

  # Total N
  N_unadj <- buderer_total_N(Se, Sp, d_se, d_sp, prev, alpha)
  N_adj <- ceiling(max(n_adj_se / prev, n_adj_sp / (1 - prev)))
  N_adj_loss <- ceiling(N_adj / (1 - loss_rate))

  # --- bound the memory of the MC validation -------------------------
  # The validation allocates three B x N_adj matrices. Fail loudly, naming
  # the numbers, instead of letting R abort on a multi-gigabyte allocation.
  if (B > 0 && B * N_adj > max_mc_cells) {
    stop(
      sprintf(
        paste0(
          "Monte Carlo validation would allocate B * N_adjusted = %.0f * %.0f ",
          "= %.3g cells\n  (> max_mc_cells = %.3g), i.e. several matrices of ",
          "that size.\n",
          "  This is driven by VIF = %.2f (Youden_ref = %.3f) and prev = %.3f, ",
          "which give N_adjusted = %.0f.\n",
          "  Set B = 0 to skip the validation, lower B, or raise ",
          "max_mc_cells if you have the memory."
        ),
        B, N_adj, B * N_adj, max_mc_cells, VIF, youden_ref, prev, N_adj
      ),
      call. = FALSE
    )
  }

  # Sensitivity table
  sens_table <- NULL
  if (sensitivity_table) {
    se_ref_grid <- seq(0.80, 1.00, 0.05)
    sp_ref_grid <- seq(0.85, 1.00, 0.05)
    grid <- expand.grid(Se_ref = se_ref_grid, Sp_ref = sp_ref_grid)
    grid <- grid[grid$Se_ref + grid$Sp_ref > 1, , drop = FALSE]
    grid$VIF <- 1 / (grid$Se_ref + grid$Sp_ref - 1)^2
    # inflation_factor is a same-value, more descriptive alias of VIF;
    # VIF is kept so existing code that reads that column still works.
    grid$inflation_factor <- grid$VIF
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

    # Mean of the APPARENT sensitivity P(index+ | ref+). With an imperfect
    # reference this converges to a value below the true Se; the gap is the
    # verification bias that the variance inflation does NOT remove.
    se_apparent_mean_adj <- mean(se_apparent)
    se_apparent_mean_unadj <- mean(se_apparent_unadj)

    mc_validation <- data.frame(
      scenario = c("unadjusted", "adjusted"),
      N = c(N_unadj, N_adj),
      P_width_target = c(P_width_unadj, P_width_adj),
      se_apparent = c(se_apparent_mean_unadj, se_apparent_mean_adj),
      se_true = c(Se, Se),
      bias = c(se_apparent_mean_unadj - Se, se_apparent_mean_adj - Se),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      # this is the Rogan-Gladen variance inflation factor, not the
      # Staquet correction. Named for what it actually is.
      method = paste("Sample Size for Imperfect Reference Standard",
                     "(Rogan-Gladen variance inflation)"),
      n_diseased = n_adj_se,
      n_total = N_adj_loss,
      inflation_factor_se = VIF,
      youden_ref = youden_ref,
      n_diseased_unadjusted = n_unadj_se,
      n_diseased_adjusted = n_adj_se,
      n_nondiseased_unadjusted = n_unadj_sp,
      n_nondiseased_adjusted = n_adj_sp,
      N_unadjusted = N_unadj,
      N_buderer = N_unadj,
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
