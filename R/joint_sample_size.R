
#' Joint Sample Size for Sensitivity, Specificity, and AUC
#'
#' Finds the minimum total sample size \code{N} at which sensitivity and
#' specificity simultaneously achieve their target precision with joint
#' probability at least \code{target_prob}, subject to the AUC also
#' achieving its target precision. Se and Sp are evaluated by Monte Carlo;
#' AUC precision is computed deterministically via the Hanley-McNeil
#' variance approximation. \strong{Read the \code{Details}: the AUC
#' component is not probabilistic.}
#'
#' @details The parameters \code{delta_se}, \code{delta_sp}, and
#'   \code{delta_auc} are \strong{half-widths} of the confidence interval.
#'   The target full CI width checked internally is \code{2 * delta}.
#'
#'   \strong{What \code{joint_prob_se_sp} is, and what it is not.} The
#'   returned probability is the Monte Carlo probability that the Se
#'   \strong{and} Sp confidence intervals both meet their target width. It
#'   is a joint probability over \strong{two} of the three criteria, not
#'   three. The AUC criterion is applied as a \strong{deterministic gate}: at
#'   each candidate \code{N} the Hanley-McNeil variance is evaluated at the
#'   \emph{assumed} \code{AUC} and the candidate is discarded unless the
#'   resulting CI width meets \code{2 * delta_auc}. Because that width is
#'   computed at the assumed parameter rather than simulated, the AUC
#'   component carries no sampling variability: it is a
#'   \strong{median-style criterion}, satisfied in roughly 50\% of real
#'   studies (the ones whose realised AUC standard error lands below its
#'   expectation), \strong{not} an assurance of \code{target_prob}. The
#'   overall probability that all three intervals meet their targets in a
#'   real study is therefore \strong{lower} than the reported
#'   \code{joint_prob_se_sp}. Treat \code{joint_prob_se_sp} as the assurance
#'   for the Se/Sp pair, conditional on an AUC criterion that is met on
#'   average.
#'
#'   \strong{Geometric coherence of \code{AUC} with \code{Se} and \code{Sp}.}
#'   \code{Se}, \code{Sp} and \code{AUC} are not free parameters: any ROC
#'   curve that is concave and passes through the operating point
#'   \eqn{(1 - Sp,\ Se)} has area at least that of the polygon joining
#'   \eqn{(0,0)}, \eqn{(1 - Sp,\ Se)} and \eqn{(1,1)}, namely
#'   \deqn{AUC_{min} = \tfrac{1}{2}(1 - Sp)\,Se + \tfrac{1}{2} Sp\,(Se + 1).}
#'   Supplying an \code{AUC} below \eqn{AUC_{min}} describes a test that
#'   cannot exist, and the function now stops with an error rather than
#'   returning a sample size for it. At the defaults \code{Se = 0.85} and
#'   \code{Sp = 0.90} this bound is \strong{0.875}; the previous default of
#'   \code{AUC = 0.80} was therefore impossible. For reference, the binormal
#'   AUC implied by the same operating point,
#'   \eqn{\Phi[(\Phi^{-1}(Se) + \Phi^{-1}(Sp)) / \sqrt{2}]}, is 0.949.
#'   The default \code{AUC = 0.90} sits between the two: it is geometrically
#'   attainable, and being smaller than the binormal value it yields a
#'   larger Hanley-McNeil variance and hence a conservative sample size.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param AUC Expected AUC. Default 0.90. Must be geometrically compatible
#'   with \code{Se} and \code{Sp}: see \code{Details}. Raising \code{Se} or
#'   \code{Sp} raises the minimum attainable AUC, so a high-accuracy
#'   operating point may require \code{AUC} to be supplied explicitly.
#' @param delta_se Half-width target for Se. Default 0.07.
#' @param delta_sp Half-width target for Sp. Default 0.05.
#' @param delta_auc Half-width target for AUC. Default 0.05.
#' @param prev Disease prevalence. Default 0.20.
#' @param target_prob Minimum joint probability \strong{for Se and Sp}.
#'   Default 0.80. See \code{Details}.
#' @param N_range Range of total N to search. Default \code{seq(100, 800, by = 10)}.
#' @param B MC replications. Default 5000. A warning is issued when
#'   \code{0 < B < 1000}, since the Monte Carlo error of the reported joint
#'   probability may then be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{n_total}{Minimum total N achieving the joint target.}
#'     \item{n_diseased}{Number of diseased at optimal N.}
#'     \item{n_non_diseased}{Number of non-diseased at optimal N.}
#'     \item{joint_prob_se_sp}{Monte Carlo probability that the \strong{Se
#'       and Sp} intervals both meet their target width at the optimal N.
#'       This is \strong{not} a three-way joint probability: the AUC
#'       criterion is a deterministic median-style gate (see
#'       \code{Details}). \code{NA} when the AUC gate was never passed
#'       anywhere in \code{N_range}, in which case no Se/Sp probability was
#'       ever computed.}
#'     \item{joint_prob}{\strong{Deprecated} alias of
#'       \code{joint_prob_se_sp}, kept for backward compatibility with
#'       version 0.2.0. Use \code{joint_prob_se_sp}.}
#'     \item{auc_gate_passed}{Logical: whether any \code{N} in
#'       \code{N_range} satisfied the AUC precision criterion.}
#'     \item{AUC_min}{The geometric lower bound on AUC implied by
#'       \code{Se} and \code{Sp}.}
#'     \item{buderer_N}{Buderer-based total N for comparison.}
#'   }
#' @note AUC CI width is computed using the Hanley-McNeil variance
#'   approximation evaluated at the assumed AUC, rather than from simulated
#'   Mann-Whitney statistics. This makes the AUC component deterministic per
#'   \code{N} and therefore a median-style, not an assurance-style,
#'   criterion. See \code{Details}.
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
#'
#' # An AUC below the geometric minimum is refused:
#' try(joint_sample_size(Se = 0.85, Sp = 0.90, AUC = 0.80))
#' @export
joint_sample_size <- function(Se = 0.85,
                              Sp = 0.90,
                              AUC = 0.90,
                              delta_se = 0.07,
                              delta_sp = 0.05,
                              delta_auc = 0.05,
                              prev = 0.20,
                              target_prob = 0.80,
                              N_range = seq(100, 800, by = 10),
                              B = 5000,
                              seed = 2026) {
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
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1, AUC > 0.5, AUC <= 1)
  stopifnot(delta_se > 0, delta_sp > 0, delta_auc > 0)
  stopifnot(prev > 0, prev < 1)
  stopifnot(target_prob > 0, target_prob < 1)
  stopifnot(B >= 1)
  warn_small_B(B)

  # --- (a): AUC must be geometrically compatible with (Se, Sp) ---------
  # Minimum area of a concave ROC through the operating point (1 - Sp, Se):
  # the polygon (0,0) -> (1-Sp, Se) -> (1,1).
  AUC_min <- 0.5 * (1 - Sp) * Se + 0.5 * Sp * (Se + 1)
  if (AUC < AUC_min - 1e-9) {
    AUC_binormal <- stats::pnorm(
      (stats::qnorm(Se) + stats::qnorm(Sp)) / sqrt(2)
    )
    stop(
      sprintf(
        paste0(
          "AUC = %.3f is geometrically impossible for Se = %.3f and Sp = %.3f.\n",
          "  Any concave ROC curve through the operating point (1 - Sp, Se) = ",
          "(%.3f, %.3f)\n  has area at least AUC_min = %.4f ",
          "[= 0.5*(1-Sp)*Se + 0.5*Sp*(Se+1)].\n",
          "  For reference, the binormal AUC implied by this operating point ",
          "is %.4f.\n",
          "  Supply an AUC in [%.4f, 1], or lower Se / Sp."
        ),
        AUC, Se, Sp, 1 - Sp, Se, AUC_min, AUC_binormal, AUC_min
      ),
      call. = FALSE
    )
  }

  target_se_width <- 2 * delta_se
  target_sp_width <- 2 * delta_sp
  target_auc_width <- 2 * delta_auc

  optimal_N <- NA_integer_
  joint_prob_achieved <- NA_real_
  # (c): distinguish "never computed" from "computed and equal to 0".
  joint_prob <- NA_real_
  auc_gate_passed <- FALSE

  for (N in N_range) {
    n_d <- floor(N * prev)
    n_nd <- N - n_d
    if (n_d < 2 || n_nd < 2) next

    # --- AUC: DETERMINISTIC gate via Hanley-McNeil (median-style) ---
    var_auc <- hanley_mcneil_var(AUC, n_d, n_nd)
    auc_width <- 2 * stats::qnorm(0.975) * sqrt(var_auc)
    auc_pass <- auc_width <= target_auc_width

    # If AUC alone fails, skip MC for Se/Sp
    if (!auc_pass) next
    auc_gate_passed <- TRUE

    # --- Se and Sp: Monte Carlo with Wilson CI ---
    set.seed(seed)
    x_se <- stats::rbinom(B, n_d, Se)
    x_sp <- stats::rbinom(B, n_nd, Sp)

    se_width <- wilson_width(x_se, n_d)
    sp_width <- wilson_width(x_sp, n_nd)

    # Joint over Se and Sp ONLY (the AUC gate is deterministic, see @details)
    pass <- (se_width <= target_se_width) &
      (sp_width <= target_sp_width)
    joint_prob <- mean(pass)

    if (joint_prob >= target_prob) {
      optimal_N <- as.integer(N)
      joint_prob_achieved <- joint_prob
      break
    }
  }

  if (is.na(optimal_N)) {
    if (!auc_gate_passed) {
      # (c): the AUC gate blocked every candidate N, so the Se/Sp Monte
      # Carlo probability was NEVER computed. Reporting 0 here (the old
      # initialisation value) would be indistinguishable from a genuine 0.
      warning("The AUC precision target (delta_auc = ", delta_auc,
              ") was not met at any N in N_range, so no N was ever ",
              "evaluated for Se/Sp. joint_prob_se_sp is reported as NA ",
              "(not 0: it was never computed). Expand N_range or relax ",
              "delta_auc.")
      joint_prob_achieved <- NA_real_
    } else {
      warning("No N in N_range achieved the target joint probability for ",
              "Se and Sp. Consider expanding N_range.")
      joint_prob_achieved <- joint_prob  # last computed value
    }
    optimal_N <- max(N_range)
  }

  n_d_final <- floor(optimal_N * prev)
  n_nd_final <- optimal_N - n_d_final

  # Buderer comparison
  buderer_N <- buderer_total_N(Se, Sp, delta_se, delta_sp, prev)

  structure(
    list(
      method = "Joint Sample Size for Se + Sp + AUC",
      n_total = optimal_N,
      n_diseased = n_d_final,
      n_non_diseased = n_nd_final,
      joint_prob_se_sp = joint_prob_achieved,
      # Deprecated alias, kept so that 0.2.0 code keeps running.
      joint_prob = joint_prob_achieved,
      auc_gate_passed = auc_gate_passed,
      buderer_N = buderer_N,
      Se = Se,
      Sp = Sp,
      AUC = AUC,
      AUC_min = AUC_min,
      prev = prev,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
