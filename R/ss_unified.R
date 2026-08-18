
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
#'   \strong{The assurance is unconditional.} Each of the \code{B}
#'   replications draws \code{prev}, \code{Se} and \code{Sp} from their
#'   priors, then simulates the full \eqn{T \times R} table. A replication
#'   that comes out degenerate -- fewer than 5 truly diseased or
#'   non-diseased subjects, or fewer than 2 reference-positive or
#'   reference-negative subjects -- \strong{counts as a failure}, because a
#'   study that lands in that state does not deliver the target precision.
#'   The denominator of the reported assurance is therefore always \code{B}.
#'   In versions <= 0.2.0 degenerate replications were dropped from the
#'   denominator, which reported an assurance \emph{conditional} on the
#'   study not being degenerate -- not the quantity a planner needs.
#'
#'   \strong{The net-benefit criterion (\code{check_nb}) is
#'   inference-based.} At each threshold in \code{pt_range} (threshold odds
#'   \eqn{w = pt / (1 - pt)}), the replication must show, from its own
#'   simulated data, that the test beats \emph{both} default strategies with
#'   \code{(1 - 0.05)} confidence: the lower confidence limit of
#'   \eqn{\widehat{NB}} must exceed 0, \strong{and} the lower confidence
#'   limit of the difference \eqn{\widehat{D} = \widehat{NB} -
#'   \widehat{NB}_{all}} must exceed 0. Because \code{N} is fixed by design
#'   but the disease (here, reference-positive) status is random, the
#'   relevant variances are the multinomial ones used by
#'   \code{\link{ss_net_benefit}} with \code{design = "cohort"}: with
#'   \eqn{p_1 = a/N} (test-positive and reference-positive) and
#'   \eqn{p_2 = b/N} (test-positive, reference-negative),
#'   \deqn{\widehat{Var}(\widehat{NB}) = \{p_1(1-p_1) + w^2 p_2(1-p_2) + 2 w p_1 p_2\}/N,}
#'   and, since \eqn{\widehat{D} = (w \cdot d - c)/N} depends only on the
#'   false negatives \eqn{c} and true negatives \eqn{d},
#'   \deqn{\widehat{Var}(\widehat{D}) = \{q_1(1-q_1) + w^2 q_2(1-q_2) + 2 w q_1 q_2\}/N,}
#'   with \eqn{q_1 = c/N} and \eqn{q_2 = d/N}. Versions <= 0.2.0 instead
#'   required only that the \emph{point estimate} of NB exceed 0 and the
#'   treat-all NB -- a criterion satisfied with near-certainty whenever the
#'   test is useful at all, and one that \code{\link{ss_net_benefit}} itself
#'   documents as abandoned.
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
#' @param check_nb Logical, require a conclusive (CI-based) net benefit.
#'   Default \code{FALSE}. See \code{Details}.
#' @param pt_range Threshold range for NB if \code{check_nb}.
#'   Default \code{c(0.15, 0.40)}.
#' @param target_assurance Joint assurance target. Default 0.80.
#' @param N_range Search range for total N. Default \code{seq(200, 1000, by = 20)}.
#' @param B MC replications. Default 5000. A warning is issued when
#'   \code{0 < B < 1000}, since the Monte Carlo error of the reported joint
#'   assurance may then be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @param full_grid Logical. Default \code{FALSE}, which stops the search
#'   over \code{N_range} at the first N that reaches \code{target_assurance}
#'   -- identical behaviour and cost to versions <= 0.3.0. If \code{TRUE},
#'   the search does not stop early: every N in \code{N_range} is evaluated
#'   so that \code{grid_results} holds the complete assurance curve.
#'   \code{optimal_N} and \code{joint_assurance} are unaffected by this
#'   switch -- they always refer to the \strong{first} N that reached
#'   \code{target_assurance}, never the last.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{joint_assurance}{Achieved joint assurance at \code{N_effective}:
#'       the probability that \strong{all} active targets -- Se, Sp, AUC
#'       (when \code{delta_auc > 0}) and net benefit (when
#'       \code{check_nb = TRUE}) -- are reached \strong{simultaneously} in
#'       the same replication, not the probability that each is reached
#'       marginally. The denominator is \code{B}: degenerate replications
#'       count as failures.}
#'     \item{comparison}{Data frame comparing methods.}
#'     \item{N_buderer}{Total N from the classical Buderer formula, i.e.
#'       the "Buderer (classical)" row of \code{comparison}.}
#'     \item{N_imperfect}{Total N from the Rogan-Gladen imperfect-reference
#'       inflation, i.e. the "Imperfect ref" row of \code{comparison}.}
#'     \item{seed}{The \code{seed} argument used for the search.}
#'     \item{target_assurance}{The \code{target_assurance} argument used for
#'       the search, echoed back for downstream use (e.g. by
#'       \code{\link{plot_assurance_curve}}).}
#'     \item{grid_results}{Data frame with columns \code{N} and
#'       \code{assurance}, one row per N evaluated while searching
#'       \code{N_range}. With the default \code{full_grid = FALSE} the
#'       search stops at \code{N_effective}, so the grid is truncated there;
#'       set \code{full_grid = TRUE} for the full curve over \code{N_range}.}
#'   }
#' @note \strong{The sensitivity being sized is the \emph{apparent}
#'   sensitivity.} This framework estimates Se and Sp against the imperfect
#'   reference standard \eqn{R}, not against true disease status \eqn{D}:
#'   the quantity whose confidence interval is being narrowed is
#'   \eqn{P(T+\,|\,R+)}, not \eqn{P(T+\,|\,D+)}. Under conditional
#'   independence of \eqn{T} and \eqn{R} given \eqn{D}, the apparent
#'   sensitivity converges to a value \strong{below} the true Se whenever
#'   the reference standard is imperfect. Sizing the study by this framework
#'   buys \emph{precision} about the apparent sensitivity; it does not
#'   remove that \emph{bias}. Recovering an unbiased estimate of the true Se
#'   additionally requires a bias correction or a latent-class analysis at
#'   the analysis stage. See \code{\link{ss_imperfect_ref}}, whose
#'   \code{mc_validation} table quantifies the gap.
#' @references
#' O'Hagan A, Stevens JW, Campbell MJ (2005). Assurance in clinical trial
#' design. \emph{Pharm Stat} 4:187-201. \doi{10.1002/pst.175}
#'
#' Wilson KJ et al. (2022). Bayesian sample size determination for
#' diagnostic accuracy studies. \emph{Stat Med} 41:2908-2922.
#' \doi{10.1002/sim.9393}
#'
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. \emph{Am J Epidemiol} 107:71-76.
#' \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' Hanley JA, McNeil BJ (1982). The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. \emph{Radiology}
#' 143:29-36. \doi{10.1148/radiology.143.1.7063747}
#' @examples
#' \donttest{
#' # Unified, assurance-based sample size (small B and coarse grid for speed;
#' # delta_auc = 0 skips the AUC constraint).
#' result <- suppressWarnings(ss_unified(
#'   prior_se = c(17, 3), prior_sp = c(18, 2),
#'   Se_ref = 0.95, Sp_ref = 0.98,
#'   delta_se = 0.07, delta_sp = 0.06, delta_auc = 0,
#'   N_range = seq(300, 900, by = 100), B = 100
#' ))
#' print(result)
#' result$comparison
#' }
#' @seealso \code{\link{ss_net_benefit}}, \code{\link{ss_imperfect_ref}}
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
                       seed = 2026,
                       full_grid = FALSE) {
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
  stopifnot(length(prior_se) == 2, all(prior_se > 0))
  stopifnot(length(prior_sp) == 2, all(prior_sp > 0))
  stopifnot(length(prior_prev) == 2, all(prior_prev > 0))
  stopifnot(Se_ref > 0, Se_ref <= 1, Sp_ref > 0, Sp_ref <= 1)
  stopifnot(Se_ref + Sp_ref > 1)
  stopifnot(loss_rate >= 0, loss_rate < 1)
  stopifnot(delta_se > 0, delta_sp > 0, delta_auc >= 0)
  stopifnot(all(pt_range > 0 & pt_range < 1))
  stopifnot(target_assurance > 0, target_assurance < 1)
  stopifnot(B >= 1)
  warn_small_B(B)

  z <- stats::qnorm(0.975)
  target_se_width <- 2 * delta_se
  target_sp_width <- 2 * delta_sp
  target_auc_width <- 2 * delta_auc

  # Prior means for comparison table
  E_prev <- prior_prev[1] / sum(prior_prev)
  E_se <- prior_se[1] / sum(prior_se)
  E_sp <- prior_sp[1] / sum(prior_sp)  # Sp arm for the Buderer row

  optimal_N <- NA_integer_
  joint_assurance_achieved <- NA_real_
  joint_assurance <- 0  # defensive init (in case N_range is degenerate)

  # one (N, assurance) pair per N evaluated; with the default
  # full_grid = FALSE the loop below still breaks at the optimum, so this
  # is only as long as the grid actually searched (see @param full_grid).
  grid_N <- vector("integer", length(N_range))
  grid_assurance <- vector("numeric", length(N_range))
  grid_idx <- 0L

  for (N in N_range) {
    set.seed(seed)
    pass_count <- 0L

    for (b in seq_len(B)) {
      # Draw parameters from priors
      prev_b <- stats::rbeta(1, prior_prev[1], prior_prev[2])
      Se_true <- stats::rbeta(1, prior_se[1], prior_se[2])
      Sp_true <- stats::rbeta(1, prior_sp[1], prior_sp[2])

      # Generate true disease status
      n_d <- stats::rbinom(1, N, prev_b)
      n_nd <- N - n_d
      # a degenerate replication is a FAILURE, not an exclusion. The
      # denominator of the assurance stays B (see @details).
      if (n_d < 5 || n_nd < 5) next

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

      # likewise a failure, not an exclusion.
      if (n_ref_pos < 2 || n_ref_neg < 2) next

      # Observed Se_obs = P(T+|R+), Sp_obs = P(T-|R-). NOTE: these are the
      # APPARENT accuracies, measured against R and not against D. See @note.
      Se_obs <- a / n_ref_pos
      Sp_obs <- d_cell / n_ref_neg

      # Wilson CI for Se_obs
      denom_se <- 1 + z^2 / n_ref_pos
      margin_se <- z * sqrt((Se_obs * (1 - Se_obs) +
        z^2 / (4 * n_ref_pos)) / n_ref_pos) / denom_se
      width_se <- 2 * margin_se
      se_pass <- (width_se <= target_se_width)

      # Wilson CI for Sp_obs
      denom_sp <- 1 + z^2 / n_ref_neg
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

      # --- Net Benefit check, now INFERENCE-based ---
      # N is fixed by design but reference-positive status is random, so the
      # multinomial (cohort) variances apply, exactly as in
      # ss_net_benefit(design = "cohort"). The treat-all net benefit is also
      # estimated, so the treat-all comparison is made on the difference
      # D = NB - NB_all = (w * d_cell - c_cell) / N.
      nb_pass <- TRUE
      if (check_nb) {
        for (pt in pt_range) {
          w <- pt / (1 - pt)

          # NB vs treat-none (0)
          p1 <- a / N          # T+ & R+  (apparent true positives)
          p2 <- b_cell / N     # T+ & R-  (apparent false positives)
          nb_hat <- p1 - w * p2
          var_nb <- (p1 * (1 - p1) + w^2 * p2 * (1 - p2) +
                       2 * w * p1 * p2) / N
          nb_lower <- nb_hat - z * sqrt(max(var_nb, 0))

          # NB vs treat-all, via the difference (depends only on c and d)
          q1 <- c_cell / N     # T- & R+  (apparent false negatives)
          q2 <- d_cell / N     # T- & R-  (apparent true negatives)
          d_hat <- w * q2 - q1
          var_d <- (q1 * (1 - q1) + w^2 * q2 * (1 - q2) +
                      2 * w * q1 * q2) / N
          d_lower <- d_hat - z * sqrt(max(var_d, 0))

          if (!(nb_lower > 0 && d_lower > 0)) {
            nb_pass <- FALSE
            break
          }
        }
      }

      all_pass <- se_pass & sp_pass & auc_pass & nb_pass
      if (all_pass) pass_count <- pass_count + 1L
    }

    # denominator is B, not the number of non-degenerate replications.
    joint_assurance <- pass_count / B

    grid_idx <- grid_idx + 1L
    grid_N[grid_idx] <- N
    grid_assurance[grid_idx] <- joint_assurance

    # optimal_N is always the FIRST N to reach the target: once set, later
    # N (only reachable with full_grid = TRUE) must not overwrite it.
    if (is.na(optimal_N) && joint_assurance >= target_assurance) {
      optimal_N <- as.integer(N)
      joint_assurance_achieved <- joint_assurance
      if (!isTRUE(full_grid)) break
    }
  }

  grid_results <- data.frame(
    N = grid_N[seq_len(grid_idx)],
    assurance = grid_assurance[seq_len(grid_idx)],
    stringsAsFactors = FALSE
  )

  if (is.na(optimal_N)) {
    warning("No N in N_range achieved target assurance. ",
            "Consider expanding N_range.")
    optimal_N <- max(N_range)
    joint_assurance_achieved <- joint_assurance
  }

  N_enrolled <- ceiling(optimal_N / (1 - loss_rate))

  # Comparison table
  buderer_N_total <- buderer_total_N(E_se, E_sp, delta_se, delta_sp, E_prev)
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
    method = c("Buderer (classical)",
               "Imperfect ref (Rogan-Gladen inflation)",
               "Unified (this method)"),
    N = c(buderer_N_total, imperfect_res$N_adjusted_loss, N_enrolled),
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
      N_buderer = buderer_N_total,
      N_imperfect = imperfect_res$N_adjusted_loss,
      seed = seed,
      target_assurance = target_assurance,
      grid_results = grid_results,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
