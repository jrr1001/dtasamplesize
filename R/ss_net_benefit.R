
#' Sample Size for Net Benefit (Decision Curve Analysis)
#'
#' Finds the minimum total sample size such that a study can
#' \emph{conclude}, with \code{(1 - alpha)} confidence, that the index
#' test has positive net benefit and outperforms the default strategies
#' (treat-all and treat-none) across a range of threshold probabilities.
#'
#' @details For each threshold probability \code{pt} (write
#'   \eqn{w = pt / (1 - pt)} for the threshold odds), the population net
#'   benefit of the test and of the treat-all strategy are
#'   \deqn{NB = prev \cdot Se - (1 - prev)(1 - Sp)\,w,}
#'   \deqn{NB_{all} = prev - (1 - prev)\,w.}
#'   The treat-none strategy has \eqn{NB = 0} by definition. A test is
#'   useful at \code{pt} only when \eqn{NB > 0} \strong{and}
#'   \eqn{NB > NB_{all}}. When the test is \strong{not} useful at a
#'   threshold under the assumed parameters, no sample size can demonstrate
#'   superiority; the threshold is flagged \code{feasible = FALSE} and
#'   \code{N_required = NA}.
#'
#'   The sample-size criterion is \strong{inference-based}: a simulated
#'   study is a "success" when the lower confidence limit of the relevant
#'   estimator excludes the competing strategy. \code{N_required} is the
#'   smallest \code{N} at which this happens in at least
#'   \code{target_prob} of the Monte Carlo studies ("assurance"). The
#'   lower limit used is that of the two-sided \code{(1 - alpha)} normal
#'   confidence interval, i.e. a one-sided test at level \code{alpha / 2}.
#'
#'   \strong{The two sampling designs differ in what is random, and
#'   therefore in the variance of the estimator.}
#'
#'   \code{design = "cohort"} (default) --- \emph{prospective cohort /
#'   consecutive series}. Disease status is \strong{random}: a single
#'   subject falls into one of four cells with probabilities
#'   \eqn{p_{TP} = prev \cdot Se}, \eqn{p_{FN} = prev (1 - Se)},
#'   \eqn{p_{FP} = (1 - prev)(1 - Sp)} and \eqn{p_{TN} = (1 - prev) Sp},
#'   so \eqn{(TP, FN, FP, TN)} is multinomial with total \code{N}. Each
#'   replicate draws \eqn{n_d \sim Bin(N, prev)},
#'   \eqn{TP \sim Bin(n_d, Se)}, \eqn{FP \sim Bin(N - n_d, 1 - Sp)}, and
#'   \eqn{\widehat{NB} = TP/N - (FP/N)\,w}. Writing \eqn{p_1 = TP/N} and
#'   \eqn{p_2 = FP/N}, the (multinomial, plug-in) variance is
#'   \deqn{\widehat{Var}(\widehat{NB}) = \{p_1(1 - p_1) + w^2 p_2 (1 - p_2) + 2 w p_1 p_2\}/N,}
#'   the last term arising because \eqn{TP} and \eqn{FP} are negatively
#'   correlated when \code{N} (not \eqn{n_d}) is fixed.
#'
#'   Under this design the prevalence -- and hence \eqn{NB_{all}} -- is
#'   \strong{also estimated}, so \eqn{NB_{all}} is not a known constant and
#'   the treat-all comparison must be made on the \strong{difference}
#'   \eqn{D = NB - NB_{all}}. Substituting the estimators gives the exact
#'   algebraic identity
#'   \deqn{\widehat{D} = \widehat{NB} - \widehat{NB}_{all} = (w \cdot TN - FN)/N,}
#'   i.e. the difference depends \strong{only} on the false negatives and
#'   the true negatives (the treat-all and test strategies agree on every
#'   test-positive subject, so those subjects cancel). With
#'   \eqn{q_1 = FN/N} and \eqn{q_2 = TN/N} its plug-in variance is, by the
#'   same multinomial argument,
#'   \deqn{\widehat{Var}(\widehat{D}) = \{q_1(1 - q_1) + w^2 q_2 (1 - q_2) + 2 w q_1 q_2\}/N.}
#'   A cohort study is counted a success when the lower confidence limit of
#'   \eqn{\widehat{NB}} exceeds 0 \strong{and} the lower confidence limit
#'   of \eqn{\widehat{D}} exceeds 0. Both closed-form variances were
#'   verified against the Monte Carlo sampling distribution of the
#'   corresponding estimator under the cohort data-generating process.
#'
#'   \code{design = "fixed"} --- \emph{fixed disease-status margins} (the
#'   behaviour of package versions <= 0.2.0). \eqn{n_d = \lfloor N \cdot prev \rfloor}
#'   and \eqn{n_{nd} = N - n_d} are fixed by design, so \eqn{NB_{all}} is a
#'   known constant and the variance is the conditional one,
#'   \deqn{Var(\widehat{NB}) = (n_d/N)^2 Se(1 - Se)/n_d + (n_{nd}/N)^2 w^2 Sp(1 - Sp)/n_{nd}.}
#'   A study is a success when the lower confidence limit of
#'   \eqn{\widehat{NB}} exceeds both 0 and \eqn{NB_{all}}. This variance
#'   ignores the sampling variability of the disease prevalence and is
#'   therefore \strong{substantially too small for a prospective cohort};
#'   see the \code{Note}.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param prev Disease prevalence. Default 0.20.
#' @param design Sampling design, \code{"cohort"} (default) or
#'   \code{"fixed"}. \code{"cohort"} treats disease status as random, as in
#'   a prospective cohort or consecutive series, and is the appropriate
#'   choice for almost all diagnostic accuracy studies. \code{"fixed"}
#'   conditions on \eqn{n_d = \lfloor N \cdot prev \rfloor} diseased
#'   subjects and reproduces the behaviour of versions <= 0.2.0; it applies
#'   only to a design that recruits the two disease groups separately with
#'   pre-specified sizes, and it is \strong{not} valid for a prospective
#'   cohort. See \code{Details}.
#' @param pt_range Threshold probabilities to evaluate.
#'   Default \code{seq(0.10, 0.50, 0.05)}.
#' @param target_prob Assurance: probability that the lower confidence
#'   limit excludes the competing default strategies. Default 0.80.
#' @param target_assurance Alias for \code{target_prob}, matching the
#'   \code{target_assurance} naming used elsewhere in the package (e.g.
#'   \code{\link{ss_unified}}, \code{\link{bam_sample_size}}). If supplied
#'   (non-\code{NULL}), it takes precedence over \code{target_prob}.
#'   Default \code{NULL} (use \code{target_prob}).
#' @param alpha Confidence-interval significance level. Default 0.05. The
#'   lower limit of the two-sided \code{(1 - alpha)} interval is used, i.e.
#'   a one-sided test at level \code{alpha / 2}.
#' @param N_range Range of total N to search.
#'   Default \code{seq(50, 2500, by = 10)}.
#' @param B MC replications. Default 5000. A warning is issued when
#'   \code{0 < B < 1000}, since the Monte Carlo error of the reported
#'   assurance may then be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{design}{The sampling design used.}
#'     \item{N_by_pt}{Data frame with columns \code{pt}, \code{feasible},
#'       \code{N_required}, \code{prob_achieved}, \code{NB_true},
#'       \code{NB_treat_all}.}
#'     \item{N_conservative}{Maximum required N across feasible thresholds
#'       (\code{NA} if none feasible, or if a feasible threshold did not
#'       converge within \code{N_range}).}
#'   }
#' @note \strong{Design matters, a lot.} Under the default cohort design the
#'   variance of the net-benefit estimator is roughly twice that of the
#'   \code{"fixed"} design at typical parameters, because the latter
#'   conditions away the sampling variability of the disease prevalence. At
#'   \code{Se = 0.85}, \code{Sp = 0.90}, \code{prev = 0.20}, \code{pt = 0.20}
#'   and \code{N = 500} the fixed-design standard error is 0.0077 whereas
#'   the true cohort sampling standard deviation is 0.0174. Sizing a
#'   prospective cohort with \code{design = "fixed"} therefore
#'   \strong{overstates the assurance} and yields a sample size that is far
#'   too small. Versions <= 0.2.0 offered only the fixed-margin variance.
#'
#'   The sampling distributions are approximated by the normal
#'   distribution; at very small \eqn{n_d} or \eqn{n_{nd}} this
#'   approximation degrades. Replicates that are degenerate (no diseased or
#'   no non-diseased subjects) are counted as \strong{failures}, not
#'   discarded, so the reported assurance is unconditional.
#' @references
#' Vickers AJ, Elkin EB (2006). Decision curve analysis: a novel method
#' for evaluating prediction models. \emph{Med Decis Making} 26:565-574.
#' \doi{10.1177/0272989X06295361}
#'
#' Vickers AJ, van Calster B, Steyerberg EW (2019). A simple,
#' step-by-step guide to interpreting decision curve analysis.
#' \emph{Diagn Progn Res} 3:18. \doi{10.1186/s41512-019-0064-7}
#' @examples
#' \donttest{
#' # N_range and B are bounded here to keep the example fast; the default
#' # search range is seq(50, 2500, by = 10) and default B is 5000 (see
#' # @param B).
#' result <- suppressWarnings(ss_net_benefit(pt_range = c(0.10, 0.30),
#'                          N_range = seq(50, 400, by = 10), B = 500))
#' print(result)
#'
#' # Legacy fixed-margin design (not valid for a prospective cohort)
#' old <- suppressWarnings(ss_net_benefit(pt_range = 0.30, design = "fixed",
#'                       N_range = seq(50, 400, by = 10), B = 500))
#' }
#' @export
ss_net_benefit <- function(Se = 0.85,
                           Sp = 0.90,
                           prev = 0.20,
                           design = c("cohort", "fixed"),
                           pt_range = seq(0.10, 0.50, 0.05),
                           target_prob = 0.80,
                           alpha = 0.05,
                           N_range = seq(50, 2500, by = 10),
                           B = 5000,
                           seed = 2026,
                           target_assurance = NULL) {
  # target_assurance is an alias of target_prob, matching the naming used
  # elsewhere in the package; default behaviour (NULL) is unchanged.
  if (!is.null(target_assurance)) target_prob <- target_assurance

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
  design <- match.arg(design)
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(prev > 0, prev < 1)
  stopifnot(all(pt_range > 0 & pt_range < 1))
  stopifnot(target_prob > 0, target_prob < 1)
  stopifnot(alpha > 0, alpha < 1)
  stopifnot(B >= 1)
  warn_small_B(B)

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
        set.seed(seed)

        if (design == "cohort") {
          # Disease status is RANDOM: (TP, FN, FP, TN) ~ Multinomial(N, .)
          n_d <- stats::rbinom(B, N, prev)
          n_nd <- N - n_d
          TP <- stats::rbinom(B, n_d, Se)
          FP <- stats::rbinom(B, n_nd, 1 - Sp)
          FN <- n_d - TP
          TN <- n_nd - FP

          # --- NB vs treat-none (0) ---
          p1 <- TP / N
          p2 <- FP / N
          NB_hat <- p1 - p2 * odds_pt
          var_nb <- (p1 * (1 - p1) + odds_pt^2 * p2 * (1 - p2) +
                       2 * odds_pt * p1 * p2) / N
          nb_lower <- NB_hat - z * sqrt(pmax(var_nb, 0))

          # --- NB vs treat-all: NB_all is ALSO estimated, so the criterion
          # is applied to the difference D = NB - NB_all = (w*TN - FN)/N,
          # which depends only on FN and TN.
          q1 <- FN / N
          q2 <- TN / N
          D_hat <- odds_pt * q2 - q1
          var_d <- (q1 * (1 - q1) + odds_pt^2 * q2 * (1 - q2) +
                      2 * odds_pt * q1 * q2) / N
          d_lower <- D_hat - z * sqrt(pmax(var_d, 0))

          # Degenerate replicates (no diseased or no non-diseased) are NOT
          # discarded: their lower limits simply fail, so they count as
          # failures and the assurance stays unconditional.
          success <- (nb_lower > 0) & (d_lower > 0)

        } else {
          # design == "fixed": disease-status margins fixed by design.
          n_d <- floor(N * prev)
          n_nd <- N - n_d
          if (n_d < 2 || n_nd < 2) next

          TP <- stats::rbinom(B, n_d, Se)
          FP <- stats::rbinom(B, n_nd, 1 - Sp)
          Se_hat <- TP / n_d
          Sp_hat <- 1 - FP / n_nd

          w_d <- n_d / N
          w_nd <- n_nd / N
          NB_hat <- w_d * Se_hat - w_nd * (1 - Sp_hat) * odds_pt

          # Variance CONDITIONAL on (n_d, n_nd); TP, FP independent.
          var_nb <- w_d^2 * Se_hat * (1 - Se_hat) / n_d +
            w_nd^2 * odds_pt^2 * Sp_hat * (1 - Sp_hat) / n_nd
          nb_lower <- NB_hat - z * sqrt(pmax(var_nb, 0))

          # NB_treat_all is a known constant under this design.
          success <- (nb_lower > 0) & (nb_lower > NB_treat_all)
        }

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

  # n_diseased / n_total derived from the conservative N (if any). Under the
  # cohort design this is the EXPECTED number of diseased subjects.
  n_d_final <- if (is.na(N_conservative)) NA_integer_ else floor(N_conservative * prev)

  structure(
    list(
      method = "Sample Size for Net Benefit (Decision Curve Analysis)",
      design = design,
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
