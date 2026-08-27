
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
#'   \strong{Internal pilot design.} Following Stark and Zapf (2020), stage 1
#'   is a genuine \emph{internal} pilot, not a throw-away sample used only to
#'   estimate the prevalence and then discarded: the \code{n_stage1}
#'   subjects recruited in stage 1 are retained and folded into the final
#'   analysed sample. Concretely, stage 2 recruits only the
#'   \strong{additional} subjects needed to reach the re-estimated total,
#'   \code{N_final - n_stage1} (never \code{N_final} again from scratch),
#'   and the final diseased count is \code{D_stage1} (already observed in
#'   stage 1) \strong{plus} the diseased subjects drawn in stage 2. This
#'   preserves the correlation between the interim prevalence estimate and
#'   the analysed sample that a real internal pilot has by construction (the
#'   same patients feed both). In versions <= 0.5.0 stage 1 was simulated
#'   and then discarded entirely: the reported \code{N_final_*} undercounted
#'   the true recruitment by \code{n_stage1} subjects (about 20-36% at the
#'   package defaults), and \code{prev_hat} was statistically independent of
#'   the analysed sample -- which is not what an internal-pilot design does.
#'
#'   \strong{Losses to follow-up.} Losses are modelled on the stage-2
#'   increment only. Each replicate recruits
#'   \code{n_stage2_adj = ceiling((N_final - n_stage1) / (1 - loss_rate))}
#'   additional subjects and \strong{loses} a random
#'   \eqn{Bin(n\_stage2\_adj, loss\_rate)} of them; the stage-1 subjects are
#'   not re-subjected to loss because \code{D_stage1} is, by construction,
#'   already fully observed (that is what makes it usable to re-estimate the
#'   prevalence at the interim look) -- mirroring how stage 1 was already
#'   exempt from the loss adjustment before this fix. The accuracy estimates
#'   and their confidence intervals are computed on the surviving (analysed)
#'   sample, \code{n_stage1 + (n_stage2_adj - n_lost)}, only. The reported
#'   \code{N_final_*} columns are the \strong{recruited} sample sizes across
#'   both stages (what an investigator must enrol in total) and
#'   \code{N_analysed_median} is the median analysed sample size.
#'   Consequently \code{precision_achieved} is approximately \strong{invariant}
#'   to \code{loss_rate}: inflating the stage-2 recruitment for losses
#'   restores the precision that the losses destroy, it does not improve on
#'   it. In versions <= 0.2.0 the inflated sample was analysed in full (the
#'   loss never happened), so the reported precision \emph{rose} with the
#'   loss rate -- an artefact.
#'
#'   \strong{Interim prevalence truncation.} The stage-1 prevalence estimate
#'   is truncated to \eqn{[0.05, 0.95]} before it is used to re-estimate the
#'   sample size. The lower bound prevents a small or zero stage-1 numerator
#'   from demanding an unbounded number of diseased subjects
#'   (\code{n_se / prev_hat}); the upper bound prevents the symmetric
#'   failure on the non-diseased arm (\code{n_sp / (1 - prev_hat)}), which
#'   previously returned \code{Inf} whenever every stage-1 subject was
#'   diseased. The truncation is a deliberate, documented cap on the
#'   adaptive rule, not a statistical estimate: with a stage-1 sample this
#'   small, prevalence estimates at the boundary are not credible, and the
#'   cap bounds the re-estimated N at \code{max(n_se / 0.05, n_sp / 0.05)}.
#'
#'   \strong{Stage-1 fraction.} For the one-time re-estimation design, Stark
#'   and Zapf (2020) recommend an internal pilot of \strong{50\%} of the
#'   initially calculated sample size: \dQuote{The appropriate size of the
#'   internal pilot study in the one-time re-estimation design is 50\% of
#'   the initially calculated sample size.} The package default,
#'   \code{fraction_stage1 = 0.40}, is smaller than that recommendation and
#'   is left unchanged here; pass \code{fraction_stage1 = 0.5} to follow the
#'   published recommendation.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param d_se Precision for Se (half-width). Default 0.07.
#' @param d_sp Precision for Sp (half-width). Default 0.05.
#' @param prev_initial Initially assumed prevalence. Default 0.30.
#' @param prev_true_range Numeric vector of true prevalence scenarios.
#'   Default \code{c(0.18, 0.25, 0.30, 0.35, 0.42)}.
#' @param fraction_stage1 Fraction of initial N recruited as the stage-1
#'   internal pilot (see \code{Details}). Default 0.40. Stark and Zapf
#'   (2020) recommend 0.50 for the one-time re-estimation design; the
#'   default here is left at 0.40 for backward compatibility -- pass
#'   \code{fraction_stage1 = 0.5} to follow the published recommendation.
#' @param loss_rate Expected loss-to-follow-up rate. Default 0.10. Losses
#'   are simulated: the analysed sample is the recruited sample minus a
#'   random binomial number of losses. See \code{Details}.
#' @param B MC replications. Default 5000. A warning is issued when
#'   \code{0 < B < 1000}, since the Monte Carlo error of
#'   \code{precision_achieved} may then be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{N_initial}{Initial planned N (before adaptation).}
#'     \item{n_stage1}{Number of subjects in stage 1.}
#'     \item{results}{Data frame with one row per true prevalence scenario:
#'       \code{prev_true}, the mean/median/P75 \strong{recruited} sample size
#'       (\code{N_final_mean}, \code{N_final_median}, \code{N_final_P75}),
#'       the median \strong{analysed} sample size after losses
#'       (\code{N_analysed_median}), and \code{precision_achieved}, the
#'       probability that both the Se and the Sp confidence intervals meet
#'       their target width \emph{in the analysed sample}.}
#'   }
#' @note \code{precision_achieved} is computed on the post-loss (analysed)
#'   sample, so it does not increase with \code{loss_rate}. Any residual
#'   variation across loss rates is Monte Carlo noise plus the effect of
#'   rounding the inflated sample size up to an integer.
#' @references
#' Stark M, Zapf A (2020). Sample size calculation and re-estimation based
#' on the prevalence in a single-arm confirmatory diagnostic accuracy
#' study. \emph{Stat Methods Med Res} 29:2958-2971.
#' \doi{10.1177/0962280220913588}
#' @examples
#' # B kept small here for a fast example; see @param B.
#' result <- suppressWarnings(
#'   ss_adaptive_prevalence(B = 500, prev_true_range = c(0.20, 0.30))
#' )
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
  # --- preserve the caller's RNG state (kind AND seed) ------------------
  # See save_rng_state()/restore_rng_state(): restoring only .Random.seed's
  # VALUE is not enough, because set.seed() called later by unrelated code
  # with no explicit `kind` argument reuses whichever kind is CURRENTLY
  # ACTIVE. The set.seed() call below names its kind explicitly
  # (Mersenne-Twister, R's own default), so the adaptive-design simulation
  # reproduces the same numbers regardless of the caller's own RNG
  # configuration.
  old_rng_state <- save_rng_state()
  on.exit(restore_rng_state(old_rng_state), add = TRUE)

  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(d_se > 0, d_sp > 0)
  stopifnot(prev_initial > 0, prev_initial < 1)
  stopifnot(all(prev_true_range > 0 & prev_true_range < 1))
  stopifnot(fraction_stage1 > 0, fraction_stage1 < 1)
  stopifnot(loss_rate >= 0, loss_rate < 1)
  stopifnot(B >= 1)
  warn_small_B(B)

  # Required sample sizes for Se and Sp (Buderer)
  n_se <- buderer_n(Se, d_se)
  n_sp <- buderer_n(Sp, d_sp)

  # Initial N
  N_initial <- buderer_total_N(Se, Sp, d_se, d_sp, prev_initial)
  N_initial_adj <- ceiling(N_initial / (1 - loss_rate))

  # Stage 1 sample size
  n_stage1 <- ceiling(fraction_stage1 * N_initial)

  # bounds on the interim prevalence estimate (see @details)
  PREV_HAT_MIN <- 0.05
  PREV_HAT_MAX <- 0.95

  z <- stats::qnorm(1 - 0.05 / 2)
  target_se_width <- 2 * d_se
  target_sp_width <- 2 * d_sp

  results_list <- vector("list", length(prev_true_range))

  for (i in seq_along(prev_true_range)) {
    prev_true <- prev_true_range[i]
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

    N_finals <- numeric(B)
    N_analysed <- numeric(B)
    precision_met <- logical(B)
    n_truncated <- 0L

    for (b in seq_len(B)) {
      # --- Stage 1: estimate prevalence ---
      D_stage1 <- stats::rbinom(1, n_stage1, prev_true)
      # truncate to [0.05, 0.95]. Without the upper bound, prev_hat = 1
      # makes n_sp / (1 - prev_hat) infinite. When the truncation binds it is
      # counted so the caller can be warned (an interim estimate outside the
      # bounds can under-size the diseased arm).
      prev_raw <- D_stage1 / n_stage1
      prev_hat <- min(max(prev_raw, PREV_HAT_MIN), PREV_HAT_MAX)
      if (prev_raw < PREV_HAT_MIN || prev_raw > PREV_HAT_MAX) {
        n_truncated <- n_truncated + 1L
      }

      # --- Re-estimate N ---
      N_final <- max(buderer_total_N(Se, Sp, d_se, d_sp, prev_hat), N_initial)

      # --- Stage 2: recruit only the ADDITIONAL subjects needed ---------
      # n_stage1 is an internal pilot (see @details "Internal pilot
      # design"): its subjects are already part of the analysed sample, so
      # stage 2 tops up to N_final rather than recruiting N_final again.
      # n_stage2_needed is >= 0 by construction (N_final >= N_initial >
      # n_stage1 whenever 0 < fraction_stage1 < 1), but max(..., 0) guards
      # the edge case fraction_stage1 -> 1 where N_final == n_stage1.
      n_stage2_needed <- max(N_final - n_stage1, 0)
      n_stage2_adj <- if (n_stage2_needed > 0) {
        ceiling(n_stage2_needed / (1 - loss_rate))
      } else {
        0L
      }

      # Recruited sample size (what the investigator must enrol, in total,
      # across both stages)
      N_final_adj <- n_stage1 + n_stage2_adj
      N_finals[b] <- N_final_adj

      # --- the losses actually HAPPEN, on the stage-2 increment only ----
      # (see @details "Losses to follow-up" for why stage 1 is exempt)
      n_lost <- if (loss_rate > 0 && n_stage2_adj > 0) {
        stats::rbinom(1, n_stage2_adj, loss_rate)
      } else {
        0L
      }
      n_stage2_analysis <- n_stage2_adj - n_lost
      N_analysis <- n_stage1 + n_stage2_analysis
      N_analysed[b] <- N_analysis

      if (N_analysis < 2) {
        precision_met[b] <- FALSE
        next
      }

      # --- Analyse the POST-LOSS sample: the D_stage1 diseased subjects
      # already observed in stage 1, PLUS a fresh draw among the surviving
      # stage-2 subjects ---------------------------------------------------
      n_d_stage2 <- if (n_stage2_analysis > 0) {
        stats::rbinom(1, n_stage2_analysis, prev_true)
      } else {
        0L
      }
      n_d <- D_stage1 + n_d_stage2
      n_nd <- N_analysis - n_d

      # A degenerate arm cannot meet the precision target: count as failure
      # rather than silently pretending one subject was observed.
      if (n_d < 1 || n_nd < 1) {
        precision_met[b] <- FALSE
        next
      }

      se_obs <- stats::rbinom(1, n_d, Se) / n_d
      sp_obs <- stats::rbinom(1, n_nd, Sp) / n_nd

      # Wilson CI for Se
      denom_se <- 1 + z^2 / n_d
      margin_se <- z * sqrt((se_obs * (1 - se_obs) +
        z^2 / (4 * n_d)) / n_d) / denom_se
      width_se <- 2 * margin_se

      # Wilson CI for Sp
      denom_sp <- 1 + z^2 / n_nd
      margin_sp <- z * sqrt((sp_obs * (1 - sp_obs) +
        z^2 / (4 * n_nd)) / n_nd) / denom_sp
      width_sp <- 2 * margin_sp

      precision_met[b] <- (width_se <= target_se_width) &
        (width_sp <= target_sp_width)
    }

    if (n_truncated > 0L) {
      warning(sprintf(
        paste0("ss_adaptive_prevalence: at true prevalence %.3f, the interim ",
               "estimate was truncated to [%.2f, %.2f] in %d of %d replicates ",
               "(%.0f%%); the diseased-arm sample size may be under-sized when ",
               "the true prevalence lies outside that range."),
        prev_true, PREV_HAT_MIN, PREV_HAT_MAX, n_truncated, B,
        100 * n_truncated / B), call. = FALSE)
    }

    results_list[[i]] <- data.frame(
      prev_true = prev_true,
      N_final_mean = mean(N_finals),
      N_final_median = stats::median(N_finals),
      N_final_P75 = stats::quantile(N_finals, 0.75, names = FALSE),
      N_analysed_median = stats::median(N_analysed),
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
      loss_rate = loss_rate,
      results = results,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
