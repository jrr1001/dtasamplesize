#' Sample Size for Time-Dependent ROC Analysis
#'
#' Monte Carlo-based sample size estimation for achieving a target
#' precision of the area under the time-dependent ROC curve, AUC(t),
#' at a given time horizon. Uses the \pkg{timeROC} package for
#' AUC(t) estimation with inverse-probability-of-censoring weighting.
#'
#' @param mu_case Mean biomarker value in cases. Default 4.5.
#' @param mu_control Mean biomarker value in controls. Default 2.5.
#' @param sigma SD of biomarker (both groups). Default 2.5.
#' @param lambda_event Event rate per year (exponential). Default 0.25.
#' @param t_horizon Time point for AUC(t) evaluation in years. Default 2.
#' @param censoring_rates Numeric vector giving, for each scenario, the
#'   target marginal probability that the \emph{censoring time} precedes
#'   \code{t_horizon}, i.e. \eqn{P(C < t\_horizon)}. Each entry calibrates
#'   an exponential censoring-time model with rate
#'   \code{-log(1 - censoring_rates) / t_horizon}; it is \strong{not} the
#'   fraction of subjects who will actually be observed as censored. See
#'   Details. Default \code{c(0.10, 0.20, 0.30)}.
#' @param delta_auc Target precision (half-width) for AUC(t). Default 0.06.
#' @param target_prob Probability of achieving precision. Default 0.80.
#' @param N_range Range of N to search. Default \code{seq(100, 500, by = 20)}.
#' @param B MC replications. Default 500 (lower due to timeROC cost). A
#'   warning is issued when \code{0 < B < 1000}, since the Monte Carlo error
#'   of the reported probabilities may then be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @details \code{censoring_rates} parameterizes the exponential model for
#'   the censoring time \eqn{C} alone: each entry \code{cens_rate} sets
#'   \eqn{P(C < t\_horizon) = cens\_rate} exactly, via
#'   \code{lambda_censor <- -log(1 - cens_rate) / t_horizon}. It does
#'   \strong{not} equal the proportion of subjects observed as censored in
#'   the simulated data, because censoring competes with the event: a
#'   subject is recorded as censored only when \eqn{C} occurs before both
#'   the event time \eqn{T} and \code{t_horizon}. Since some subjects who
#'   would have been censored experience the event first, the observed
#'   censoring proportion is systematically \strong{lower} than
#'   \code{censoring_rates}.
#'
#'   For the package defaults (\code{lambda_event = 0.25},
#'   \code{t_horizon = 2}, under which the event alone would be observed
#'   before the horizon for 39.35\% of subjects), Monte Carlo simulation
#'   (2,000,000 subjects per scenario) gives:
#'   \tabular{lll}{
#'     \code{censoring_rates} \tab \eqn{P(C < t\_horizon)} \tab observed censoring proportion \cr
#'     0.10 \tab 0.0996 \tab 0.0787 \cr
#'     0.20 \tab 0.2002 \tab 0.1591 \cr
#'     0.30 \tab 0.3000 \tab 0.2396 \cr
#'     0.50 \tab 0.5003 \tab 0.4054 \cr
#'   }
#'   The observed proportion runs at roughly 80\% of the nominal
#'   \code{censoring_rates} value across this range; the ratio is not
#'   universal and shifts with \code{lambda_event} and \code{t_horizon}
#'   (a higher event rate leaves less "room" for censoring to be
#'   observed, which lowers the ratio further). Treat 80\% as a rough
#'   guide for the package defaults, not a general conversion factor.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{results}{Data frame with columns \code{censoring_rate} (the
#'       nominal \eqn{P(C < t\_horizon)} value from \code{censoring_rates},
#'       not the observed censoring proportion; see Details),
#'       \code{N_required}, \code{prob_achieved}.}
#'   }
#' @note Parameters in the default example are HYPOTHETICAL. No published
#'   AUC values exist for ctDNA as a continuous discriminator in DLBCL.
#'   Requires the \pkg{timeROC} package (listed in Suggests).
#' @references
#' Heagerty PJ, Lumley T, Pepe MS (2000). Time-dependent ROC curves for
#' censored survival data and a diagnostic marker. \emph{Biometrics}
#' 56:337-344. \doi{10.1111/j.0006-341X.2000.00337.x}
#'
#' Blanche P, Dartigues J-F, Jacqmin-Gadda H (2013). Estimating and
#' comparing time-dependent areas under receiver operating characteristic
#' curves for censored event times with competing risks. \emph{Stat Med}
#' 32:5381-5397. \doi{10.1002/sim.5958}
#' @examples
#' \donttest{
#' # B kept small here for a fast example; see @param B.
#' result <- suppressWarnings(
#'   ss_time_dependent_roc(B = 50, N_range = seq(100, 300, by = 50))
#' )
#' print(result)
#' }
#' @export
ss_time_dependent_roc <- function(mu_case = 4.5,
                                  mu_control = 2.5,
                                  sigma = 2.5,
                                  lambda_event = 0.25,
                                  t_horizon = 2,
                                  censoring_rates = c(0.10, 0.20, 0.30),
                                  delta_auc = 0.06,
                                  target_prob = 0.80,
                                  N_range = seq(100, 500, by = 20),
                                  B = 500,
                                  seed = 2026) {
  # --- preserve the caller's RNG state (kind AND seed) ------------------
  # See save_rng_state()/restore_rng_state(): restoring only .Random.seed's
  # VALUE is not enough, because set.seed() called later by unrelated code
  # with no explicit `kind` argument reuses whichever kind is CURRENTLY
  # ACTIVE. The set.seed() call below names ALL THREE kinds explicitly
  # (Mersenne-Twister / Inversion / Rejection, R's own defaults) -- not
  # just the uniform generator. This function is the one place in the
  # package that draws from the NORMAL generator (stats::rnorm(), for the
  # simulated biomarker below), so naming `kind` alone was not enough: a
  # caller on a different `normal.kind` (e.g. "Box-Muller" instead of the
  # default "Inversion") got different simulated biomarker values, and
  # therefore a different N_required, from the identical seed and
  # arguments -- e.g. N_total = 420 under Inversion vs 400 under
  # Box-Muller for otherwise identical calls. Naming all three kinds makes
  # the simulated AUC(t) precision reproduce the same numbers regardless
  # of the caller's own RNG configuration.
  old_rng_state <- save_rng_state()
  on.exit(restore_rng_state(old_rng_state), add = TRUE)

  if (!requireNamespace("timeROC", quietly = TRUE)) {
    stop("Package 'timeROC' is required. Install with: ",
         "install.packages('timeROC')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required. Install with: ",
         "install.packages('survival')")
  }

  # Validate inputs
  stopifnot(sigma > 0, lambda_event > 0, t_horizon > 0)
  stopifnot(all(censoring_rates >= 0 & censoring_rates < 1))
  stopifnot(delta_auc > 0, target_prob > 0, target_prob < 1)
  stopifnot(B >= 1)
  warn_small_B(B)

  # Attach survival so timeROC can find Surv() in formula evaluation
  if (!("package:survival" %in% search())) {
    attachNamespace("survival")
    on.exit(detach("package:survival"), add = TRUE)
  }

  target_width <- 2 * delta_auc
  results_list <- vector("list", length(censoring_rates))

  for (cr_idx in seq_along(censoring_rates)) {
    cens_rate <- censoring_rates[cr_idx]

    # Compute censoring rate parameter
    if (cens_rate > 0) {
      lambda_censor <- -log(1 - cens_rate) / t_horizon
    } else {
      lambda_censor <- 0
    }

    found_N <- NA_integer_
    found_prob <- NA_real_
    best_prob <- 0
    prob <- NA_real_  # assurance at the largest N tried (set in the loop)

    for (N in N_range) {
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
      success_count <- 0L

      for (b in seq_len(B)) {
        ok <- tryCatch({
          # Generate event times
          T_event <- stats::rexp(N, lambda_event)

          # Generate censoring times
          if (lambda_censor > 0) {
            C_time <- stats::rexp(N, lambda_censor)
          } else {
            C_time <- rep(Inf, N)
          }

          # Observed time and event indicator
          Y <- pmin(T_event, C_time)
          delta <- as.integer(T_event <= C_time)

          # Biomarker conditional on case status at t_horizon
          is_case <- (T_event <= t_horizon)
          marker <- ifelse(is_case,
                           stats::rnorm(N, mu_case, sigma),
                           stats::rnorm(N, mu_control, sigma))

          # Need at least some events and non-events before t_horizon
          n_events <- sum(delta == 1 & Y <= t_horizon)
          n_nonevents <- sum(Y > t_horizon)
          if (n_events < 5 || n_nonevents < 5) {
            FALSE
          } else {
            # timeROC triggers the R-devel (R >= 4.6) warning "object length is
            # not a multiple of subscript length" from its internal indexing;
            # it is muffled here, only for this exact message, because it is
            # upstream and does not change the result under current R releases.
            roc_obj <- withCallingHandlers(
              timeROC::timeROC(
                T = Y, delta = delta, marker = marker,
                cause = 1, times = t_horizon, iid = TRUE
              ),
              warning = function(w) {
                if (identical(conditionMessage(w),
                              "object length is not a multiple of subscript length")) {
                  invokeRestart("muffleWarning")
                }
              }
            )
            auc_hat <- roc_obj$AUC[2]
            se_auc <- roc_obj$inference$vect_sd_1[2]

            if (is.na(auc_hat) || is.na(se_auc) || se_auc <= 0) {
              FALSE
            } else {
              ci_width <- 2 * stats::qnorm(0.975) * se_auc
              ci_width <= target_width
            }
          }
        }, error = function(e) FALSE)

        if (isTRUE(ok)) success_count <- success_count + 1L
      }

      prob <- success_count / B
      best_prob <- max(best_prob, prob)
      if (prob >= target_prob) {
        found_N <- as.integer(N)
        found_prob <- prob
        break
      }
    }

    if (is.na(found_N)) {
      warning("Nominal censoring rate ", cens_rate, ": target precision not ",
              "reached within N_range (best = ", round(best_prob, 3),
              "). Consider expanding N_range or increasing B.")
      # Report the largest N tried together with the assurance achieved
      # *at that N* (a consistent (N, prob) pair); best_prob over the grid
      # is conveyed in the warning above.
      found_N <- max(N_range)
      found_prob <- prob
    }

    results_list[[cr_idx]] <- data.frame(
      censoring_rate = cens_rate,
      N_required = found_N,
      prob_achieved = found_prob,
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  # Use the result for the middle censoring rate (or max) as primary
  N_primary <- max(results$N_required)
  n_events_expected <- floor(N_primary * (1 - exp(-lambda_event * t_horizon)))

  structure(
    list(
      method = "Sample Size for Time-Dependent ROC (AUC(t))",
      n_diseased = n_events_expected,
      n_total = N_primary,
      results = results,
      mu_case = mu_case,
      mu_control = mu_control,
      sigma = sigma,
      lambda_event = lambda_event,
      t_horizon = t_horizon,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
