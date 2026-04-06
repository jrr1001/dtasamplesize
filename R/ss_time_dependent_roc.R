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
#' @param censoring_rates Numeric vector of censoring proportions.
#'   Default \code{c(0.10, 0.20, 0.30)}.
#' @param delta_auc Target precision (half-width) for AUC(t). Default 0.06.
#' @param target_prob Probability of achieving precision. Default 0.80.
#' @param N_range Range of N to search. Default \code{seq(100, 500, by = 20)}.
#' @param B MC replications. Default 500 (lower due to timeROC cost).
#' @param seed Random seed. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{results}{Data frame with columns \code{censoring_rate},
#'       \code{N_required}, \code{prob_achieved}.}
#'   }
#' @note Parameters in the default example are HYPOTHETICAL. No published
#'   AUC values exist for ctDNA as a continuous discriminator in DLBCL.
#'   Requires the \pkg{timeROC} package (listed in Suggests).
#' @examples
#' \donttest{
#' result <- ss_time_dependent_roc(B = 50, N_range = seq(100, 300, by = 50))
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

    for (N in N_range) {
      set.seed(seed)
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
            roc_obj <- timeROC::timeROC(
              T = Y, delta = delta, marker = marker,
              cause = 1, times = t_horizon, iid = TRUE
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
