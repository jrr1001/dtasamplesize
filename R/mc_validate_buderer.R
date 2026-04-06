#' Monte Carlo Validation of the Buderer Sample Size Formula
#'
#' Generates simulated datasets under the exact assumptions of the
#' Buderer formula and evaluates how often the resulting confidence
#' interval achieves the planned precision.
#'
#' @param Se Expected sensitivity (0 < Se < 1). Default 0.85.
#' @param d Desired precision (half-width of CI). Default 0.07.
#' @param n_diseased Number of diseased individuals. If \code{NULL}
#'   (default), computed from the Buderer formula.
#' @param B Number of Monte Carlo replications. Default 5000.
#' @param ci_method CI method: \code{"wald"}, \code{"wilson"}, or
#'   \code{"both"} (default).
#' @param alpha Significance level. Default 0.05.
#' @param seed Random seed for reproducibility. Default 2026.
#' @return Object of class \code{"dtasamplesize"} with:
#'   \describe{
#'     \item{method}{Character string identifying the method.}
#'     \item{n_diseased}{Sample size used.}
#'     \item{n_total}{Same as \code{n_diseased} (no prevalence adjustment).}
#'     \item{Se}{Expected sensitivity used.}
#'     \item{d}{Desired precision used.}
#'     \item{results}{\code{data.frame} with columns: \code{ci_method},
#'       \code{P_width_target}, \code{coverage}, \code{median_width},
#'       \code{P75_width}, \code{P90_width}.}
#'   }
#' @references
#' Buderer NMF (1996). Statistical methodology: I. Incorporating the
#' prevalence of disease into the sample size calculation for sensitivity
#' and specificity. \emph{Acad Emerg Med} 3:895-900.
#' \doi{10.1111/j.1553-2712.1996.tb03538.x}
#'
#' Newcombe RG (1998). Two-sided confidence intervals for the single
#' proportion: comparison of seven methods. \emph{Stat Med} 17:857-872.
#' \doi{10.1002/(SICI)1097-0258(19980430)17:8<857::AID-SIM777>3.0.CO;2-E}
#' @examples
#' result <- mc_validate_buderer(Se = 0.85, d = 0.07, B = 1000)
#' print(result)
#' @export
mc_validate_buderer <- function(Se = 0.85,
                                d = 0.07,
                                n_diseased = NULL,
                                B = 5000,
                                ci_method = "both",
                                alpha = 0.05,
                                seed = 2026) {
  # Validate inputs
  stopifnot(Se > 0, Se < 1)
  stopifnot(d > 0, d < 0.5)
  stopifnot(B >= 1)
  stopifnot(ci_method %in% c("wald", "wilson", "both"))
  stopifnot(alpha > 0, alpha < 1)


  # Compute n_diseased from Buderer if not supplied
  if (is.null(n_diseased)) {
    n_diseased <- buderer_n(Se, d, alpha)
  }

  target_width <- 2 * d
  z <- stats::qnorm(1 - alpha / 2)

  # Simulate all B binomial samples at once
  set.seed(seed)
  x <- stats::rbinom(B, n_diseased, Se)
  se_hat <- x / n_diseased

  methods <- if (ci_method == "both") c("wald", "wilson") else ci_method
  results_list <- list()

  for (m in methods) {
    if (m == "wald") {
      se_se <- sqrt(se_hat * (1 - se_hat) / n_diseased)
      ci_lower <- pmax(0, se_hat - z * se_se)
      ci_upper <- pmin(1, se_hat + z * se_se)
    } else {
      # Wilson score (vectorized)
      denom <- 1 + z^2 / n_diseased
      center <- (se_hat + z^2 / (2 * n_diseased)) / denom
      margin <- z * sqrt((se_hat * (1 - se_hat) +
        z^2 / (4 * n_diseased)) / n_diseased) / denom
      ci_lower <- pmax(0, center - margin)
      ci_upper <- pmin(1, center + margin)
    }

    ci_width <- ci_upper - ci_lower
    coverage <- mean(ci_lower <= Se & Se <= ci_upper)
    P_width_target <- mean(ci_width <= target_width)

    results_list[[m]] <- data.frame(
      ci_method = m,
      P_width_target = P_width_target,
      coverage = coverage,
      median_width = stats::median(ci_width),
      P75_width = stats::quantile(ci_width, 0.75, names = FALSE),
      P90_width = stats::quantile(ci_width, 0.90, names = FALSE),
      stringsAsFactors = FALSE
    )
  }

  results <- do.call(rbind, results_list)
  rownames(results) <- NULL

  structure(
    list(
      method = "Monte Carlo Validation of Buderer Formula",
      n_diseased = n_diseased,
      n_total = n_diseased,
      Se = Se,
      d = d,
      B = B,
      alpha = alpha,
      results = results,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
