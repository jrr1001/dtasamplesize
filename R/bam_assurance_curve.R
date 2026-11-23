#' Exact Assurance Curve A(N) for the Bayesian Assurance Method
#'
#' Computes the \strong{exact} assurance \eqn{A(N)} (no simulation, no
#' sampling error) at every total sample size \code{N} in \code{N_range},
#' with the same model, conventions and arguments as
#' \code{\link{bam_sample_size}} with \code{method = "exact"}. It is the
#' validated way to draw the assurance-versus-\code{N} curve (see
#' \code{\link{plot_assurance_curve}}) and to inspect how the assurance
#' behaves around the first \code{N} that reaches the target.
#'
#' @details \strong{Same calculation as \code{bam_sample_size()}.} For each
#'   \code{N} the function evaluates the closed form described in
#'   \code{\link{bam_sample_size}}: the number of diseased subjects is
#'   Beta-Binomial(\code{N}, \code{prior_prev}), and, given the split, the Se
#'   and Sp arms are independent, each with its own Beta-Binomial predictive
#'   count under its \emph{design} prior and a posterior width computed under
#'   its \emph{analysis} prior. The per-arm probabilities \eqn{P_{se}(n)} and
#'   \eqn{P_{sp}(n)} are built once, up to \code{max(N_range)}, and reused for
#'   every \code{N}; the empty-arm convention \eqn{P_{se}(0) = P_{sp}(0) = 0}
#'   is applied exactly as in \code{bam_sample_size()}. The internal helpers
#'   are the same functions, so \code{A(N)} is identical, bit for bit, to
#'   \code{bam_sample_size(..., method = "exact", N_range = N)$max_assurance_evaluated}
#'   at each \code{N}.
#'
#'   \strong{Conventions, stated explicitly.}
#'   \itemize{
#'     \item \emph{Width criterion} (\code{width_type}): \code{"full_width"}
#'       (default) is the full width of the \strong{equal-tailed two-sided}
#'       \code{alpha_ci} posterior credible interval; \code{"lower_half_width"}
#'       is the posterior median minus the \strong{one-sided} lower
#'       \code{alpha_ci} bound (convention of Williamson et al. 2023).
#'       \code{delta_se} and \code{delta_sp} are on that scale; widths are
#'       compared exactly (\code{width <= delta}).
#'     \item \emph{Objective} (\code{assurance_for}): \code{"joint"} (default)
#'       requires both widths in the same study; \code{"se"} or \code{"sp"}
#'       give the \strong{marginal} assurance of one arm (the \code{assurance}
#'       column is then a marginal assurance, not a joint one).
#'     \item \emph{Design and analysis priors}: \code{prior_se},
#'       \code{prior_sp}, \code{prior_prev} are design priors;
#'       \code{analysis_prior_se} and \code{analysis_prior_sp} (default
#'       \code{NULL}: same as the design prior) define the posterior. There is
#'       no analysis prior for prevalence. If \code{priors} is given, all five
#'       come from it and passing any of them as well is an error.
#'   }
#'
#'   \strong{Monotonicity.} The exact assurance is not guaranteed to be
#'   non-decreasing in \code{N}: because the posterior width is a discrete
#'   function of the count, small local decreases can occur. The returned
#'   object reports the number of decreases between consecutive \code{N_range}
#'   values (\code{n_decreases}) and the largest one (\code{max_decrease}).
#'   \code{N_first} is the \emph{first} \code{N} in \code{N_range} (ascending)
#'   with \eqn{A(N) \ge} \code{target_assurance}, as in
#'   \code{bam_sample_size()}; it says nothing about later values.
#'
#'   \strong{Cost.} The per-arm cache costs \eqn{O(\max(N)^2)} \code{qbeta}
#'   evaluations (a few seconds for \code{max(N_range)} around 700); a warning
#'   is issued above 3000.
#'
#' @param N_range Integer vector of total sample sizes, each \eqn{\ge 2}
#'   (\code{N = 1} is not evaluated: with the joint objective it always has
#'   assurance 0, and with a marginal objective it is not accepted either).
#'   Required. It
#'   is evaluated as \code{sort(unique(as.integer(N_range)))}; non-integer
#'   values are an error.
#' @param prior_se,prior_sp,prior_prev,delta_se,delta_sp,target_assurance,alpha_ci,width_type,assurance_for,analysis_prior_se,analysis_prior_sp,priors
#'   As in \code{\link{bam_sample_size}} (same defaults, same validation:
#'   \code{alpha_ci} must be a single number in (0, 1), and greater than 0.5
#'   with \code{width_type = "lower_half_width"}; the same warning when the
#'   effective analysis prior of a used arm already meets its target width
#'   with no data). \code{target_assurance} only defines the horizontal
#'   reference and \code{N_first}; it does not alter \code{A(N)}.
#' @return An object of class \code{"dta_bam_curve"}, a list with:
#'   \describe{
#'     \item{curve}{Data frame with columns \code{N} (integer, ascending) and
#'       \code{assurance} (exact \eqn{A(N)}).}
#'     \item{N_first}{First \code{N} with \eqn{A(N) \ge} \code{target_assurance};
#'       \code{NA_integer_} if none.}
#'     \item{assurance_at_N_first}{\eqn{A} at \code{N_first} (\code{NA_real_}
#'       if none).}
#'     \item{target_reached}{Logical.}
#'     \item{target_assurance}{The target used.}
#'     \item{n_decreases, max_decrease, monotone}{Number of decreases between
#'       consecutive points of the curve, the largest decrease (0 if none) and
#'       whether the curve is non-decreasing over \code{N_range}.}
#'     \item{width_type, alpha_ci, assurance_for, delta_se, delta_sp}{
#'       Conventions used.}
#'     \item{prior_se, prior_sp, prior_prev}{Design priors.}
#'     \item{analysis_prior_se, analysis_prior_sp, analysis_prior_given}{
#'       Effective analysis priors and whether each was given explicitly.}
#'     \item{assurance_method}{Always \code{"exact"}.}
#'     \item{call}{The matched call.}
#'   }
#' @references
#' Wilson KJ et al. (2022). Bayesian sample size determination for
#' diagnostic accuracy studies. \emph{Stat Med} 41:2908-2922.
#' \doi{10.1002/sim.9393}
#'
#' Williamson SF, Williams CJ, Lendrem BC, Wilson KJ (2023). Sample size
#' determination for point-of-care COVID-19 diagnostic tests: a Bayesian
#' approach. \emph{Diagn Progn Res} 7:17. \doi{10.1186/s41512-023-00153-1}
#' @seealso \code{\link{bam_sample_size}}, \code{\link{plot_assurance_curve}}
#' @examples
#' cv <- bam_assurance_curve(
#'   prior_se = c(17, 3), prior_sp = c(2, 2), prior_prev = c(4, 16),
#'   delta_se = 0.14, delta_sp = 0.10, N_range = 670:685
#' )
#' cv
#' cv$N_first
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_assurance_curve(cv)
#' }
#' }
#' @export
bam_assurance_curve <- function(prior_se = c(17, 3),
                                prior_sp = c(2, 2),
                                delta_se = 0.14,
                                delta_sp = 0.10,
                                target_assurance = 0.80,
                                prior_prev = c(6, 14),
                                N_range,
                                alpha_ci = 0.95,
                                width_type = c("full_width", "lower_half_width"),
                                assurance_for = c("joint", "se", "sp"),
                                analysis_prior_se = NULL,
                                analysis_prior_sp = NULL,
                                priors = NULL) {
  cl <- match.call()
  width_type <- tryCatch(
    match.arg(width_type),
    error = function(e) stop(
      "`width_type` must be one of \"full_width\" or \"lower_half_width\".",
      call. = FALSE
    )
  )
  assurance_for <- tryCatch(
    match.arg(assurance_for),
    error = function(e) stop(
      "`assurance_for` must be one of \"joint\", \"se\" or \"sp\".",
      call. = FALSE
    )
  )

  # --- optional "dta_bam_priors" object: same handling as bam_sample_size()
  if (!is.null(priors)) {
    if (!inherits(priors, "dta_bam_priors")) {
      stop("`priors` must be NULL or an object of class \"dta_bam_priors\" ",
           "(see bam_prior_from_pilot()).", call. = FALSE)
    }
    clash <- c(
      prior_se = !missing(prior_se), prior_sp = !missing(prior_sp),
      prior_prev = !missing(prior_prev),
      analysis_prior_se = !missing(analysis_prior_se),
      analysis_prior_sp = !missing(analysis_prior_sp)
    )
    if (any(clash)) {
      stop("`priors` supplies prior_se, prior_sp, prior_prev, ",
           "analysis_prior_se and analysis_prior_sp; do not also pass: ",
           paste(names(clash)[clash], collapse = ", "), ".", call. = FALSE)
    }
    prior_se <- priors$prior_se
    prior_sp <- priors$prior_sp
    prior_prev <- priors$prior_prev
    analysis_prior_se <- priors$analysis_prior_se
    analysis_prior_sp <- priors$analysis_prior_sp
  }

  # --- validation (same checks as bam_sample_size(), plus N_range, alpha_ci)
  if (missing(N_range)) {
    stop("`N_range` is required: an integer vector of total sample sizes, ",
         "each >= 2.", call. = FALSE)
  }
  stopifnot(length(prior_se) == 2, all(prior_se > 0))
  stopifnot(length(prior_sp) == 2, all(prior_sp > 0))
  stopifnot(length(prior_prev) == 2, all(prior_prev > 0))
  stopifnot(delta_se > 0, delta_sp > 0)
  stopifnot(target_assurance > 0, target_assurance < 1)
  .bam_check_analysis_prior(analysis_prior_se, "analysis_prior_se")
  .bam_check_analysis_prior(analysis_prior_sp, "analysis_prior_sp")
  .bam_check_alpha_ci(alpha_ci, width_type)
  if (!is.numeric(N_range) || length(N_range) < 1L || anyNA(N_range) ||
        any(!is.finite(N_range)) || any(N_range != round(N_range)) ||
        any(N_range < 2)) {
    stop("`N_range` must be a non-empty vector of whole numbers, all >= 2 ",
         "(N = 1 is not evaluated; with the joint objective it always ",
         "has assurance 0).", call. = FALSE)
  }

  a_se <- prior_se[1]
  b_se <- prior_se[2]
  a_sp <- prior_sp[1]
  b_sp <- prior_sp[2]

  analysis_given <- c(se = !is.null(analysis_prior_se),
                      sp = !is.null(analysis_prior_sp))
  analysis_prior_se_eff <- if (is.null(analysis_prior_se)) prior_se else analysis_prior_se
  analysis_prior_sp_eff <- if (is.null(analysis_prior_sp)) prior_sp else analysis_prior_sp
  aA_se <- analysis_prior_se_eff[1]
  bA_se <- analysis_prior_se_eff[2]
  aA_sp <- analysis_prior_sp_eff[1]
  bA_sp <- analysis_prior_sp_eff[2]

  ci_lower_q <- (1 - alpha_ci) / 2
  ci_upper_q <- 1 - ci_lower_q

  use_se_arm <- assurance_for != "sp"
  use_sp_arm <- assurance_for != "se"

  # Same safeguard as bam_sample_size(): the effective analysis prior of each
  # used arm already meets its target width with no data.
  .bam_warn_prior_only(prior_se, prior_sp, analysis_prior_se_eff,
                       analysis_prior_sp_eff, delta_se, delta_sp,
                       ci_lower_q, ci_upper_q, width_type,
                       use_se_arm, use_sp_arm)

  N_used <- sort(unique(as.integer(N_range)))
  N_max_exact <- max(N_used)
  if (N_max_exact > 3000L) {
    warning(
      "bam_assurance_curve() caches per-arm probabilities up to N = ",
      N_max_exact, ", which is O(N^2) qbeta evaluations and may be very ",
      "slow.", call. = FALSE
    )
  }

  # Per-arm caches, built once up to max(N_range), exactly as in the
  # user-N_range branch of bam_sample_size().
  P_se_exact <- NULL
  P_sp_exact <- NULL
  if (use_se_arm) {
    P_se_exact <- .bam_exact_width_prob(
      N_max_exact, aA_se, bA_se, delta_se, ci_lower_q, ci_upper_q,
      a_d = a_se, b_d = b_se, width_type = width_type
    )
  }
  if (use_sp_arm) {
    P_sp_exact <- .bam_exact_width_prob(
      N_max_exact, aA_sp, bA_sp, delta_sp, ci_lower_q, ci_upper_q,
      a_d = a_sp, b_d = b_sp, width_type = width_type
    )
  }
  # Empty-arm convention: an arm of size 0 is always a failure.
  if (use_se_arm) P_se_exact[1] <- 0
  if (use_sp_arm) P_sp_exact[1] <- 0

  A <- numeric(length(N_used))
  for (i in seq_along(N_used)) {
    A[i] <- .bam_exact_joint_assurance(
      N_used[i], prior_prev[1], prior_prev[2], P_se_exact, P_sp_exact,
      assurance_for = assurance_for
    )
  }

  hit <- which(A >= target_assurance)
  target_reached <- length(hit) > 0L
  N_first <- if (target_reached) N_used[hit[1]] else NA_integer_
  A_first <- if (target_reached) A[hit[1]] else NA_real_

  d <- diff(A)
  n_decreases <- sum(d < 0)
  max_decrease <- max(0, -d)

  structure(
    list(
      curve = data.frame(N = N_used, assurance = A),
      N_first = N_first,
      assurance_at_N_first = A_first,
      target_reached = target_reached,
      target_assurance = target_assurance,
      n_decreases = n_decreases,
      max_decrease = max_decrease,
      monotone = n_decreases == 0L,
      width_type = width_type,
      alpha_ci = alpha_ci,
      assurance_for = assurance_for,
      delta_se = delta_se,
      delta_sp = delta_sp,
      prior_se = prior_se,
      prior_sp = prior_sp,
      prior_prev = prior_prev,
      analysis_prior_se = analysis_prior_se_eff,
      analysis_prior_sp = analysis_prior_sp_eff,
      analysis_prior_given = analysis_given,
      assurance_method = "exact",
      call = cl
    ),
    class = "dta_bam_curve"
  )
}

# Short description of the objective (used in print() and plot labels).
.bam_curve_objective <- function(assurance_for) {
  switch(assurance_for, joint = "joint assurance",
         se = "marginal Se assurance", sp = "marginal Sp assurance")
}

#' Print Method for dta_bam_curve Objects
#'
#' @param x An object of class \code{"dta_bam_curve"}, as returned by
#'   \code{\link{bam_assurance_curve}}.
#' @param ... Currently unused.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{bam_assurance_curve}}
#' @export
print.dta_bam_curve <- function(x, ...) {
  fmt <- function(ab) paste0("Beta(", format(ab[1]), ", ", format(ab[2]), ")")
  cat("Exact assurance curve A(N) (Bayesian Assurance Method)\n")
  cat("------------------------------------------------------\n")
  cat("Objective      :", .bam_curve_objective(x$assurance_for), "\n")
  cat("Width criterion:", .bam_width_label(x$width_type, x$alpha_ci), "\n")
  cat("Targets (delta): Se ", x$delta_se, ", Sp ", x$delta_sp, "\n", sep = "")
  cat("Design priors  : Se ", fmt(x$prior_se), ", Sp ", fmt(x$prior_sp),
      ", prevalence ", fmt(x$prior_prev), "\n", sep = "")
  # Compare VALUES: an analysis prior that was given but equals the design
  # prior (explicitly, or through a `priors` object with use_in = "both") is
  # "same as design".
  ana <- function(ab, db) {
    paste0(fmt(ab), if (identical(as.numeric(ab), as.numeric(db)))
      " (same as design)" else " (given)")
  }
  cat("Analysis priors: Se ", ana(x$analysis_prior_se, x$prior_se), ", Sp ",
      ana(x$analysis_prior_sp, x$prior_sp), "\n", sep = "")
  cat("N evaluated    : ", nrow(x$curve), " values, ", min(x$curve$N), " to ",
      max(x$curve$N), "\n", sep = "")
  cat("Target         : A(N) >= ", x$target_assurance, "\n", sep = "")
  if (isTRUE(x$target_reached)) {
    cat("First N        : ", x$N_first, " (A = ",
        sprintf("%.6f", x$assurance_at_N_first), ")\n", sep = "")
  } else {
    cat("First N        : none in the range (max A = ",
        sprintf("%.6f", max(x$curve$assurance)), ")\n", sep = "")
  }
  if (x$monotone) {
    cat("Monotone       : yes (non-decreasing over the range)\n")
  } else {
    cat("Monotone       : no (", x$n_decreases, " decrease(s); largest ",
        sprintf("%.6f", x$max_decrease), ")\n", sep = "")
  }
  invisible(x)
}
