# --- Internal helpers for the exact (non-Monte-Carlo) joint search --------
# Not exported. Used only by bam_sample_size() when method = "exact"; see
# that function's @details for the closed-form derivation these implement.

# For a single arm governed by a Beta(a, b) prior and a full-width credible
# interval target of delta, return a vector P of length n_max + 1 where
# P[n + 1] = P(width <= delta | arm size = n), marginalizing over the arm's
# own success count x ~ BetaBinomial(n, a, b). The posterior credible
# interval width after x successes in n trials is a deterministic function
# of (x, n) alone (Beta(a + x, b + n - x)), so this sum is a finite,
# non-random enumeration -- no simulation involved.
#
# This depends only on (n_max, a, b, delta), never on N or on which
# candidate N is being evaluated, so bam_sample_size() calls it exactly
# once per arm (up to max(N_range)) and reuses the result across every
# candidate N in the search -- the caching the exact mode relies on to stay
# fast despite its O(N^2) cost.
.bam_exact_width_prob <- function(n_max, a, b, delta, ci_lower_q, ci_upper_q) {
  P <- numeric(n_max + 1L)
  lbeta_ab <- lbeta(a, b)
  for (n in 0:n_max) {
    x <- 0:n
    post_a <- a + x
    post_b <- b + n - x
    width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    log_pmf <- lchoose(n, x) + lbeta(post_a, post_b) - lbeta_ab
    P[n + 1L] <- sum(exp(log_pmf)[width <= delta])
  }
  P
}

# Exact joint assurance at a single total N: the number of diseased subjects
# n_d is Beta-Binomial(N, a_p, b_p), and, conditional on the resulting split
# (n_d, N - n_d), the Se and Sp arms are independent (each depends only on
# its own true-value prior and its own arm size). P_se and P_sp must already
# hold P_arm(n) at index n + 1 for every n from 0 to at least N (see
# .bam_exact_width_prob() above); this function performs no qbeta calls of
# its own, so evaluating it for many candidate N once the cache exists is
# cheap (O(N) per candidate).
.bam_exact_joint_assurance <- function(N, a_p, b_p, P_se, P_sp) {
  k <- 0:N
  log_pmf <- lchoose(N, k) + lbeta(a_p + k, b_p + N - k) - lbeta(a_p, b_p)
  p_nd <- exp(log_pmf)
  sum(p_nd * P_se[k + 1L] * P_sp[N - k + 1L])
}

#' Bayesian Assurance Method for DTA Sample Size
#'
#' Finds the minimum total sample size at which the posterior credible
#' intervals for sensitivity AND specificity \strong{simultaneously} reach
#' their target width, with probability at least \code{target_assurance},
#' under a single cohort study of size \code{N}.
#'
#' @details The parameters \code{delta_se} and \code{delta_sp} are the
#'   \strong{full width} of the posterior credible interval (not half-width).
#'   For example, \code{delta_se = 0.14} corresponds to a half-width of 0.07.
#'   This differs from other functions in this package that use half-width.
#'
#'   \strong{The headline result (\code{N_total}) is a JOINT assurance under
#'   cohort design.} At each candidate total \code{N}, every one of the
#'   \code{B} Monte Carlo replications draws \code{prev ~ Beta(prior_prev)},
#'   \code{Se_true ~ Beta(prior_se)} and \code{Sp_true ~ Beta(prior_sp)}, and
#'   then, because the number of diseased subjects in a fixed-\code{N} cohort
#'   study is itself random, draws \code{n_d ~ Binomial(N, prev)} and sets
#'   \code{n_nd = N - n_d} (the same cohort convention used by
#'   \code{\link{ss_net_benefit}} and \code{\link{ss_unified}}). Counts
#'   \code{x_se ~ Binomial(n_d, Se_true)} and
#'   \code{x_sp ~ Binomial(n_nd, Sp_true)} are then drawn, and the posterior
#'   Beta credible interval width is computed for each arm using the same
#'   \code{alpha_ci} and the same full-width convention described above. A
#'   replication is a \strong{success} only if \strong{both} widths meet
#'   their targets \strong{in that same replication}; \code{joint_assurance}
#'   is the proportion of successes, with denominator always \code{B}. A
#'   degenerate replication (\code{n_d = 0} or \code{n_nd = 0}, possible when
#'   \code{N} is small or the prevalence draw is extreme) is always counted
#'   as a \strong{failure}, never dropped from the denominator -- the same
#'   convention already adopted by \code{\link{ss_unified}}.
#'
#'   This corrects a defect present through version 0.4.0, where \code{n_se}
#'   and \code{n_sp} were found \strong{independently}, each against its own
#'   \strong{marginal} assurance, and the previous \code{n_total} (now
#'   \code{N_total_median}, see below) was only the median of
#'   \code{pmax(n_se / prev, n_sp / (1 - prev))} over draws of \code{prev}.
#'   That quantity never verified that both widths hold in the same study,
#'   and its joint assurance is systematically lower than
#'   \code{target_assurance} -- e.g., for the defaults with
#'   \code{prior_prev = c(4, 16)}, the old \code{n_total} achieved a true
#'   joint assurance of about 0.72 against a target of 0.80.
#'
#'   \strong{Two ways to evaluate the joint assurance: \code{method}.} The
#'   generative model above -- draw \code{prev}, \code{Se_true}, \code{Sp_true};
#'   split \code{N} into \code{n_d} and \code{n_nd}; draw the two arm counts;
#'   check both widths -- has a closed form, because \code{n_d} is
#'   Beta-Binomial given \code{prior_prev}, and, given an arm size, the arm's
#'   own success count is Beta-Binomial given its own prior, with a
#'   deterministic (hence exactly enumerable) credible-interval width for
#'   every possible count. \code{method = "exact"} (the default) evaluates
#'   this closed form directly, with no simulation and hence no sampling
#'   noise: for a candidate \code{N},
#'   \deqn{
#'     \mathrm{assurance}(N) = \sum_{k=0}^{N} P(n_d = k) \, P_{se}(k) \, P_{sp}(N - k),
#'   }
#'   where \code{P(n_d = k)} is the Beta-Binomial(\code{N}, \code{prior_prev})
#'   mass at \code{k}, and, for an arm of size \code{n},
#'   \deqn{
#'     P_{se}(n) = \sum_{x=0}^{n} \mathrm{BetaBinom}(x; n, a_{se}, b_{se}) \,
#'       \mathbf{1}[w_{se}(x, n) \le \delta_{se}],
#'   }
#'   with \code{w_se(x, n)} the posterior credible interval width after
#'   \code{x} successes in \code{n} trials (the same width formula and the
#'   same \code{alpha_ci} used everywhere else in this function), and
#'   \code{P_sp(n)} defined analogously from \code{prior_sp} and
#'   \code{delta_sp}. By convention \code{P_se(0) = P_sp(0) = 0}: an arm of
#'   size zero (\code{n_d = 0} or \code{n_nd = 0}) is the same degenerate
#'   replication described above and is always a failure, regardless of how
#'   narrow the prior-only credible interval happens to be, matching
#'   \code{method = "monte_carlo"} exactly. The two arms are conditionally
#'   independent given \code{(n_d, n_nd)}, which licenses the product
#'   \code{P_se(k) * P_sp(N - k)} inside the sum. \code{N} is accepted as
#'   soon as this exact
#'   probability reaches \code{target_assurance}; there is no sampling error
#'   to guard against, so \code{B} plays no role in this calculation (it is
#'   still used for the legacy per-arm searches and heuristics described
#'   below) and \code{seed} does not affect the result either. This is the
#'   preferred mode: it is deterministic, reproducible bit-for-bit, and free
#'   of the dependence on \code{B} that affects \code{method = "monte_carlo"}
#'   (e.g., for the defaults with \code{prior_prev = c(4, 16)}, the exact
#'   first \code{N} with joint assurance \eqn{\ge} 0.80 is 678, with assurance
#'   0.800349, whereas \code{method = "monte_carlo"} with \code{B = 20000}
#'   returns 683 because of Monte Carlo noise near the crossing point -- see
#'   \strong{Anti-noise acceptance rule} below). \code{P_se} and \code{P_sp}
#'   depend only on arm size, not on \code{N}, so they are cached once, up to
#'   \code{max(N_range)}, and reused across every candidate \code{N}; the
#'   cache itself costs \eqn{O(N_{\max}^2)} \code{qbeta} evaluations (about
#'   500,000 for \code{N_max} around 700, a few seconds). If
#'   \code{max(N_range)} is large enough that this becomes impractical, a
#'   warning suggests \code{method = "monte_carlo"} instead.
#'
#'   \code{method = "monte_carlo"} instead estimates the same probability by
#'   simulation, as described above, and is kept for continuity with
#'   versions \verb{<= 0.4.x} of this joint search (introduced mid-cycle in
#'   0.5.0) and as a fallback for \code{N} ranges too large for the exact
#'   cache.
#'
#'   \strong{Anti-noise acceptance rule (\code{method = "monte_carlo"} only).}
#'   Because \code{joint_assurance} at each \code{N} is itself a Monte Carlo
#'   estimate with sampling error, an \code{N} is accepted only if the
#'   \strong{lower bound of its one-sided 95\% confidence interval} also
#'   reaches the target:
#'   \code{joint_assurance - qnorm(0.95) * assurance_mcse >= target_assurance}
#'   (where \code{assurance_mcse = sqrt(joint_assurance * (1 - joint_assurance) / B)}
#'   is the Monte Carlo standard error of a proportion estimated from
#'   \code{B} replications). Without this rule, the search could accept the
#'   first \code{N} where the Monte Carlo estimate happens to cross
#'   \code{target_assurance} by chance, on a noisy grid, rather than an
#'   \code{N} where the target is genuinely reached. This rule makes the
#'   search slightly conservative (it can select an \code{N} a few grid
#'   steps above the raw crossing point), especially at small \code{B}. This
#'   rule does not apply under \code{method = "exact"}, which has no sampling
#'   error to guard against: the first \code{N} whose exact assurance reaches
#'   \code{target_assurance} is accepted outright.
#'
#'   \strong{\code{n_range} vs. \code{N_range}.} \code{n_range} keeps its
#'   original meaning: the grid of \strong{per-arm} sample sizes searched
#'   independently for the diagnostic (legacy) fields \code{n_diseased},
#'   \code{n_non_diseased}, \code{assurance_se}, and \code{assurance_sp}. It
#'   no longer determines the headline result. \code{N_range} is a new
#'   argument that controls the grid of \strong{total} sample sizes searched
#'   for \code{N_total}, the joint-assurance result described above. When
#'   \code{N_range = NULL} (the default), it is built automatically as
#'   \code{max(20, n_se, n_sp):ceiling(N_total_P90)}, i.e. from the smallest
#'   total that could possibly supply both per-arm requirements up to the
#'   90th percentile of the legacy prevalence-uncertainty heuristic (see
#'   \code{N_total_P90} below), which in practice comfortably brackets the
#'   true joint requirement. No existing argument name was removed or had
#'   its default behavior changed; \code{N_range} is purely additive.
#'
#' @param prior_se Numeric vector \code{c(alpha, beta)} for Beta prior on
#'   sensitivity. Default \code{c(17, 3)} (E[Se]=0.85, ~20 pseudo-observations).
#' @param prior_sp Numeric vector \code{c(alpha, beta)} for Beta prior on
#'   specificity. Default \code{c(2, 2)} (vague).
#' @param delta_se Target full width for Se credible interval. Default 0.14.
#' @param delta_sp Target full width for Sp credible interval. Default 0.10.
#' @param target_assurance Minimum assurance probability. Default 0.80.
#' @param prior_prev Numeric vector \code{c(alpha, beta)} for Beta prior on
#'   prevalence. Default \code{c(6, 14)} (E[prev]=0.30).
#' @param n_range Integer vector of candidate \strong{per-arm} n values
#'   searched independently for the legacy/diagnostic fields
#'   \code{n_diseased}, \code{n_non_diseased}, \code{assurance_se}, and
#'   \code{assurance_sp}. Default \code{20:500}. See \code{Details} for how
#'   this differs from \code{N_range}.
#' @param N_range Integer vector of candidate \strong{total} N values
#'   searched for the headline joint-assurance result \code{N_total}.
#'   Default \code{NULL}, which builds the grid automatically from
#'   \code{n_se}, \code{n_sp}, and \code{N_total_P90} -- see \code{Details}.
#' @param B Number of MC replications. Default 5000. Used for the legacy
#'   per-arm searches and heuristics (\code{n_diseased}, \code{n_non_diseased},
#'   \code{assurance_se}, \code{assurance_sp}, \code{N_total_median},
#'   \code{N_total_P75}, \code{N_total_P90}) regardless of \code{method}. Also
#'   governs the headline joint search (\code{N_total}, \code{joint_assurance},
#'   \code{assurance_mcse}) when \code{method = "monte_carlo"}, but
#'   \strong{ignored by that headline joint search when \code{method =
#'   "exact"}} (the default), since that calculation has no Monte Carlo error
#'   to control -- see \code{Details}. A warning is issued when
#'   \code{0 < B < 1000}, since Monte Carlo error may then be substantial in
#'   whichever results \code{B} actually affects: under \code{method =
#'   "monte_carlo"} this is the package-wide small-\code{B} warning (it
#'   applies to the headline result too); under \code{method = "exact"} the
#'   warning is reworded to name only the legacy fields above, since the
#'   headline result has no Monte Carlo error to warn about there. Either
#'   version is silenced with \code{options(dtasamplesize.warn_small_B =
#'   FALSE)}.
#' @param alpha_ci Credible interval level. Default 0.95.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit. Affects the legacy per-arm searches and
#'   heuristics regardless of \code{method}, but has \strong{no effect on
#'   \code{N_total} / \code{joint_assurance} when \code{method = "exact"}}
#'   (the default), since that calculation is deterministic.
#' @param method Either \code{"exact"} (the default) or \code{"monte_carlo"}.
#'   Selects how the headline joint search (\code{N_total},
#'   \code{joint_assurance}, \code{assurance_mcse}) is computed; see
#'   \code{Details}. \code{"exact"} evaluates the joint assurance in closed
#'   form via Beta-Binomial enumeration, with no simulation and no sampling
#'   error. \code{"monte_carlo"} reproduces the simulation-based search
#'   (including the anti-noise acceptance rule) unchanged.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{n_diseased}{Minimum diseased sample for the \strong{marginal}
#'       Se assurance (legacy diagnostic field; searched over \code{n_range}).}
#'     \item{n_non_diseased}{Minimum non-diseased sample for the
#'       \strong{marginal} Sp assurance (legacy diagnostic field; searched
#'       over \code{n_range}).}
#'     \item{n_total}{Alias of \code{N_total} (kept so that code written
#'       against the printed/generic \code{n_total} field, used across this
#'       package's \code{"dtasamplesize"} objects, keeps working). Note this
#'       is a \strong{behavior change} from versions <= 0.4.0, where
#'       \code{n_total} held \code{N_total_median}: see \code{Details}.}
#'     \item{N_total}{The minimum total N whose \strong{joint} assurance
#'       (Se and Sp both within target width, in the same replication, under
#'       cohort design) reaches \code{target_assurance}. Under \code{method =
#'       "exact"} this is the exact probability itself; under \code{method =
#'       "monte_carlo"} it is accepted at the 95\% Monte Carlo confidence
#'       lower bound (the anti-noise rule). This is the headline result of
#'       this function; see \code{Details}.}
#'     \item{joint_assurance}{The joint assurance achieved at \code{N_total}
#'       (point estimate, not the lower confidence bound used to accept it
#'       under \code{method = "monte_carlo"}; under \code{method = "exact"}
#'       this already \emph{is} the exact probability, with no further bound
#'       to distinguish it from).}
#'     \item{assurance_mcse}{Monte Carlo standard error of
#'       \code{joint_assurance}, i.e.
#'       \code{sqrt(joint_assurance * (1 - joint_assurance) / B)}, under
#'       \code{method = "monte_carlo"}. Always exactly \code{0} under
#'       \code{method = "exact"}, chosen over \code{NA} because
#'       \code{joint_assurance} under that method genuinely has zero sampling
#'       error (it is not "unknown"), and because \code{0} keeps downstream
#'       arithmetic such as the anti-noise formula well-defined without a
#'       special case.}
#'     \item{assurance_method}{The \code{method} argument actually used
#'       (\code{"exact"} or \code{"monte_carlo"}). Not named \code{method}
#'       because that field is already used, across every
#'       \code{"dtasamplesize"} object in this package (see
#'       \code{print.dtasamplesize}), to hold the human-readable method
#'       \emph{name} printed as this object's header (here, always
#'       \code{"Bayesian Assurance Method (BAM) for DTA Sample Size"},
#'       regardless of \code{method}); reusing it for the exact/Monte Carlo
#'       switch would silently break that shared convention.}
#'     \item{assurance_se}{Achieved \strong{marginal} assurance for Se at
#'       \code{n_diseased} (legacy diagnostic field).}
#'     \item{assurance_sp}{Achieved \strong{marginal} assurance for Sp at
#'       \code{n_non_diseased} (legacy diagnostic field).}
#'     \item{N_total_median}{\strong{Deprecated heuristic}, kept for backward
#'       compatibility with versions <= 0.4.0: the median of
#'       \code{pmax(n_se / prev, n_sp / (1 - prev))} over draws of
#'       \code{prev ~ Beta(prior_prev)}, where \code{n_se} and \code{n_sp}
#'       are the independent marginal requirements above. This quantity
#'       does \strong{not} guarantee a joint assurance of
#'       \code{target_assurance}; use \code{N_total} instead. Retained
#'       unchanged so that existing code and reports keep working.}
#'     \item{N_total_P75}{75th percentile of the same legacy heuristic
#'       distribution. Also used to size the default \code{N_range}.}
#'     \item{N_total_P90}{90th percentile of the same legacy heuristic
#'       distribution. Also used to size the default \code{N_range}.}
#'     \item{buderer_n_se}{Buderer sample size for comparison.}
#'   }
#' @references
#' Wilson KJ et al. (2022). Bayesian sample size determination for
#' diagnostic accuracy studies. \emph{Stat Med} 41:2908-2922.
#' \doi{10.1002/sim.9393}
#' @examples
#' # B and N_range kept small here for a fast example; see @param B.
#' result <- suppressWarnings(bam_sample_size(
#'   B = 500, n_range = 20:500, N_range = seq(100, 700, by = 20)
#' ))
#' print(result)
#' @export
bam_sample_size <- function(prior_se = c(17, 3),
                            prior_sp = c(2, 2),
                            delta_se = 0.14,
                            delta_sp = 0.10,
                            target_assurance = 0.80,
                            prior_prev = c(6, 14),
                            n_range = 20:500,
                            B = 5000,
                            alpha_ci = 0.95,
                            seed = 2026,
                            N_range = NULL,
                            method = c("exact", "monte_carlo")) {
  method <- match.arg(method)

  # --- preserve the caller's RNG state (kind AND seed) ------------------
  # Restoring only .Random.seed's VALUE is not enough: set.seed() called
  # later by unrelated code with no explicit `kind` argument reuses
  # whichever kind is CURRENTLY ACTIVE, not whatever kind .Random.seed
  # happens to encode. Every set.seed() call below therefore names its kind
  # explicitly (Mersenne-Twister, R's own default) rather than inheriting
  # whatever generator the caller happened to have active, so the legacy
  # per-arm searches and heuristics reproduce the same numbers regardless
  # of the caller's own RNG configuration; see save_rng_state().
  old_rng_state <- save_rng_state()
  on.exit(restore_rng_state(old_rng_state), add = TRUE)

  # Validate inputs
  stopifnot(length(prior_se) == 2, all(prior_se > 0))
  stopifnot(length(prior_sp) == 2, all(prior_sp > 0))
  stopifnot(length(prior_prev) == 2, all(prior_prev > 0))
  stopifnot(delta_se > 0, delta_sp > 0)
  stopifnot(target_assurance > 0, target_assurance < 1)
  stopifnot(B >= 1)
  stopifnot(is.null(N_range) || (length(N_range) >= 1 && all(N_range > 0)))

  # Small-B advisory, scoped to what B actually governs under each mode. B
  # always affects the legacy per-arm fields (n_diseased, n_non_diseased,
  # assurance_se, assurance_sp, N_total_median, N_total_P75, N_total_P90),
  # regardless of method -- see the searches below. Under method =
  # "monte_carlo" it ALSO governs the headline joint search (N_total /
  # joint_assurance / assurance_mcse), so the package-wide warn_small_B()
  # message (which speaks generically of "the assurance estimate") applies
  # as-is. Under method = "exact" the headline result has no Monte Carlo
  # error at all -- B plays no role in it (see @details) -- so reusing that
  # generic message would wrongly imply sampling error in a result that has
  # none; this reworded version names only the legacy fields it actually
  # affects.
  if (method == "monte_carlo") {
    warn_small_B(B)
  } else if (B > 0 && B < 1000 &&
               isTRUE(getOption("dtasamplesize.warn_small_B", TRUE))) {
    warning(
      "B = ", B, " is small. Under method = \"exact\" this affects only ",
      "the legacy per-arm diagnostic fields (n_diseased, n_non_diseased, ",
      "assurance_se, assurance_sp, N_total_median, N_total_P75, ",
      "N_total_P90), where Monte Carlo error may be substantial. The ",
      "headline N_total / joint_assurance are unaffected: method = ",
      "\"exact\" computes them in closed form, with no Monte Carlo error ",
      "and no dependence on B. B >= 1000 (ideally 5000) is recommended if ",
      "you rely on the legacy fields. Set options(dtasamplesize.warn_small_B ",
      "= FALSE) to silence this.",
      call. = FALSE
    )
  }

  a_se <- prior_se[1]
  b_se <- prior_se[2]
  a_sp <- prior_sp[1]
  b_sp <- prior_sp[2]

  ci_lower_q <- (1 - alpha_ci) / 2
  ci_upper_q <- 1 - ci_lower_q

  # --- Search for n_se (legacy marginal / diagnostic field) ---
  n_se <- NA_integer_
  assurance_se_achieved <- NA_real_

  for (n in n_range) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    se_true <- stats::rbeta(B, a_se, b_se)
    x <- stats::rbinom(B, n, se_true)
    # Posterior: Beta(a_se + x, b_se + n - x)
    post_a <- a_se + x
    post_b <- b_se + n - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance <- mean(ci_width <= delta_se)
    if (assurance >= target_assurance) {
      n_se <- as.integer(n)
      assurance_se_achieved <- assurance
      break
    }
  }

  if (is.na(n_se)) {
    warning("No n in n_range achieved target assurance for Se. ",
            "Consider expanding n_range.")
    n_se <- max(n_range)
    # Compute assurance at max
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    se_true <- stats::rbeta(B, a_se, b_se)
    x <- stats::rbinom(B, n_se, se_true)
    post_a <- a_se + x
    post_b <- b_se + n_se - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance_se_achieved <- mean(ci_width <= delta_se)
  }

  # --- Search for n_sp (legacy marginal / diagnostic field) ---
  n_sp <- NA_integer_
  assurance_sp_achieved <- NA_real_

  for (n in n_range) {
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    sp_true <- stats::rbeta(B, a_sp, b_sp)
    x <- stats::rbinom(B, n, sp_true)
    post_a <- a_sp + x
    post_b <- b_sp + n - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance <- mean(ci_width <= delta_sp)
    if (assurance >= target_assurance) {
      n_sp <- as.integer(n)
      assurance_sp_achieved <- assurance
      break
    }
  }

  if (is.na(n_sp)) {
    warning("No n in n_range achieved target assurance for Sp. ",
            "Consider expanding n_range.")
    n_sp <- max(n_range)
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
    sp_true <- stats::rbeta(B, a_sp, b_sp)
    x <- stats::rbinom(B, n_sp, sp_true)
    post_a <- a_sp + x
    post_b <- b_sp + n_sp - x
    ci_width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    assurance_sp_achieved <- mean(ci_width <= delta_sp)
  }

  # --- Legacy heuristic: total N distribution accounting for prevalence
  # uncertainty only, via the INDEPENDENT marginal requirements n_se / n_sp.
  # Kept unchanged (byte-identical) for backward compatibility -- see
  # @details and N_total_median in @return. This is NOT the joint-assurance
  # result; it is used below only to size the default N_range.
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  prev_draws <- stats::rbeta(B, prior_prev[1], prior_prev[2])
  N_total_draws <- pmax(n_se / prev_draws, n_sp / (1 - prev_draws))
  # Guard against prevalence draws at the 0/1 boundary (possible with very
  # diffuse priors), which would make a draw non-finite and contaminate the
  # upper quantiles.
  N_total_draws <- N_total_draws[is.finite(N_total_draws)]
  N_total_median <- stats::median(N_total_draws)
  N_total_P75 <- stats::quantile(N_total_draws, 0.75, names = FALSE)
  N_total_P90 <- stats::quantile(N_total_draws, 0.90, names = FALSE)

  # --- Joint search for total N (headline result, v0.5.0) -----------------
  # Unlike the two searches above, this one evaluates BOTH arms within the
  # SAME replication of a single cohort of size N, with n_d ~ Binomial(N,
  # prev) (the disease count is random, not fixed at floor(N * prev)). See
  # @details for the full generative model and the anti-noise acceptance
  # rule.
  if (is.null(N_range)) {
    N_lo <- max(20L, n_se, n_sp)
    N_hi <- max(N_lo + 1L, ceiling(N_total_P90))
    N_range_used <- N_lo:N_hi
  } else {
    N_range_used <- N_range
  }

  N_total <- NA_integer_
  joint_assurance_achieved <- NA_real_
  assurance_mcse_achieved <- NA_real_

  if (method == "monte_carlo") {
    z_mcse <- stats::qnorm(0.95)

    joint_assurance_last <- NA_real_
    assurance_mcse_last <- NA_real_

    for (N in N_range_used) {
      set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
      prev_b <- stats::rbeta(B, prior_prev[1], prior_prev[2])
      se_true_b <- stats::rbeta(B, a_se, b_se)
      sp_true_b <- stats::rbeta(B, a_sp, b_sp)

      n_d <- stats::rbinom(B, N, prev_b)
      n_nd <- N - n_d

      x_se <- stats::rbinom(B, n_d, se_true_b)
      x_sp <- stats::rbinom(B, n_nd, sp_true_b)

      post_a_se <- a_se + x_se
      post_b_se <- b_se + n_d - x_se
      width_se <- stats::qbeta(ci_upper_q, post_a_se, post_b_se) -
        stats::qbeta(ci_lower_q, post_a_se, post_b_se)

      post_a_sp <- a_sp + x_sp
      post_b_sp <- b_sp + n_nd - x_sp
      width_sp <- stats::qbeta(ci_upper_q, post_a_sp, post_b_sp) -
        stats::qbeta(ci_lower_q, post_a_sp, post_b_sp)

      success <- (width_se <= delta_se) & (width_sp <= delta_sp)
      # A degenerate replication (no diseased, or no non-diseased subjects)
      # always counts as a FAILURE and is never dropped from the denominator.
      degenerate <- (n_d == 0L) | (n_nd == 0L)
      success[degenerate] <- FALSE

      joint_assurance <- mean(success)
      assurance_mcse <- sqrt(joint_assurance * (1 - joint_assurance) / B)

      joint_assurance_last <- joint_assurance
      assurance_mcse_last <- assurance_mcse

      # Anti-noise rule: accept N only if the lower bound of the one-sided
      # 95% CI of the Monte Carlo estimate itself still reaches the target.
      if ((joint_assurance - z_mcse * assurance_mcse) >= target_assurance) {
        N_total <- as.integer(N)
        joint_assurance_achieved <- joint_assurance
        assurance_mcse_achieved <- assurance_mcse
        break
      }
    }

    if (is.na(N_total)) {
      warning("No N in N_range achieved the target JOINT assurance (at the ",
              "95% Monte Carlo confidence lower bound). Consider expanding ",
              "N_range or increasing B.")
      N_total <- as.integer(max(N_range_used))
      joint_assurance_achieved <- joint_assurance_last
      assurance_mcse_achieved <- assurance_mcse_last
    }
  } else {
    # --- method == "exact": closed-form Beta-Binomial enumeration, no
    # simulation, no seed/B dependence, no sampling error to guard against
    # (hence no anti-noise rule). See @details for the derivation.
    N_max_exact <- max(N_range_used)

    # O(N^2) qbeta calls to build the cache below; warn rather than silently
    # grinding for minutes on an oversized N_range.
    exact_n_max_practical <- 3000L
    if (N_max_exact > exact_n_max_practical) {
      warning(
        "method = \"exact\" caches per-arm probabilities up to N = ",
        N_max_exact, ", which is O(N^2) qbeta evaluations and may be very ",
        "slow. Consider method = \"monte_carlo\" for N_range this large."
      )
    }

    # Cached once per arm, up to N_max_exact, and reused across every
    # candidate N below -- P_se(k) / P_sp(m) depend only on arm size.
    P_se_exact <- .bam_exact_width_prob(
      N_max_exact, a_se, b_se, delta_se, ci_lower_q, ci_upper_q
    )
    P_sp_exact <- .bam_exact_width_prob(
      N_max_exact, a_sp, b_sp, delta_sp, ci_lower_q, ci_upper_q
    )

    # A degenerate replication (n_d = 0 or n_nd = 0) is always a FAILURE,
    # never a success -- the same convention already enforced under method =
    # "monte_carlo" (see @details). In the joint sum inside
    # .bam_exact_joint_assurance(), n_d = 0 occurs only at k = 0 and n_nd = 0
    # only at k = N, and both terms read the arm-size-0 entry of the
    # relevant cache (P_se_exact[1] for k = 0, P_sp_exact[1] for k = N),
    # regardless of which candidate N is being evaluated. Zeroing that entry
    # once, here, therefore removes exactly the degenerate contribution from
    # every candidate N's sum, with no change needed to
    # .bam_exact_joint_assurance() itself. Without this, P_arm(0) is simply
    # the credible-interval width evaluated at the PRIOR (no data), which is
    # not 0 in general and can even be 1 for tight informative priors --
    # wrongly crediting an arm that received no subjects at all.
    P_se_exact[1] <- 0
    P_sp_exact[1] <- 0

    joint_assurance_last <- NA_real_

    for (N in N_range_used) {
      joint_assurance <- .bam_exact_joint_assurance(
        N, prior_prev[1], prior_prev[2], P_se_exact, P_sp_exact
      )
      joint_assurance_last <- joint_assurance

      if (joint_assurance >= target_assurance) {
        N_total <- as.integer(N)
        joint_assurance_achieved <- joint_assurance
        # No Monte Carlo sampling error under the exact calculation -- see
        # @return for why 0 (not NA) is used here.
        assurance_mcse_achieved <- 0
        break
      }
    }

    if (is.na(N_total)) {
      warning("No N in N_range achieved the target JOINT assurance under ",
              "the exact calculation. Consider expanding N_range.")
      N_total <- as.integer(max(N_range_used))
      joint_assurance_achieved <- joint_assurance_last
      assurance_mcse_achieved <- 0
    }
  }

  # --- Buderer comparison ---
  buderer_n_se <- buderer_n(a_se / (a_se + b_se), delta_se / 2)

  structure(
    list(
      method = "Bayesian Assurance Method (BAM) for DTA Sample Size",
      n_diseased = n_se,
      n_non_diseased = n_sp,
      n_total = N_total,
      N_total = N_total,
      joint_assurance = joint_assurance_achieved,
      assurance_mcse = assurance_mcse_achieved,
      assurance_method = method,
      assurance_se = assurance_se_achieved,
      assurance_sp = assurance_sp_achieved,
      N_total_median = ceiling(N_total_median),
      N_total_P75 = ceiling(N_total_P75),
      N_total_P90 = ceiling(N_total_P90),
      buderer_n_se = buderer_n_se,
      prior_se = prior_se,
      prior_sp = prior_sp,
      prior_prev = prior_prev,
      delta_se = delta_se,
      delta_sp = delta_sp,
      target_assurance = target_assurance,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
