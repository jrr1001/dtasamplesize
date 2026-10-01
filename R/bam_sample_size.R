# --- Internal helpers for the exact (non-Monte-Carlo) joint search --------
# Not exported. Used only by bam_sample_size() when method = "exact"; see
# that function's @details for the closed-form derivation these implement.

# For a single arm governed by a Beta(a, b) prior and a full-width credible
# interval target of delta, return a vector P of length (n_hi - n_lo + 1)
# where P[n - n_lo + 1] = P(width <= delta | arm size = n), for n running
# from n_lo to n_hi, marginalizing over the arm's own success count
# x ~ BetaBinomial(n, a, b). The posterior credible interval width after x
# successes in n trials is a deterministic function of (x, n) alone
# (Beta(a + x, b + n - x)), so this sum is a finite, non-random enumeration
# -- no simulation involved.
#
# Each entry depends only on (n, a, b, delta), never on N, on any OTHER
# arm size, or on which candidate N is being evaluated -- entries for
# different n are mutually independent. This is what lets
# bam_sample_size() extend its exact-mode cache block by block (see
# .bam_exact_extend_cache() below): a call with (n_lo, n_hi) = (501, 1000)
# returns exactly the same values for n = 501..1000 as a call with
# (n_lo, n_hi) = (0, 1000) would, so re-querying only the NEW arm sizes
# when a block grows is safe and produces byte-identical numbers to
# recomputing the whole range from scratch.
.bam_exact_width_prob <- function(n_hi, a, b, delta, ci_lower_q, ci_upper_q,
                                   n_lo = 0L) {
  n_lo <- as.integer(n_lo)
  n_hi <- as.integer(n_hi)
  P <- numeric(n_hi - n_lo + 1L)
  lbeta_ab <- lbeta(a, b)
  for (n in n_lo:n_hi) {
    x <- 0:n
    post_a <- a + x
    post_b <- b + n - x
    width <- stats::qbeta(ci_upper_q, post_a, post_b) -
      stats::qbeta(ci_lower_q, post_a, post_b)
    log_pmf <- lchoose(n, x) + lbeta(post_a, post_b) - lbeta_ab
    P[n - n_lo + 1L] <- sum(exp(log_pmf)[width <= delta])
  }
  P
}

# Extend an exact-mode per-arm cache (as built by .bam_exact_width_prob())
# from covering arm sizes 0:old_hi to covering 0:new_hi, computing ONLY the
# new entries (old_hi+1):new_hi and concatenating them onto what is already
# cached. `old_cache` may be NULL (nothing cached yet) or length 0, in which
# case this simply builds 0:new_hi from scratch. Returns the extended cache
# (length new_hi + 1).
.bam_exact_extend_cache <- function(old_cache, old_hi, new_hi, a, b, delta,
                                     ci_lower_q, ci_upper_q) {
  if (is.null(old_cache) || old_hi < 0L) {
    return(.bam_exact_width_prob(new_hi, a, b, delta, ci_lower_q, ci_upper_q,
                                  n_lo = 0L))
  }
  if (new_hi <= old_hi) {
    return(old_cache[seq_len(new_hi + 1L)])
  }
  new_part <- .bam_exact_width_prob(new_hi, a, b, delta, ci_lower_q,
                                     ci_upper_q, n_lo = old_hi + 1L)
  c(old_cache, new_part)
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
#'   convention, and the same threshold, already adopted by
#'   \code{\link{ss_unified}} and \code{\link{joint_sample_size}}. An arm of
#'   size \strong{one} is \emph{not} treated as degenerate: the posterior
#'   Beta credible interval and (in \code{joint_sample_size()}) the Wilson
#'   interval are both well-defined at \eqn{n = 1}.
#'
#'   This corrects a defect present through version 0.4.0, where \code{n_se}
#'   and \code{n_sp} were found \strong{independently}, each against its own
#'   \strong{marginal} assurance, and the previous \code{n_total} (now
#'   \code{N_total_median}, see below) was only the median of
#'   \code{pmax(n_se / prev, n_sp / (1 - prev))} over draws of \code{prev}.
#'   That quantity never verified that both widths hold in the same study,
#'   and its joint assurance is systematically lower than
#'   \code{target_assurance} -- e.g., using the default \code{prior_se},
#'   \code{prior_sp}, \code{delta_se}, \code{delta_sp} and
#'   \code{target_assurance} together with \code{prior_prev = c(4, 16)}
#'   (the worked example from the accompanying article; mean prevalence
#'   0.20, \strong{not} the default \code{prior_prev}, which is
#'   \code{c(6, 14)}), the old \code{n_total} achieved a true joint
#'   assurance of about 0.72 against a target of 0.80.
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
#'   \code{P_se(k) * P_sp(N - k)} inside the sum. There is no sampling error
#'   to guard against in this exact calculation, so \code{B} plays no role
#'   in it at all (it is still used for the legacy per-arm searches and
#'   heuristics described below) and \code{seed} does not affect the result
#'   either -- this is a genuine property of the arithmetic (a finite sum
#'   with no random inputs), not an approximation. \code{method = "exact"}
#'   is deterministic, reproducible bit-for-bit, and free of the dependence
#'   on \code{B} that affects \code{method = "monte_carlo"} (e.g., using the
#'   default \code{prior_se}, \code{prior_sp}, \code{delta_se}, \code{delta_sp}
#'   and \code{target_assurance} together with \code{prior_prev = c(4, 16)}
#'   (the worked example from the accompanying article; mean prevalence
#'   0.20, \strong{not} the default \code{prior_prev}, which is
#'   \code{c(6, 14)}), the exact first \code{N} with joint assurance
#'   \eqn{\ge} 0.80 is 678, with assurance 0.800349, whereas
#'   \code{method = "monte_carlo"} with \code{B = 20000} returns 683 because
#'   of Monte Carlo noise near the crossing point -- see \strong{Anti-noise
#'   acceptance rule} below).
#'
#'   \strong{Evaluation is exact with respect to the discrete model above,
#'   up to floating-point and \code{qbeta} quantile error.} "Exact" here
#'   means the Beta-Binomial sum is evaluated by enumeration rather than by
#'   simulation -- it is not a claim that \code{qbeta()} itself is free of
#'   floating-point rounding; see \strong{N = 2 vs. N = 1} below for how
#'   this is handled at the smallest arm sizes, where it matters most.
#'
#'   \strong{Which \code{N} are searched: the integer scan (\code{N_range =
#'   NULL}) vs. a user-supplied \code{N_range}.} Versions \verb{<= 0.6.5} (a
#'   defect corrected here) built the automatic search grid from \code{n_se},
#'   \code{n_sp} and \code{N_total_P90} -- all three themselves Monte Carlo
#'   quantities depending on \code{B} and \code{seed} -- so the SAME exact,
#'   deterministic calculation could silently return different \code{N_total}
#'   values depending on arguments (\code{B}, \code{seed}) that \code{method
#'   = "exact"} otherwise has, and documents, NO dependence on. When
#'   \code{N_range = NULL} (the default), the search instead scans every
#'   integer \code{N} starting at \code{N = 2} (see \strong{N = 2 vs. N = 1}
#'   below), strictly ascending, built in doubling blocks capped at
#'   \code{N_max}: the first block covers \code{2:min(500, N_max)}; if no
#'   crossing is found there, the next block extends the scan to
#'   \code{min(1000, N_max)}, then \code{min(2000, N_max)}, doubling the
#'   block ceiling each time until either a crossing is found or the ceiling
#'   reaches \code{N_max}. This never depends on \code{n_se}, \code{n_sp},
#'   \code{N_total_P90}, \code{B} or \code{seed}, and is exhaustive over
#'   integers up to \code{N_max}: the \code{N} it returns is the smallest
#'   integer whose exact joint assurance reaches \code{target_assurance},
#'   full stop, not merely the smallest one a heuristic grid happened to
#'   include. Supplying \code{N_range} explicitly instead searches exactly
#'   \code{sort(unique(as.integer(N_range)))}, in ascending order; the result
#'   is then "the first candidate in \code{N_range} that reaches
#'   \code{target_assurance}" -- the smallest value actually searched, not
#'   necessarily the smallest INTEGER overall, since an explicit
#'   \code{N_range} may skip values. \code{search_type} in the return value
#'   records which of the two was used (\code{"integer_scan_auto"} or
#'   \code{"user_N_range"}), and \code{N_range_used} records every value
#'   actually evaluated, in the order evaluated.
#'
#'   \strong{No crossing found (both \code{method}s).} If the scan (or the
#'   supplied \code{N_range}) is exhausted with no \code{N} reaching
#'   \code{target_assurance}, versions \verb{<= 0.6.5} (a defect corrected
#'   here) silently returned \code{max(N_range)} as if it were a validated
#'   solution. This function now instead sets \code{N_total = NA_integer_},
#'   \code{n_total = NA_integer_}, \code{joint_assurance = NA_real_},
#'   \code{target_reached = FALSE}, and reports the diagnostics
#'   \code{max_assurance_evaluated} (the highest exact/Monte-Carlo-estimated
#'   assurance seen anywhere in the scan) and \code{N_at_max_assurance} (the
#'   \code{N} at which it was seen), together with a warning naming the
#'   ceiling reached and how to search further (\code{N_max} for the
#'   automatic integer scan, a wider \code{N_range} otherwise).
#'
#'   \strong{N = 2 vs. N = 1.} The integer scan starts at \code{N = 2}, not
#'   \code{N = 1}: a total of \code{N = 1} puts every subject on one arm and
#'   none on the other (\code{n_d = 1, n_nd = 0} or vice versa), so its exact
#'   joint assurance is \strong{always} exactly \code{0} regardless of the
#'   priors -- the same degenerate-arm override already described above
#'   forces the zero-size arm's contribution to \code{0} at both \code{k = 0}
#'   and \code{k = N = 1}. \code{N = 1} is therefore never a useful search
#'   target and is skipped by construction; it is still evaluated correctly
#'   (as \code{0}) if a caller explicitly includes it in \code{N_range}. An
#'   arm of size exactly \strong{one} (which first becomes possible at
#'   \code{N = 2}, when the OTHER arm gets the other one) is, by contrast,
#'   \strong{not} degenerate: its one-subject posterior is well-defined and
#'   is scored on its actual width, matching \code{method = "monte_carlo"}
#'   (see the \code{H-06} regression test).
#'
#'   \code{P_se} and \code{P_sp} depend only on arm size, not on \code{N},
#'   so they are cached once per block (see above) and reused across every
#'   candidate \code{N} in that block; growing a block only computes the
#'   NEW arm sizes (see \code{.bam_exact_extend_cache()}), not the whole
#'   cache from scratch. The cache itself costs \eqn{O(N_{\max}^2)}
#'   \code{qbeta} evaluations in the worst case (about 500,000 for an
#'   \code{N_max} around 700, a few seconds). If the final block ceiling
#'   is large enough that this becomes impractical, a warning suggests
#'   \code{method = "monte_carlo"} instead.
#'
#'   \code{method = "monte_carlo"} instead estimates the same probability by
#'   simulation, as described above, and is kept for continuity with
#'   versions \verb{<= 0.4.x} of this joint search (introduced mid-cycle in
#'   0.5.0) and as a fallback for \code{N} ranges too large for the exact
#'   cache. Under \code{N_range = NULL} it still uses the legacy heuristic
#'   grid (\code{max(20, n_se, n_sp):ceiling(N_total_P90)}) described below,
#'   since a Monte Carlo search already has its own sampling error to manage
#'   and gains nothing from the deterministic integer scan used by
#'   \code{method = "exact"}; \code{search_type} for this case is
#'   \code{"monte_carlo_auto"}.
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
#'   \code{N_range = NULL} (the default): under \code{method = "exact"} it
#'   is the deterministic integer scan described above (never a function of
#'   \code{n_se}, \code{n_sp}, \code{B} or \code{seed}); under \code{method =
#'   "monte_carlo"} it remains the legacy heuristic grid
#'   \code{max(20, n_se, n_sp):ceiling(N_total_P90)}, i.e. from the smallest
#'   total that could possibly supply both per-arm requirements up to the
#'   90th percentile of the legacy prevalence-uncertainty heuristic (see
#'   \code{N_total_P90} below), which in practice comfortably brackets the
#'   true joint requirement for that search mode. No existing argument name
#'   was removed; \code{N_range} is purely additive.
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
#'   Default \code{NULL}: under \code{method = "exact"} this triggers the
#'   deterministic integer scan described in \code{Details} (independent of
#'   \code{n_se}, \code{n_sp}, \code{B} and \code{seed}); under \code{method
#'   = "monte_carlo"} it builds the legacy heuristic grid from \code{n_se},
#'   \code{n_sp}, and \code{N_total_P90} -- see \code{Details}. When supplied
#'   explicitly, searched as \code{sort(unique(as.integer(N_range)))}.
#' @param N_max Integer. Upper ceiling for the automatic integer scan used
#'   by \code{method = "exact"} when \code{N_range = NULL}. Default
#'   \code{3000L}, the same practical limit the exact cache already warned
#'   about. Ignored when \code{N_range} is supplied explicitly, or under
#'   \code{method = "monte_carlo"}. See \code{Details}.
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
#'       cohort design) reaches \code{target_assurance}, among the candidate
#'       N actually searched (see \code{N_range_used}). Under \code{method =
#'       "monte_carlo"} acceptance is at the 95\% Monte Carlo confidence
#'       lower bound (the anti-noise rule). \code{NA_integer_} when no
#'       candidate reached the target (\code{target_reached = FALSE}); see
#'       \code{Details}, "No crossing found". This is the headline result of
#'       this function when a solution exists; see \code{Details}.}
#'     \item{target_reached}{Logical. \code{TRUE} if some candidate \code{N}
#'       reached \code{target_assurance}; \code{FALSE} otherwise, in which
#'       case \code{N_total}, \code{n_total} and \code{joint_assurance} are
#'       all \code{NA} and \code{max_assurance_evaluated} /
#'       \code{N_at_max_assurance} carry the diagnostics instead. See
#'       \code{Details}, "No crossing found".}
#'     \item{N_range_used}{Integer vector, the candidate \code{N} values
#'       actually evaluated, in the order evaluated. Under the automatic
#'       integer scan (\code{method = "exact"}, \code{N_range = NULL}) this
#'       is \code{2:N_total} when a crossing was found, or the full scanned
#'       range up to the final block ceiling (at most \code{N_max}) when it
#'       was not.}
#'     \item{search_type}{One of \code{"integer_scan_auto"} (exact method,
#'       \code{N_range = NULL}), \code{"user_N_range"} (an explicit
#'       \code{N_range} was supplied, either method), or
#'       \code{"monte_carlo_auto"} (Monte Carlo method, \code{N_range =
#'       NULL}). See \code{Details}.}
#'     \item{N_max}{The \code{N_max} argument actually in effect for the
#'       automatic integer scan (\code{search_type == "integer_scan_auto"});
#'       \code{NA_integer_} otherwise (an explicit \code{N_range} or
#'       \code{method = "monte_carlo"} does not use \code{N_max}).}
#'     \item{joint_assurance}{The joint assurance achieved at \code{N_total}
#'       (point estimate, not the lower confidence bound used to accept it
#'       under \code{method = "monte_carlo"}; under \code{method = "exact"}
#'       this already \emph{is} the exact probability, with no further bound
#'       to distinguish it from). \code{NA_real_} when
#'       \code{target_reached = FALSE} -- see \code{max_assurance_evaluated}
#'       for the diagnostic value in that case.}
#'     \item{max_assurance_evaluated}{The highest assurance (exact
#'       probability, or Monte Carlo point estimate, according to
#'       \code{method}) observed anywhere among \code{N_range_used}. Always
#'       populated, including when \code{target_reached = TRUE} (where it
#'       equals \code{joint_assurance}, since assurance is evaluated in
#'       ascending \code{N} and the search stops at the first crossing).
#'       This is the field to read for "how close did the search get" when
#'       \code{target_reached = FALSE}.}
#'     \item{N_at_max_assurance}{The \code{N} (one element of
#'       \code{N_range_used}) at which \code{max_assurance_evaluated} was
#'       observed. Equals \code{N_total} when \code{target_reached = TRUE}.}
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
                            N_max = 3000L,
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
  stopifnot(length(N_max) == 1, N_max >= 2)

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

  # --- Joint search for total N (headline result, v0.5.0; M-01 contract,
  # v0.6.6) ------------------------------------------------------------
  # Unlike the two searches above, this one evaluates BOTH arms within the
  # SAME replication of a single cohort of size N, with n_d ~ Binomial(N,
  # prev) (the disease count is random, not fixed at floor(N * prev)). See
  # @details for the full generative model, the anti-noise acceptance rule,
  # the deterministic integer scan (method = "exact", N_range = NULL), and
  # the "no crossing found" contract shared by both methods.
  exact_n_max_practical <- 3000L

  N_total <- NA_integer_
  joint_assurance_achieved <- NA_real_
  assurance_mcse_achieved <- NA_real_
  target_reached <- FALSE
  max_assurance_evaluated <- NA_real_
  N_at_max_assurance <- NA_integer_
  N_max_out <- NA_integer_

  if (method == "monte_carlo") {
    # N_range = NULL keeps the legacy heuristic grid for this mode (a Monte
    # Carlo search already manages its own sampling error and gains nothing
    # from the deterministic integer scan used by method = "exact"; see
    # @details).
    if (is.null(N_range)) {
      N_lo <- max(20L, n_se, n_sp)
      N_hi <- max(N_lo + 1L, ceiling(N_total_P90))
      N_range_used <- N_lo:N_hi
      search_type <- "monte_carlo_auto"
    } else {
      N_range_used <- N_range
      search_type <- "user_N_range"
    }

    z_mcse <- stats::qnorm(0.95)

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

      if (is.na(max_assurance_evaluated) || joint_assurance > max_assurance_evaluated) {
        max_assurance_evaluated <- joint_assurance
        N_at_max_assurance <- as.integer(N)
      }

      # Anti-noise rule: accept N only if the lower bound of the one-sided
      # 95% CI of the Monte Carlo estimate itself still reaches the target.
      if ((joint_assurance - z_mcse * assurance_mcse) >= target_assurance) {
        N_total <- as.integer(N)
        joint_assurance_achieved <- joint_assurance
        assurance_mcse_achieved <- assurance_mcse
        target_reached <- TRUE
        break
      }
    }

    if (!target_reached) {
      # M-01 (v0.6.6): no longer silently returns max(N_range) as if it
      # were a validated solution (the defect this corrects; see
      # @details, "No crossing found"). N_total / n_total /
      # joint_assurance are NA_integer_/NA_real_; the diagnostics
      # max_assurance_evaluated and N_at_max_assurance report what the
      # search actually saw, and N_range covers the fix: expanding it is
      # the only way to search further under this mode.
      warning(
        "No N in N_range achieved the target JOINT assurance (at the 95% ",
        "Monte Carlo confidence lower bound): the highest assurance seen ",
        "was ", sprintf("%.6f", max_assurance_evaluated), " at N = ",
        N_at_max_assurance, ", against target_assurance = ",
        target_assurance, ". N_total, n_total and joint_assurance are NA ",
        "(target_reached = FALSE). Expand N_range or increase B to search ",
        "further; see max_assurance_evaluated / N_at_max_assurance for ",
        "diagnostics.",
        call. = FALSE
      )
    }
  } else {
    # --- method == "exact": closed-form Beta-Binomial enumeration, no
    # simulation, no sampling error to guard against (hence no anti-noise
    # rule), and (since v0.6.6) no dependence on seed/B for the headline
    # N_total / joint_assurance either way -- see @details.
    if (is.null(N_range)) {
      # --- deterministic integer scan, N = 2, 3, 4, ..., in doubling
      # blocks capped at N_max (see @details). Never depends on n_se,
      # n_sp, N_total_P90, B or seed.
      search_type <- "integer_scan_auto"
      N_max_out <- as.integer(N_max)

      block_cap <- min(500L, as.integer(N_max))
      start_N <- 2L
      cache_hi <- -1L
      P_se_exact <- NULL
      P_sp_exact <- NULL
      warned_slow <- FALSE

      repeat {
        if (!warned_slow && block_cap > exact_n_max_practical) {
          warning(
            "method = \"exact\" caches per-arm probabilities up to N = ",
            block_cap, ", which is O(N^2) qbeta evaluations and may be ",
            "very slow. Consider method = \"monte_carlo\" for an N_max ",
            "this large."
          )
          warned_slow <- TRUE
        }
        # Extend (not rebuild) the per-arm caches: only arm sizes
        # (cache_hi + 1):block_cap are newly computed; see
        # .bam_exact_extend_cache().
        P_se_exact <- .bam_exact_extend_cache(
          P_se_exact, cache_hi, block_cap, a_se, b_se, delta_se,
          ci_lower_q, ci_upper_q
        )
        P_sp_exact <- .bam_exact_extend_cache(
          P_sp_exact, cache_hi, block_cap, a_sp, b_sp, delta_sp,
          ci_lower_q, ci_upper_q
        )
        cache_hi <- block_cap
        # Degenerate arm-size-0 override (see the long comment below, kept
        # with the final non-auto branch); idempotent, so reapplying it on
        # every block extension is harmless.
        P_se_exact[1] <- 0
        P_sp_exact[1] <- 0

        for (N in start_N:block_cap) {
          joint_assurance <- .bam_exact_joint_assurance(
            N, prior_prev[1], prior_prev[2], P_se_exact, P_sp_exact
          )
          if (is.na(max_assurance_evaluated) || joint_assurance > max_assurance_evaluated) {
            max_assurance_evaluated <- joint_assurance
            N_at_max_assurance <- N
          }
          if (joint_assurance >= target_assurance) {
            N_total <- as.integer(N)
            joint_assurance_achieved <- joint_assurance
            assurance_mcse_achieved <- 0
            target_reached <- TRUE
            break
          }
        }

        if (target_reached || block_cap >= as.integer(N_max)) break
        start_N <- block_cap + 1L
        block_cap <- min(block_cap * 2L, as.integer(N_max))
      }

      N_range_used <- if (target_reached) 2L:N_total else 2L:block_cap

      if (!target_reached) {
        warning(
          "No integer N in 2:", N_max, " (N_max) achieved the target ",
          "JOINT assurance under the exact calculation: the highest ",
          "assurance seen was ", sprintf("%.6f", max_assurance_evaluated),
          " at N = ", N_at_max_assurance, ", against target_assurance = ",
          target_assurance, ". N_total, n_total and joint_assurance are NA ",
          "(target_reached = FALSE). Increase N_max to search further; ",
          "see max_assurance_evaluated / N_at_max_assurance for ",
          "diagnostics.",
          call. = FALSE
        )
      }
    } else {
      # --- user-supplied N_range: searched as sort(unique(as.integer(.))),
      # ascending -- see @details, "Which N are searched".
      search_type <- "user_N_range"
      N_range_used <- sort(unique(as.integer(N_range)))
      N_max_exact <- max(N_range_used)

      # O(N^2) qbeta calls to build the cache below; warn rather than
      # silently grinding for minutes on an oversized N_range.
      if (N_max_exact > exact_n_max_practical) {
        warning(
          "method = \"exact\" caches per-arm probabilities up to N = ",
          N_max_exact, ", which is O(N^2) qbeta evaluations and may be ",
          "very slow. Consider method = \"monte_carlo\" for N_range this ",
          "large."
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
      # never a success -- the same convention already enforced under
      # method = "monte_carlo" (see @details). In the joint sum inside
      # .bam_exact_joint_assurance(), n_d = 0 occurs only at k = 0 and
      # n_nd = 0 only at k = N, and both terms read the arm-size-0 entry of
      # the relevant cache (P_se_exact[1] for k = 0, P_sp_exact[1] for
      # k = N), regardless of which candidate N is being evaluated. Zeroing
      # that entry once, here, therefore removes exactly the degenerate
      # contribution from every candidate N's sum, with no change needed to
      # .bam_exact_joint_assurance() itself. Without this, P_arm(0) is
      # simply the credible-interval width evaluated at the PRIOR (no
      # data), which is not 0 in general and can even be 1 for tight
      # informative priors -- wrongly crediting an arm that received no
      # subjects at all.
      P_se_exact[1] <- 0
      P_sp_exact[1] <- 0

      for (N in N_range_used) {
        joint_assurance <- .bam_exact_joint_assurance(
          N, prior_prev[1], prior_prev[2], P_se_exact, P_sp_exact
        )
        if (is.na(max_assurance_evaluated) || joint_assurance > max_assurance_evaluated) {
          max_assurance_evaluated <- joint_assurance
          N_at_max_assurance <- N
        }

        if (joint_assurance >= target_assurance) {
          N_total <- as.integer(N)
          joint_assurance_achieved <- joint_assurance
          # No Monte Carlo sampling error under the exact calculation --
          # see @return for why 0 (not NA) is used here.
          assurance_mcse_achieved <- 0
          target_reached <- TRUE
          break
        }
      }

      if (!target_reached) {
        warning(
          "No N in N_range achieved the target JOINT assurance under the ",
          "exact calculation: the highest assurance seen was ",
          sprintf("%.6f", max_assurance_evaluated), " at N = ",
          N_at_max_assurance, ", against target_assurance = ",
          target_assurance, ". N_total, n_total and joint_assurance are NA ",
          "(target_reached = FALSE). Consider expanding N_range; see ",
          "max_assurance_evaluated / N_at_max_assurance for diagnostics.",
          call. = FALSE
        )
      }
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
      target_reached = target_reached,
      N_range_used = N_range_used,
      search_type = search_type,
      N_max = N_max_out,
      joint_assurance = joint_assurance_achieved,
      max_assurance_evaluated = max_assurance_evaluated,
      N_at_max_assurance = N_at_max_assurance,
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
