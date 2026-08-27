
#' Achievable Ceiling of the Net-Benefit Criterion as N -> Infinity
#'
#' Computes the largest joint assurance that \code{check_nb}'s CI-based
#' net-benefit criterion (see \code{\link{ss_unified}}) can ever reach,
#' for any total sample size \code{N}. As \code{N -> Inf} the half-widths
#' of both CI-based inequalities in the criterion shrink to 0, so the
#' criterion collapses to a statement about the marginal cell
#' probabilities alone: for every threshold \code{pt} in \code{pt_range}
#' (threshold odds \eqn{w = pt/(1-pt)}),
#' \deqn{p_1 - w\,p_2 > 0 \quad\text{and}\quad w\,q_2 - q_1 > 0,}
#' with
#' \deqn{p_1 = \mathrm{prev}\cdot Se\cdot Se_{ref} +
#'   (1-\mathrm{prev})(1-Sp)(1-Sp_{ref}),}
#' \deqn{p_2 = \mathrm{prev}\cdot Se\cdot(1-Se_{ref}) +
#'   (1-\mathrm{prev})(1-Sp)\,Sp_{ref},}
#' \deqn{q_1 = \mathrm{prev}\cdot(1-Se)\cdot Se_{ref} +
#'   (1-\mathrm{prev})\,Sp\,(1-Sp_{ref}),}
#' \deqn{q_2 = \mathrm{prev}\cdot(1-Se)\cdot(1-Se_{ref}) +
#'   (1-\mathrm{prev})\,Sp\,Sp_{ref}.}
#' \code{prev}, \code{Se} and \code{Sp} are drawn from their priors, so
#' this ceiling is a probability \strong{over the priors}, not over any
#' finite study: it is the fraction of prior draws for which a perfectly
#' precise (N = Inf) study would still conclude the test is useful at
#' every threshold in \code{pt_range}. No finite N can ever clear a
#' \code{target_assurance} above this ceiling.
#'
#' \strong{Why a wide \code{pt_range} tends to lower the ceiling.}
#' \eqn{p_1 - w\,p_2 > 0} is equivalent to \eqn{w < p_1/p_2}: raising the
#' top of \code{pt_range} raises \code{w} and makes this side of the
#' criterion (beating "treat none") harder to satisfy. \eqn{w\,q_2 - q_1 >
#' 0} is equivalent to \eqn{w > q_1/q_2}: lowering the bottom of
#' \code{pt_range} lowers \code{w} and makes the other side (beating
#' "treat all") harder to satisfy. A wide range therefore squeezes the
#' criterion from both ends at once, since \strong{all} thresholds in
#' \code{pt_range} must pass simultaneously.
#'
#' The Monte Carlo draws use their own generator call (see below) and do
#' not touch the caller's or \code{\link{ss_unified}}'s RNG stream --
#' both the RNG kind and \code{.Random.seed} are saved and restored on
#' exit, including on error (\code{\link{save_rng_state}}). Fully
#' vectorized -- no per-N or per-observation simulation is needed for an
#' N -> Inf limit -- so \code{B_ceiling} can be large without materially
#' adding to the cost of the \code{check_nb} search it guards.
#'
#' \strong{The generator kind is fixed, not inherited.} Versions <= 0.6.0
#' saved and restored \code{.Random.seed} but called \code{set.seed(seed)}
#' with no \code{kind}, so the draws above -- and hence the reported
#' ceiling -- silently depended on whatever generator the \strong{caller}
#' had active, exactly the defect 0.6.0 fixed in \code{\link{ss_unified}}
#' itself but left uncorrected in this helper. Measured with identical
#' arguments and seed, the reported ceiling ranged from 0.28253
#' (Marsaglia-Multicarry) to 0.28539 (Knuth-TAOCP-2002); with a
#' \code{pt_range} whose true ceiling sits just above 0.80, some
#' generators had the ceiling check correctly skip the (unreachable)
#' search while others let it run. \code{set.seed()} below now fixes the
#' kind explicitly to \code{"Mersenne-Twister"} (R's own default), so the
#' ceiling no longer depends on the caller's RNG configuration.
#'
#' \strong{Independent streams per \code{pt_range}, not a shared sample.}
#' A caller that evaluates this ceiling at several different
#' \code{pt_range} values under the same \code{seed} -- e.g. one row per
#' threshold range in a table -- would, under versions <= 0.6.0, draw the
#' \strong{exact same} \code{(prev, Se, Sp)} triplet for every one of
#' them: \code{pt_range} only entered the calculation \emph{after} the
#' draws, in the loop over thresholds below, so nothing about it affected
#' \code{set.seed()}. Whichever direction that one shared sample's noise
#' happened to point, every reported ceiling inherited it, in the same
#' direction and to a similar relative magnitude, since they were all
#' evaluated on literally the same draws -- a table of "independent"
#' ceiling estimates whose errors were, in fact, perfectly correlated
#' with each other and did not average out across rows. \code{seed} is
#' therefore now combined with a deterministic fingerprint of
#' \code{pt_range} (order-independent, since the joint criterion below is
#' itself symmetric in \code{pt_range}'s order) before seeding, so that
#' two calls sharing a \code{seed} but differing in \code{pt_range} draw
#' statistically independent samples, while a given \code{(seed,
#' pt_range)} pair remains exactly reproducible.
#'
#' @param prior_se,prior_sp,prior_prev Beta priors, as in
#'   \code{\link{ss_unified}}.
#' @param Se_ref,Sp_ref Reference standard accuracy, as in
#'   \code{\link{ss_unified}}.
#' @param pt_range Net-benefit thresholds, as in \code{\link{ss_unified}}.
#' @param seed Seed for the internal Monte Carlo draws. Combined with a
#'   fingerprint of \code{pt_range} before seeding; see \code{Details}.
#' @param B_ceiling Number of prior draws. Default 200000000 (2e8). The
#'   ceiling is a sample proportion, so its Monte Carlo standard error is
#'   \eqn{\sqrt{p(1-p)/B_{ceiling}}} -- a real bound on the estimator's
#'   sampling variability, not an after-the-fact observation -- which is
#'   below \eqn{2.6 \times 10^{-5}} at this default for any \eqn{p} in the
#'   0.7-0.95 range this function typically returns (worst case at the
#'   low end, \eqn{p = 0.7}: \eqn{\sqrt{0.7 \times 0.3 / 2\times 10^8}
#'   \approx 2.6\times10^{-5}}). What that bound has to overcome, for the
#'   package's own default \code{pt_range = c(0.15, 0.40)}, is a
#'   genuinely close call: the true ceiling there (independently derived
#'   by tensor-product Gauss-Legendre quadrature, not by this function)
#'   is approximately 0.879606, only \eqn{1.06 \times 10^{-4}} above
#'   0.8795 -- the boundary at which the value published to three
#'   decimals (\code{round(x, 3)}) flips from 0.879 to 0.880. At the
#'   previous default (\code{B_ceiling = 30000000}, whose standard error
#'   is close to \eqn{6\times10^{-5}} there), that boundary sat under two
#'   standard errors away, putting a non-trivial fraction of seeds
#'   (empirically, close to 3\%) on the wrong side of the rounded value.
#'   \code{B_ceiling = 200000000} pushes the standard error down far
#'   enough that the boundary sits more than four standard errors from
#'   the true value, making the published third decimal robust to the
#'   seed. \code{B_ceiling = 200000} (the default in versions <= 0.6.0)
#'   additionally reused one shared sample across every \code{pt_range} a
#'   caller evaluated (see \code{Details}), so its errors did not average
#'   out across a table of ceilings the way independent-sample Monte
#'   Carlo error normally does. An independent cross-check against the
#'   Gauss-Legendre reference across six \code{pt_range} scenarios
#'   (ceilings from 0.695 to 0.969) found the largest observed absolute
#'   error at the current default to be well inside this bound -- an
#'   empirical confirmation, not itself the guarantee; the guarantee is
#'   the standard-error bound above. Fully vectorized draws keep the
#'   runtime of this default at roughly two minutes per call on ordinary
#'   hardware (a few seconds at the more modest \code{B_ceiling} values
#'   sufficient for exploratory use), well under the cost of the
#'   \code{check_nb} grid search it guards at realistic \code{B} and
#'   \code{N_range}.
#' @return Numeric scalar in \eqn{[0, 1]}: the achievable ceiling.
#' @keywords internal
#' @noRd
nb_assurance_ceiling <- function(prior_se, prior_sp, prior_prev,
                                  Se_ref, Sp_ref, pt_range,
                                  seed, B_ceiling = 200000000L) {
  old_rng_state <- save_rng_state()
  on.exit(restore_rng_state(old_rng_state), add = TRUE)

  # --- seed salted by pt_range: see @details ------------------------------
  # A simple polynomial (Horner) hash of pt_range's sorted, microsecond-
  # rounded values, folded into `seed` by addition. Sorting first makes the
  # fingerprint -- and hence the stream -- depend only on the SET of
  # thresholds, not their order, matching the order-independence of the
  # joint criterion itself (see the loop below). The modulus 2147483647
  # (2^31 - 1, itself prime) keeps every intermediate value within double
  # precision's exact-integer range for the B_ceiling and pt_range sizes
  # this function is ever called with, and keeps the final seed a valid
  # (positive) argument to set.seed().
  pt_int <- sort(round(pt_range * 1e6))
  fp <- 0
  for (v in pt_int) fp <- (fp * 1000003 + v) %% 2147483647
  seed_eff <- as.integer((as.numeric(seed) + fp) %% 2147483647)
  # All three RNG kinds are named explicitly, not just the uniform
  # generator: see @details, "The generator kind is fixed, not
  # inherited," and the package-wide defect 0.6.2 fixes (NEWS.md) --
  # naming `kind` alone leaves `normal.kind`/`sample.kind` inherited from
  # the caller. This function only ever draws via stats::rbeta(), which
  # does not consume the normal or discrete-sampling streams, so the
  # ceiling itself does not depend on those two; both are still named
  # here, for the same completeness reason the guarantee is made general
  # elsewhere rather than function-specific.
  set.seed(seed_eff, kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")

  prev_b <- stats::rbeta(B_ceiling, prior_prev[1], prior_prev[2])
  Se_b   <- stats::rbeta(B_ceiling, prior_se[1],   prior_se[2])
  Sp_b   <- stats::rbeta(B_ceiling, prior_sp[1],   prior_sp[2])

  p1 <- prev_b * Se_b * Se_ref +
    (1 - prev_b) * (1 - Sp_b) * (1 - Sp_ref)
  p2 <- prev_b * Se_b * (1 - Se_ref) +
    (1 - prev_b) * (1 - Sp_b) * Sp_ref
  q1 <- prev_b * (1 - Se_b) * Se_ref +
    (1 - prev_b) * Sp_b * (1 - Sp_ref)
  q2 <- prev_b * (1 - Se_b) * (1 - Se_ref) +
    (1 - prev_b) * Sp_b * Sp_ref

  ok <- rep(TRUE, B_ceiling)
  for (pt in pt_range) {
    w <- pt / (1 - pt)
    ok <- ok & (p1 - w * p2 > 0) & (w * q2 - q1 > 0)
  }
  mean(ok)
}

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
#'   \strong{The CI-based net-benefit criterion has a hard ceiling that no
#'   N can cross.} As \eqn{N \to \infty} the half-widths of both
#'   inequalities above vanish, so the criterion collapses to a statement
#'   about the drawn \eqn{(\mathrm{prev}, Se, Sp)} alone (see
#'   \code{nb_assurance_ceiling}, internal). When \code{check_nb = TRUE},
#'   this ceiling is computed \strong{before} the grid search (from the
#'   priors, \code{Se_ref}, \code{Sp_ref} and \code{pt_range} only, by a
#'   fast vectorized Monte Carlo independent of \code{B} and \code{N_range})
#'   and returned as \code{nb_ceiling}. If \code{target_assurance} exceeds
#'   this ceiling, the criterion is unreachable at \strong{any} N, the grid
#'   search is skipped entirely (it cannot succeed, so nothing is gained by
#'   running it), and the function warns accordingly -- a message that
#'   \strong{replaces}, rather than joins, the generic "expand N_range"
#'   warning below, since expanding \code{N_range} cannot help this failure
#'   mode. A wide \code{pt_range} tends to push the ceiling down because
#'   raising its top makes the "beats treat-none" side harder while
#'   lowering its bottom makes the "beats treat-all" side harder, so both
#'   ends of a wide range squeeze the (necessarily joint, over all
#'   thresholds) criterion at once; see \code{nb_assurance_ceiling} for the
#'   full derivation.
#'
#'   \strong{The N selected by the search is itself a Monte Carlo estimate,
#'   and \code{decision} controls how that estimate is used.} At ordinary
#'   \code{B} the per-N Monte Carlo error is not small: at \code{B = 1200}
#'   the MCSE of an assurance near 0.80 is about 0.0115. Two distinct
#'   sources of noise are at play: the \emph{per-N} sampling error just
#'   described, and, across N, the \emph{level} of the whole assurance
#'   curve, which drifts up or down as one block under a fresh seed because
#'   the same priors and thresholds are shared by every N in the search.
#'   Empirically (a nested Se/Sp/AUC/net-benefit scenario checked against a
#'   B = 1e7 reference N* = 2216 across 20-25 seeds) the curve-level drift
#'   dominates the per-N noise, and it is what makes taking the
#'   \strong{first} N to cross a noisy curve a biased, high-variance
#'   stopping rule: across seeds it landed with sd 96.8, bias +100.5, RMSE
#'   138.1, and its nominal one-sided 95\% guarantee actually held in only
#'   84\% of seeds (16\% selected an N whose true assurance fell short of
#'   \code{target_assurance}).
#'
#'   \code{decision = "isotonic"} (default) instead evaluates the
#'   \strong{entire} \code{N_range} (regardless of \code{full_grid}),
#'   fits a monotone non-decreasing curve to \code{assurance} against
#'   \code{N} via \code{\link[stats]{isoreg}} (pooling adjacent violators;
#'   monotonicity itself is not in question -- assurance cannot decrease
#'   in N -- only its noisy estimate can), and inverts that fitted curve at
#'   \code{target_assurance} by linear interpolation between the two
#'   bracketing grid points, rounding up to the next integer N. Pooling
#'   the whole curve, rather than reading off a single noisy crossing,
#'   is what removes most of the curve-level bias described above: in the
#'   same benchmark this rule landed with sd 41.1, bias +10.4, RMSE 41.3.
#'   \code{joint_assurance} for this rule is the fitted (not raw) assurance
#'   at \code{N_effective}. \code{assurance_mcse} is not defined for a
#'   pooled, possibly off-grid estimate and is reported as \code{NA};
#'   \code{assurance_lower} is \strong{not} \code{NA} here -- it is the
#'   margin-adjusted (block-pooled Wilson) curve actually inverted to pick
#'   \code{N_effective}, described below and returned as a real number
#'   whenever \code{N_effective} is found. See \code{Return}.
#'
#'   Achieving this requires that the \code{B} replications at different N
#'   be genuinely independent draws, not the same draws replayed with a
#'   shift. Versions <= 0.5.0 called \code{set.seed(seed)} once per N,
#'   which is only \emph{partial} common random numbers: \code{rbinom(1,
#'   N, prev_b)} consumes a number of uniforms that depends on N, so the
#'   streams for different N agree for the first replication and diverge
#'   unpredictably from the second onward -- worse than either genuine CRN
#'   (which would help isolate a change between two specific N) or genuine
#'   independence (which the isotonic smoother needs, since it relies on
#'   neighbouring N's noise being uncorrelated). This version instead seeds
#'   the L'Ecuyer-CMRG generator once from \code{seed} and advances to a
#'   fresh stream for each N via \code{\link[parallel]{nextRNGStream}} --
#'   each N's B replications are statistically independent of every other
#'   N's, while the whole search remains exactly reproducible from the
#'   single \code{seed} argument, and the caller's own RNG state and kind
#'   are restored on exit exactly as before. This change alone moved the
#'   stopping-rule benchmark above from sd 96.8 (matching the behaviour of
#'   genuine common random numbers) toward the sd 53.9 expected of fully
#'   independent streams; common random numbers were tried deliberately and
#'   found to make the isotonic rule \emph{worse} (sd 91.7 vs 47.0), because
#'   they make the dominant curve-level drift perfectly correlated across
#'   N instead of averaging it out.
#'
#'   \code{decision = "lower_bound"} accepts the first N at which the lower
#'   limit of a one-sided 95\% Monte Carlo confidence bound on the
#'   assurance clears the target,
#'   \deqn{\mathrm{WilsonLower}_{0.95}(\hat{a}, B) \geq \text{target\_assurance},}
#'   where \eqn{\hat{a}} is \code{joint_assurance} at that N and
#'   \eqn{\mathrm{WilsonLower}} is the one-sided Wilson score lower bound
#'   (see below). This guards against per-N sampling error specifically --
#'   it selects a larger (or equal) N than \code{"point"} for the same
#'   draws -- but, being still a first-crossing rule, remains exposed to
#'   the curve-level drift above; it is the benchmark against which
#'   \code{"isotonic"}'s improvement was measured. Kept for a cheap,
#'   single-N-at-a-time search (it can still stop early when
#'   \code{full_grid = FALSE}) and for compatibility with versions
#'   0.4.1-0.5.0.
#'
#'   \code{decision = "point"} reproduces the behaviour of versions <= 0.4.0:
#'   an N is accepted as soon as \eqn{\hat{a} \geq} \code{target_assurance},
#'   with no correction for the fact that \eqn{\hat{a}} is itself noisy.
#'
#'   \strong{The lower confidence bound is Wilson, not Wald, so it cannot
#'   report a vacuous \code{1.000} at small B.} The textbook normal-
#'   approximation (Wald) bound, \eqn{\hat{a} - z_{0.95}\sqrt{\hat{a}(1 -
#'   \hat{a})/B}}, has a standard error that is \strong{exactly} 0 whenever
#'   \eqn{\hat{a} \in \{0, 1\}}, regardless of B: at \code{B = 1}, a single
#'   replicate that happens to pass every active target gives \eqn{\hat{a}
#'   = 1} with a Wald standard error of \eqn{\sqrt{1 \cdot 0 / 1} = 0}, so
#'   the Wald bound itself is \code{1.000} -- both \code{joint_assurance}
#'   and \code{assurance_lower} report total certainty from one replicate,
#'   under \code{decision \%in\% c("point", "lower_bound")} \strong{and}
#'   under \code{"isotonic"}'s block-pooled margin (see below), and every
#'   decision rule accepts the first N tried. This is an artifact of
#'   approximating the sampling distribution of \eqn{\hat{a}} by a normal
#'   centred at \eqn{\hat{a}} itself, which collapses to a point mass right
#'   where \eqn{\hat{a}} sits at a boundary; it says nothing about how
#'   confident one should actually be. \code{assurance_lower} is instead
#'   computed with the Wilson score bound (\code{wilson_lower()}, internal;
#'   see its documentation for the formula and references), which inverts
#'   the normal approximation to the \emph{score} rather than to
#'   \eqn{\hat{a}} and therefore stays strictly inside \eqn{(0, 1)} for any
#'   finite B -- e.g. at \code{B = 1} and \eqn{\hat{a} = 1} it gives
#'   approximately \code{0.27}, not \code{1.000}. It agrees with the Wald
#'   bound to \eqn{O(z^2/B)}, which is already on the order of
#'   \eqn{10^{-4}} or smaller at the \code{B = 20000} used for this
#'   package's own published results, so \strong{no result computed at a
#'   realistic B is affected by this change}; only the small-B behaviour
#'   is. \code{assurance_mcse} itself is left as the ordinary Wald standard
#'   error, \eqn{\sqrt{\hat{a}(1 - \hat{a})/B}} -- a purely descriptive
#'   quantity, documented as such in \code{Return} -- so it is still
#'   exactly 0 at \eqn{\hat{a} \in \{0, 1\}}; it is \code{assurance_lower}
#'   that no longer collapses.
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
#'   over \code{N_range} for \code{decision \%in\% c("lower_bound", "point")}
#'   as soon as the \strong{smallest} N in \code{N_range} that reaches
#'   \code{target_assurance} has been found (see below for why "smallest",
#'   not "first-in-traversal-order"). If \code{TRUE}, the search does not
#'   stop early: every N in \code{N_range} is evaluated so that
#'   \code{grid_results} holds the complete assurance curve. \code{optimal_N}
#'   and \code{joint_assurance} are unaffected by this switch for
#'   \code{decision \%in\% c("lower_bound", "point")} -- they always refer to
#'   the \strong{smallest} N that reached \code{target_assurance}, never a
#'   larger one, \strong{regardless of the order \code{N_range} was given
#'   in}. Versions <= 0.6.1 instead accepted whichever N reached
#'   \code{target_assurance} \emph{first in \code{N_range}'s own order},
#'   which is only the smallest N when \code{N_range} happens to already be
#'   sorted ascending: e.g. the identical set of candidate N reached
#'   \code{N_effective = 850} when \code{N_range} was given ascending but
#'   \code{1200} (the largest candidate, +41\%) when given descending, and
#'   yet another value when shuffled, with no warning of any kind, because
#'   descending/shuffled orders are evaluated largest-or-arbitrary-N-first
#'   and \code{decision \%in\% c("lower_bound", "point")} stopped at
#'   whichever N it reached first. \code{N_unique_sorted} (the distinct
#'   candidate N, ascending) is now always the traversal order for these two
#'   rules, independent of \code{N_range}'s own order, so the result -- and,
#'   under \code{full_grid = FALSE}, the search cost -- no longer depends on
#'   how \code{N_range} was arranged. \code{decision = "isotonic"} always
#'   evaluates the complete \code{N_range}, in the order given (one row of
#'   \code{grid_results} per element, including any repeats), regardless of
#'   \code{full_grid}, since the isotonic fit needs the whole curve; it
#'   already sorts internally before fitting and was already invariant to
#'   \code{N_range}'s order (see \code{Details}), so it is unaffected by
#'   this change and is more expensive than the default
#'   \code{full_grid = FALSE} cost of the other two rules.
#' @param decision How to decide the required N from the assurance curve:
#'   \code{"isotonic"} (default), \code{"lower_bound"} or \code{"point"}.
#'   See \code{Details}.
#' @param nb_B_ceiling Number of prior draws used by the internal
#'   \code{check_nb} ceiling calculation (\code{nb_assurance_ceiling()};
#'   see \code{Details}, "The CI-based net-benefit criterion has a hard
#'   ceiling"). Ignored when \code{check_nb = FALSE}. Default 200000000,
#'   which keeps the ceiling's Monte Carlo standard error below about
#'   \eqn{2.6 \times 10^{-5}} for ceilings in the typical 0.7-0.95 range (see
#'   \code{nb_assurance_ceiling}'s own documentation for how this was
#'   chosen and verified). Exposed mainly so that a caller who only needs an
#'   approximate ceiling -- e.g. while exploring \code{pt_range} choices,
#'   or in a test that checks the ceiling's presence rather than its exact
#'   value -- can pass a much smaller value for speed; the calculation is
#'   fully vectorized and independent of \code{B} and \code{N_range}, so
#'   its cost does not otherwise scale with the rest of the search.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{joint_assurance}{Achieved joint assurance at \code{N_effective}:
#'       the probability that \strong{all} active targets -- Se, Sp, AUC
#'       (when \code{delta_auc > 0}) and net benefit (when
#'       \code{check_nb = TRUE}) -- are reached \strong{simultaneously} in
#'       the same replication, not the probability that each is reached
#'       marginally. The denominator is \code{B}: degenerate replications
#'       count as failures. Under \code{decision = "isotonic"} this is the
#'       fitted (smoothed), not raw, assurance at \code{N_effective} --
#'       the actual point estimate, \strong{not} margin-adjusted downward
#'       (see \code{assurance_lower} and \code{Details} for the margin that
#'       \code{N_effective} was actually selected against).}
#'     \item{assurance_mcse}{Monte Carlo standard error of
#'       \code{joint_assurance} at \code{N_effective}, i.e.
#'       \code{sqrt(joint_assurance * (1 - joint_assurance) / B)}.
#'       \code{NA} under \code{decision = "isotonic"} (not defined for a
#'       pooled, possibly off-grid estimate).}
#'     \item{assurance_lower}{Lower limit of the one-sided 95\% Wilson score
#'       confidence bound on \code{joint_assurance} (see \code{Details},
#'       "The lower confidence bound is Wilson, not Wald"), for
#'       \code{decision \%in\% c("point", "lower_bound")}; this is the
#'       quantity compared against \code{target_assurance} when
#'       \code{decision = "lower_bound"} (see \code{Details}). Under
#'       \code{decision = "isotonic"} it is instead the same Wilson bound
#'       applied to the fitted (smoothed) curve at its block-pooled
#'       effective sample size, evaluated at \code{N_effective} (see
#'       \code{Details}): the quantity \code{"isotonic"} actually inverts
#'       to pick \code{N_effective}, and so is guaranteed \code{>=
#'       target_assurance} whenever \code{N_effective} was found (never
#'       \code{NA} in that case). Unlike the normal-approximation (Wald)
#'       bound used in versions <= 0.6.0, this cannot report a vacuous
#'       \code{1.000} when a tiny B happens to pass every replicate.}
#'     \item{nb_ceiling}{The achievable ceiling of the \code{check_nb}
#'       criterion as \eqn{N \to \infty} (see \code{Details}); \code{NA}
#'       when \code{check_nb = FALSE}.}
#'     \item{comparison}{Data frame comparing methods (columns \code{method}
#'       and \code{N}). The imperfect-reference correction contributes two
#'       rows, not one: \code{estimand = "apparent"} (the exact closed-form
#'       Buderer-style calculation, no inflation factor needed) and
#'       \code{estimand = "corrected"} (the delta-method misclassification
#'       correction for the true Se/Sp), since a single row cannot
#'       represent both target quantities at once. See
#'       \code{\link{ss_imperfect_ref}}.}
#'     \item{N_buderer}{Total N from the classical Buderer formula, i.e.
#'       the "Buderer (classical)" row of \code{comparison}.}
#'     \item{N_imperfect}{Total N for the imperfect-reference correction's
#'       \strong{apparent} estimand (\code{\link{ss_imperfect_ref}}'s
#'       default), i.e. the "apparent estimand" row of \code{comparison}.
#'       Kept under this name for backward compatibility; the (typically
#'       much larger) misclassification-corrected N is the "corrected
#'       estimand" row of \code{comparison}.}
#'     \item{seed}{The \code{seed} argument used for the search.}
#'     \item{target_assurance}{The \code{target_assurance} argument used for
#'       the search, echoed back for downstream use (e.g. by
#'       \code{\link{plot_assurance_curve}}).}
#'     \item{decision}{The \code{decision} argument used for the search.}
#'     \item{grid_results}{Data frame with columns \code{N} and
#'       \code{assurance}, one row per N evaluated while searching
#'       \code{N_range}. Under \code{decision = "isotonic"} this is one row
#'       per \strong{element} of \code{N_range} (repeats included, in the
#'       order given) whenever the search runs at all, since the isotonic
#'       fit always evaluates the complete range. Under
#'       \code{decision \%in\% c("lower_bound", "point")} it is one row per
#'       \strong{distinct} N actually evaluated -- ascending, and with
#'       \code{full_grid = TRUE} covering every unique value in
#'       \code{N_range} -- since a repeated N is looked up rather than
#'       re-evaluated (see \code{@param full_grid}); with the default
#'       \code{full_grid = FALSE} the search instead stops as soon as the
#'       smallest N reaching \code{target_assurance} is found, so the grid
#'       is truncated there.}
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
                       full_grid = FALSE,
                       decision = c("isotonic", "lower_bound", "point"),
                       nb_B_ceiling = 200000000L) {
  # --- preserve the caller's RNG state -------------------------------
  # This function switches the generator to L'Ecuyer-CMRG (see below) to
  # draw an independent stream per N via parallel::nextRNGStream(), so
  # restoring just the .Random.seed VALUE is not enough: set.seed() called
  # later by anyone else with no explicit `kind` argument reuses whatever
  # kind is currently active, not whatever kind .Random.seed happens to
  # encode, so a caller left on L'Ecuyer-CMRG after this function returns
  # gets L'Ecuyer-CMRG draws even when they set their own seed. Both the
  # RNG kind (RNGkind(), i.e. the uniform, normal and discrete-sampling
  # generators) and the .Random.seed vector are therefore saved here and
  # restored on exit -- including on error, via on.exit().
  old_rng_kind <- RNGkind()
  old_seed_present <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_present) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    # Setting the kind (even back to what it already was) can itself
    # create .Random.seed with an arbitrary value when none existed
    # before, so the exact seed vector -- or its absence -- is restored
    # AFTER the kind, overwriting whatever RNGkind() just did.
    suppressWarnings(RNGkind(
      kind = old_rng_kind[1],
      normal.kind = old_rng_kind[2],
      sample.kind = old_rng_kind[3]
    ))
    if (old_seed_present) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  # Validate inputs
  decision <- match.arg(decision)
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
  stopifnot(nb_B_ceiling >= 1)

  # --- validate N_range explicitly, with a diagnostic message -----------
  # Versions <= 0.6.1 did not validate N_range at all: a negative, NA,
  # infinite or non-integer element was not rejected here, so it instead
  # reached rbinom()/rmultinom() deep inside the search loop and died with
  # a generic, uninformative "missing value where TRUE/FALSE needed" --
  # a crash with no indication of which argument, or which element of it,
  # was at fault. Checked in an order that reports the most fundamental
  # problem first (each check assumes the ones before it already passed).
  if (!is.numeric(N_range) || length(N_range) < 1) {
    stop("N_range must be a numeric vector of length >= 1 (got ",
         if (length(N_range) == 0) "a zero-length vector" else class(N_range)[1],
         ").", call. = FALSE)
  }
  if (anyNA(N_range)) {
    stop("N_range contains NA (at position(s) ",
         paste(which(is.na(N_range)), collapse = ", "),
         "). Every element must be a finite positive integer.", call. = FALSE)
  }
  if (any(!is.finite(N_range))) {
    stop("N_range contains a non-finite value (Inf or -Inf) at position(s) ",
         paste(which(!is.finite(N_range)), collapse = ", "),
         ". Every element must be a finite positive integer.", call. = FALSE)
  }
  if (any(N_range < 1)) {
    stop("N_range contains a value below 1 (minimum given: ", min(N_range),
         "). Every element must be a positive integer total sample size.",
         call. = FALSE)
  }
  if (any(N_range != round(N_range))) {
    stop("N_range contains a non-integer value (e.g. ",
         N_range[which(N_range != round(N_range))[1]],
         "). Every element must be a whole number (total sample size).",
         call. = FALSE)
  }

  warn_small_B(B)

  z <- stats::qnorm(0.975)
  z_decision <- stats::qnorm(0.95)  # one-sided 95% lower bound on assurance
  target_se_width <- 2 * delta_se
  target_sp_width <- 2 * delta_sp
  target_auc_width <- 2 * delta_auc

  # Prior means for comparison table
  E_prev <- prior_prev[1] / sum(prior_prev)
  E_se <- prior_se[1] / sum(prior_se)
  E_sp <- prior_sp[1] / sum(prior_sp)  # Sp arm for the Buderer row

  # --- net-benefit ceiling: computed BEFORE the search, see @details ---
  nb_ceiling_val <- NA_real_
  if (check_nb) {
    nb_ceiling_val <- nb_assurance_ceiling(
      prior_se = prior_se, prior_sp = prior_sp, prior_prev = prior_prev,
      Se_ref = Se_ref, Sp_ref = Sp_ref, pt_range = pt_range, seed = seed,
      B_ceiling = nb_B_ceiling
    )
  }
  # The criterion is unreachable at ANY N: no point running the (much more
  # expensive) grid search, since it is guaranteed to fail everywhere.
  skip_search <- check_nb && !is.na(nb_ceiling_val) &&
    (target_assurance > nb_ceiling_val)

  optimal_N <- NA_integer_
  joint_assurance_achieved <- NA_real_
  assurance_mcse_achieved <- NA_real_
  assurance_lower_achieved <- NA_real_
  joint_assurance <- 0  # defensive init (in case N_range is degenerate)
  assurance_mcse <- 0
  assurance_lower <- -Inf

  # one (N, assurance) pair per N evaluated; with the default
  # full_grid = FALSE the loop below still breaks at the optimum (unless
  # decision = "isotonic", which always evaluates the whole range), so this
  # is only as long as the grid actually searched (see @param full_grid).
  grid_N <- vector("integer", length(N_range))
  grid_assurance <- vector("numeric", length(N_range))
  grid_idx <- 0L

  if (!skip_search) {
    # --- independent (not common-random-number) RNG stream per N -------
    # Seed the L'Ecuyer-CMRG generator once from `seed`, then advance to a
    # fresh, statistically independent stream for each DISTINCT N via
    # parallel::nextRNGStream(). See @details for why the previous
    # per-N set.seed(seed) pattern was only PARTIAL common random numbers,
    # and why decision = "isotonic" needs genuine independence rather than
    # either that accident or deliberate CRN.
    #
    # Streams are generated in one fixed order -- increasing N -- and then
    # looked up by the VALUE of N, not by N's position in N_range. Keying
    # by position would make the stream a given N receives depend on where
    # that N happens to sit in the N_range the caller passed in, so
    # reordering N_range (descending, shuffled, seq(..., by = -10), a
    # hand-typed c() not given in order, ...) would silently hand different
    # N values to different streams and change every downstream result even
    # though the same SET of N values and the same seed were used. Keying
    # by value instead guarantees that a given N always gets the same
    # stream for a given seed, regardless of where it appears in N_range.
    # This alone makes the SIMULATED DATA at each N invariant to N_range's
    # order, and, together with the explicit sort before the isotonic fit
    # below, is what makes decision = "isotonic" fully invariant to
    # N_range's order. It is NOT, by itself, enough to make
    # decision \%in\% c("point", "lower_bound") invariant: those two rules
    # additionally traverse N in ASCENDING order below (via
    # N_unique_sorted, not N_range's own order) so that "the first N
    # encountered that reaches target_assurance" means the SMALLEST such N
    # regardless of how N_range was arranged -- see @param full_grid.
    # Versions <= 0.6.1 traversed N_range in the order given for every
    # decision rule, so value-keyed RNG streams alone left the SIMULATED
    # DATA order-invariant while the SELECTED N for "point"/"lower_bound"
    # was not: a descending or shuffled N_range could -- and did, with no
    # warning -- select a larger N than the same set given ascending. A
    # duplicated N value in N_range is looked up the same way both times,
    # so it deterministically receives the very same stream and therefore
    # the very same simulated replications -- correct, since it is the same
    # design being asked about twice, not two independent looks at it.
    # normal.kind/sample.kind are named explicitly too, even though this
    # search only ever draws via stats::rbeta()/rbinom()/rmultinom() and so
    # does not itself consume the normal or discrete-sampling streams --
    # see the same completeness note in nb_assurance_ceiling() above and
    # NEWS.md's 0.6.2 entry.
    set.seed(seed, kind = "L'Ecuyer-CMRG",
             normal.kind = "Inversion", sample.kind = "Rejection")
    stream_state <- .Random.seed
    N_unique_sorted <- sort(unique(N_range))
    rng_stream_by_N <- vector("list", length(N_unique_sorted))
    for (j in seq_along(N_unique_sorted)) {
      rng_stream_by_N[[j]] <- stream_state
      stream_state <- parallel::nextRNGStream(stream_state)
    }

    # --- traversal order: see the note above and @param full_grid --------
    # decision = "isotonic" needs the whole curve and already sorts before
    # fitting, so it keeps evaluating N_range in the order (and with any
    # repeats) the caller gave, exactly as before -- one grid_results row
    # per element of N_range. decision \%in\% c("point", "lower_bound")
    # instead evaluates each DISTINCT N exactly once, in ascending order,
    # so that "the first N encountered that reaches target_assurance" is
    # always the smallest such N, regardless of N_range's own order or of
    # any repeated values in it (defect 1).
    N_seq_to_eval <- if (decision == "isotonic") N_range else N_unique_sorted

    for (i in seq_along(N_seq_to_eval)) {
      N <- N_seq_to_eval[i]
      stream_idx <- match(N, N_unique_sorted)
      assign(".Random.seed", rng_stream_by_N[[stream_idx]], envir = .GlobalEnv)
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
      # Monte Carlo error of that proportion (reported as-is; see @return),
      # and the resulting one-sided 95% lower confidence bound -- a WILSON
      # score bound, not the normal-approximation (Wald) bound
      # joint_assurance - z * assurance_mcse it replaced. The two agree
      # closely at the B this package's own results are computed with (see
      # wilson_lower()'s documentation), but the Wald bound is degenerate at
      # the extremes: at B = 1 a single passing replicate gives
      # joint_assurance = 1 with assurance_mcse = sqrt(1 * 0 / 1) = 0
      # EXACTLY, so the Wald bound reports assurance_lower = 1.000 --
      # manufactured certainty from one replicate, regardless of B. The
      # Wilson bound never does this: it stays strictly below 1 for any
      # finite B (see @details for why the bound, not the point estimate,
      # is the acceptance criterion under decision = "lower_bound").
      assurance_mcse <- sqrt(joint_assurance * (1 - joint_assurance) / B)
      assurance_lower <- wilson_lower(joint_assurance, B, z_decision)

      grid_idx <- grid_idx + 1L
      grid_N[grid_idx] <- N
      grid_assurance[grid_idx] <- joint_assurance

      if (decision != "isotonic") {
        accept <- if (decision == "point") {
          joint_assurance >= target_assurance
        } else {
          assurance_lower >= target_assurance
        }

        # optimal_N is always the SMALLEST N to reach the target: N_seq_to_eval
        # is N_unique_sorted (ascending) for these two rules, so the first
        # acceptance encountered here IS the smallest one, regardless of
        # N_range's own order (defect 1). Once set, a later (necessarily
        # larger) N -- only reachable with full_grid = TRUE -- must not
        # overwrite it.
        if (is.na(optimal_N) && accept) {
          optimal_N <- as.integer(N)
          joint_assurance_achieved <- joint_assurance
          assurance_mcse_achieved <- assurance_mcse
          assurance_lower_achieved <- assurance_lower
          if (!isTRUE(full_grid)) break
        }
      }
      # decision == "isotonic": never break -- the whole N_range is always
      # needed to fit the curve (see below, after the loop).
    }
  }

  grid_results <- data.frame(
    N = grid_N[seq_len(grid_idx)],
    assurance = grid_assurance[seq_len(grid_idx)],
    stringsAsFactors = FALSE
  )

  # --- decision = "isotonic": invert the pooled curve at target_assurance ---
  if (!skip_search && decision == "isotonic" && nrow(grid_results) > 0) {
    # stats::isoreg() returns $x in the ORIGINAL order of its input and
    # $yf in the SORTED order -- an internal detail of its implementation,
    # not something documented as part of its contract. Pairing them
    # positionally (as reading $x, $yf off the fit invites) silently
    # mismatches each N against a fitted value that belongs to a DIFFERENT
    # N whenever the input order is not already increasing. grid_results$N
    # follows N_range's own order (see the search loop above), so this
    # fired, with no warning of any kind, for any N_range that was not
    # itself already sorted ascending (descending, shuffled, a hand-typed
    # c() not given in order, ...). Sorting explicitly before the call
    # removes the ambiguity entirely: once the input is already increasing,
    # isoreg()'s "original order" and "sorted order" coincide by
    # construction, so $x and $yf come back aligned to that same order no
    # matter what order N_range itself was in.
    ord <- order(grid_results$N)
    N_sorted <- grid_results$N[ord]
    assurance_sorted <- grid_results$assurance[ord]

    # --- deduplicate BEFORE fitting: a repeated N carries ZERO new
    # information, not an extra independent observation (defect 3) --------
    # A duplicated N in N_range always replays the very same L'Ecuyer-CMRG
    # stream (see the search loop above), so its assurance is byte-identical
    # to the first copy's -- literally the same simulated draws, not a fresh
    # independent look at that N. Versions <= 0.6.1 fit isoreg() on the
    # UNDEDUPLICATED sorted grid, so a repeated N contributed extra,
    # perfectly-correlated "replicates" to whatever pooled block it fell
    # into below, inflating that block's `block_len` -- and hence its
    # effective sample size B * block_len in the Wilson margin -- with no
    # new evidence behind it. Measured on a 16-point base grid with each
    # point replicated 1x/2x/3x/5x/10x (otherwise identical B, seed,
    # priors): N_effective fell from 857 to 835 as the replication count
    # rose, purely from this fabricated precision, with the reported
    # assurance_lower correspondingly, and wrongly, INCREASING. Since
    # duplicate N's are guaranteed identical, deduplicating first (keeping
    # each distinct N once) changes nothing about which values are fit --
    # only how many times an unchanged value is counted -- so a grid with
    # no duplicates is completely unaffected by this step.
    dedup_keep <- !duplicated(N_sorted)
    N_fit <- N_sorted[dedup_keep]
    assurance_fit <- assurance_sorted[dedup_keep]

    iso <- stats::isoreg(N_fit, assurance_fit)
    Nseq <- iso$x        # == N_fit; verified by the defensive check below
    fitted <- iso$yf     # pooled, non-decreasing fitted assurance, aligned to Nseq

    # --- defensive checks: this alignment must never fail silently --------
    # A future change to isoreg()'s internals -- or a slip in the sorting/
    # dedup above -- that broke the Nseq/fitted correspondence would
    # otherwise reproduce exactly the silent wrong-answer failure mode this
    # fix removes. A warning would not do here: the original bug's whole
    # failure mode WAS a warning-free silent mismatch, so any violation of
    # these invariants stops the function outright instead.
    if (length(Nseq) != length(fitted) || length(Nseq) != length(N_fit)) {
      stop(
        "Internal error in ss_unified(): the isotonic fit's N and fitted-",
        "assurance vectors do not have the expected length (Nseq = ",
        length(Nseq), ", fitted = ", length(fitted), ", N_fit = ",
        length(N_fit), "). Aborting rather than risk silently ",
        "pairing an N with the wrong fitted assurance; please report this ",
        "as a bug."
      )
    }
    # isoreg() always returns $x as double, even when grid_results$N (and
    # therefore N_fit) happens to be stored as integer -- e.g. whenever
    # N_range itself was built with `:` (as in 1600:1600) rather than
    # seq(). That storage-type difference is expected and harmless, so the
    # comparison below is by VALUE (as.double() on both sides), not by
    # identical(), which would otherwise treat integer 1600L and double
    # 1600 as different and trip this check on a distinction that has
    # nothing to do with real misalignment.
    if (is.unsorted(Nseq) || !identical(as.double(Nseq), as.double(N_fit))) {
      stop(
        "Internal error in ss_unified(): stats::isoreg() did not return ",
        "$x in the expected non-decreasing, pre-sorted order, so Nseq and ",
        "fitted can no longer be relied on to be aligned. Aborting rather ",
        "than risk silently pairing an N with the wrong fitted assurance; ",
        "please report this as a bug."
      )
    }

    # --- conservative margin on top of the smoothed curve (see @details) --
    # Smoothing alone removes most of the curve-level bias but leaves NO
    # safety margin: inverting `fitted` at target_assurance exactly picks
    # the N where the SMOOTHED POINT ESTIMATE crosses the target, and a
    # point estimate is above its true value about half the time. A
    # one-sided 95% Monte Carlo margin (the same z used by
    # decision = "lower_bound") is therefore subtracted from `fitted`
    # BEFORE inversion, so the rule accepts the first N whose smoothed
    # estimate clears target_assurance by more than its own sampling noise.
    #
    # The margin's standard error is NOT each N's own raw (pre-smoothing)
    # MCSE, sqrt(raw * (1 - raw) / B): isoreg() pools each run of
    # consecutive violating N into one flat block sharing a single fitted
    # value (identifiable from `fitted` itself -- rle() finds the flat
    # runs), and under that block's own working assumption -- that the
    # true assurance is constant across it, which is exactly what pooling
    # asserts -- the fitted value for a block of `block_len` grid points is
    # the mean of `block_len` independent per-N estimates (independent
    # because each N draws its own L'Ecuyer-CMRG stream), so its true
    # standard error is sqrt(fitted * (1 - fitted) / (B * block_len)): the
    # single-N formula divided by sqrt(block_len). Using the single-N MCSE
    # everywhere ignores this pooling and is needlessly conservative on
    # wide blocks; it was tried and pushed N_effective past even the
    # pre-isotonic first-crossing rule's own bias on the step-5 benchmark
    # in @details (a regression, not an improvement), which is why the
    # block-pooled margin is used instead. It remains conservative rather
    # than exact wherever the constant-assurance working assumption is
    # only approximate within a block.
    #
    # The margin itself is a WILSON score bound on `fitted`, with effective
    # sample size B * block_len, not the Wald bound
    # fitted - z * mcse_N implied by the paragraph above. The two agree
    # closely at this effective sample size for any fitted value bounded
    # away from 0 or 1 -- which is the entire useful range of a target
    # assurance -- but the Wald bound is degenerate at fitted = 1 (its
    # standard error is exactly 0 there regardless of B or block_len), so
    # for a tiny B (B = 1 is the extreme case) a block whose single
    # replicate happened to pass everything reports fitted_lower = 1.000:
    # manufactured certainty from the smallest possible amount of evidence,
    # not a real margin. The Wilson bound stays strictly below 1 for any
    # finite effective sample size instead. mcse_N is kept, unchanged, as
    # the descriptive block-pooled Monte Carlo standard error referenced
    # above; it no longer directly forms fitted_lower by subtraction.
    block_len <- rep(rle(fitted)$lengths, rle(fitted)$lengths)

    # --- cap the pooling multiplier at B: a minimal B cannot manufacture a
    # strong margin (defect 5) -------------------------------------------
    # rle(fitted) cannot tell a block that isoreg() genuinely POOLED to
    # enforce monotonicity apart from a run of adjacent, already-distinct N
    # whose RAW assurance simply happened to coincide -- both look like an
    # identical run of `fitted` values. At a realistic B the raw assurance
    # is close to continuous (multiples of 1/B), so two distinct N's raw
    # values landing on EXACTLY the same estimate by chance is negligible,
    # and long flat runs reliably indicate genuine pooling; block_len is
    # therefore trustworthy as-is there (defect 3's dedup above handles the
    # OTHER way a block gets inflated, literal repeated N). At a tiny B,
    # though, raw assurance can only take B + 1 possible values -- at
    # B = 1, just {0, 1} -- so long flat runs of DISTINCT, non-duplicated N
    # arise routinely by chance alone, not because the underlying curve is
    # genuinely flat there; treating such a run as B * block_len
    # independent trials then manufactures a deceptively strong bound from
    # almost no evidence. Reproduced directly: a generous scenario at
    # B = 1 across a wide N_range let a single pass-or-fail replicate per N
    # (no duplicates anywhere in N_range) pool into pooled blocks dozens of
    # points wide, reporting assurance_lower close to 0.81 -- built on
    # exactly one bit of information per contributing N. Capping the
    # pooling multiplier at B (never trusting more than B * B = B^2
    # "trials" worth of pooled evidence when B itself is small, i.e.
    # `min(block_len, B)`) makes this impossible: at B = 1 the effective
    # sample size behind ANY block is at most 1 * 1 = 1, so
    # wilson_lower(fitted, 1, z_decision) tops out under 0.27 -- nowhere
    # near a "strong" bound -- regardless of how wide the (spurious) flat
    # run happens to be. The cap only bites when block_len exceeds B, which
    # requires more grid points pooled together than B itself; every
    # existing scenario in this package's own tests and manuscript uses a
    # grid with far fewer points than its B (or a realistic B in the
    # thousands to tens of thousands), so the cap is a complete no-op there
    # -- min(block_len, B) == block_len whenever block_len <= B, which
    # holds throughout ordinary use.
    block_len_capped <- pmin(block_len, B)
    mcse_N <- sqrt(fitted * (1 - fitted) / (B * block_len_capped))
    fitted_lower <- wilson_lower(fitted, B * block_len_capped, z_decision)

    idx_reach <- which(fitted_lower >= target_assurance)
    if (length(idx_reach) > 0) {
      k <- idx_reach[1]
      if (k == 1) {
        N_star <- Nseq[1]
      } else {
        g_lo <- fitted_lower[k - 1]
        g_hi <- fitted_lower[k]
        N_star <- if (g_hi > g_lo) {
          Nseq[k - 1] +
            (target_assurance - g_lo) / (g_hi - g_lo) * (Nseq[k] - Nseq[k - 1])
        } else {
          Nseq[k]
        }
      }
      optimal_N <- as.integer(ceiling(N_star))

      # A duplicated N in N_range produces a duplicated value in Nseq (with
      # an identical fitted value at each copy -- see the RNG-stream note
      # above: a repeated N always replays the same stream, hence the same
      # result). stats::approx() would still work on such ties, but it warns
      # ("collapsing to unique 'x' values") every time it sees one, which is
      # exactly the noise a merely-repeated grid point should NOT produce.
      # Collapsing the ties ourselves first keeps that expected, harmless
      # case silent; it changes nothing about the interpolation itself,
      # since the collapsed points carried identical y-values to begin with.
      keep <- !duplicated(Nseq)
      N_interp <- Nseq[keep]
      fitted_interp <- fitted[keep]
      fitted_lower_interp <- fitted_lower[keep]

      # With fewer than two usable grid points (e.g. a single-N N_range, or
      # an N_range that is a single N repeated) there is nothing to
      # interpolate -- stats::approx() errors below two x values -- and the
      # isotonic fit at one point is just that point's own raw value, so
      # both curves are read off directly.
      if (length(N_interp) >= 2) {
        joint_assurance_achieved <- stats::approx(
          x = N_interp, y = fitted_interp, xout = optimal_N, rule = 2
        )$y
        assurance_lower_achieved <- stats::approx(
          x = N_interp, y = fitted_lower_interp, xout = optimal_N, rule = 2
        )$y
      } else {
        joint_assurance_achieved <- fitted_interp[1]
        assurance_lower_achieved <- fitted_lower_interp[1]
      }
      # joint_assurance is the UNADJUSTED fitted (smoothed) curve at
      # N_effective -- the actual point estimate of the assurance, not the
      # conservative value the margin was checked against. assurance_lower
      # is that margin-adjusted curve itself: because margin >= 0 pointwise
      # and both curves are interpolated with the same weights, it is
      # algebraically guaranteed that
      #   joint_assurance_achieved >= assurance_lower_achieved
      #     >= target_assurance
      # whenever optimal_N is found here, mirroring the role
      # assurance_lower already plays under decision = "lower_bound" (the
      # quantity actually compared against target_assurance). assurance_mcse
      # remains undefined for a pooled, possibly off-grid estimate.
      assurance_mcse_achieved <- NA_real_
    }
    # else: the margin-adjusted curve never reaches target_assurance within
    # N_range -- optimal_N stays NA and falls through to the generic
    # "did not converge" branch below, exactly as for the other rules.
  }

  found <- !is.na(optimal_N)
  if (skip_search) {
    # Replaces (not joins) the generic "expand N_range" message: expanding
    # N_range cannot help when the criterion is unreachable at any N.
    warning(
      "check_nb's net-benefit criterion cannot reach target_assurance = ",
      target_assurance, " at ANY N with these priors, Se_ref = ", Se_ref,
      ", Sp_ref = ", Sp_ref, " and pt_range = [",
      paste(range(pt_range), collapse = ", "), "]: as N -> Inf the ",
      "achievable joint assurance converges to a ceiling of ",
      sprintf("%.5f", nb_ceiling_val), ", which is below target_assurance. ",
      "Expanding N_range cannot help. Lower target_assurance below the ",
      "ceiling, narrow pt_range, or use more informative priors for Se, ",
      "Sp and prevalence; see nb_ceiling in the returned object and ",
      "?ss_unified for the structural reason a wide pt_range lowers it.",
      call. = FALSE
    )
    optimal_N <- max(N_range)
    joint_assurance_achieved <- NA_real_
    assurance_mcse_achieved <- NA_real_
    assurance_lower_achieved <- NA_real_
  } else if (!found) {
    warning("No N in N_range achieved target assurance. ",
            "Consider expanding N_range.")
    optimal_N <- max(N_range)
    # Report the RAW assurance actually observed AT the N being returned
    # (max(N_range)), not whichever N happened to be evaluated LAST by the
    # search loop above (defect 2). Those coincide only when N_range's
    # traversal order happens to end at its own maximum -- true for an
    # ascending N_range under decision = "isotonic" (which always
    # traverses N_range in the order given), but false for a descending or
    # shuffled N_range, and false for decision \%in\% c("point",
    # "lower_bound") whenever N_unique_sorted's ascending traversal was cut
    # short by full_grid = FALSE. Reproduced directly: the identical
    # non-converging search (same seed, B, priors, SET of candidate N)
    # reported joint_assurance = 0.37850 given ascending, 0.06400 given
    # descending and 0.30000 given shuffled -- three different numbers
    # attached to the SAME returned N_effective, because each run reported
    # whatever N its loop happened to finish on rather than the assurance
    # at max(N_range) itself. grid_results always contains a row for
    # max(N_range): the loop above only stops early (decision \%in\%
    # c("point", "lower_bound") with full_grid = FALSE) upon acceptance,
    # and acceptance is exactly what did NOT happen here, so every branch
    # that can reach this point evaluated the complete search space.
    # Duplicate rows at max(N_range) (a repeated value in N_range) are
    # guaranteed identical (same value-keyed RNG stream), so the first is
    # used without ambiguity.
    max_N_assurance <- grid_results$assurance[grid_results$N == optimal_N][1]
    joint_assurance_achieved <- max_N_assurance
    assurance_mcse_achieved <- sqrt(max_N_assurance * (1 - max_N_assurance) / B)
    assurance_lower_achieved <- wilson_lower(max_N_assurance, B, z_decision)
  } else if (decision != "isotonic" && assurance_mcse_achieved >
               (joint_assurance_achieved - target_assurance) / 2) {
    # The selected N clears target_assurance (by whichever rule `decision`
    # applies), but only barely relative to its own Monte Carlo error: a
    # different seed, or a slightly larger B, could plausibly move the
    # verdict. This is a distinct condition from warn_small_B(), which
    # flags B on its own regardless of how close the search landed to the
    # target; this one flags a search result that landed close to the
    # target with an implausibly small margin of error behind it. Not
    # meaningful under decision = "isotonic", whose assurance_mcse is NA.
    warning(
      "The joint assurance at the selected N (N_effective = ", optimal_N,
      ") is close to its own Monte Carlo noise: assurance = ",
      round(joint_assurance_achieved, 4), ", MCSE = ",
      round(assurance_mcse_achieved, 4), ", target = ", target_assurance,
      ". Consider increasing B for a more reliable N_effective.",
      call. = FALSE
    )
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

  # ss_imperfect_ref() always computes and returns BOTH estimands (see its
  # Details): the apparent estimand needs no inflation factor at all -- it
  # is the exact closed-form calculation with Se_apparent/Sp_apparent in
  # place of Se/Sp -- and the corrected estimand recovers the true Se/Sp via
  # the delta-method variance of the exactly-identified misclassification
  # correction. They are reported as two separate rows below, rather than
  # collapsed into a single "inflation" row, because one scalar cannot
  # represent both target quantities at once.
  comparison <- data.frame(
    method = c("Buderer (classical)",
               "Imperfect ref, apparent estimand (exact)",
               "Imperfect ref, corrected estimand (misclassification)",
               "Unified (this method)"),
    N = c(buderer_N_total,
          imperfect_res$N_apparent_loss,
          imperfect_res$N_corrected_loss,
          N_enrolled),
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
      assurance_mcse = assurance_mcse_achieved,
      assurance_lower = assurance_lower_achieved,
      nb_ceiling = nb_ceiling_val,
      comparison = comparison,
      N_buderer = buderer_N_total,
      # apparent estimand explicitly (not the estimand-selected generic
      # alias): this is what "N_imperfect" has always meant here, and stays
      # so regardless of what estimand a future caller might select for
      # ss_imperfect_ref() elsewhere.
      N_imperfect = imperfect_res$N_apparent_loss,
      seed = seed,
      target_assurance = target_assurance,
      decision = decision,
      grid_results = grid_results,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
