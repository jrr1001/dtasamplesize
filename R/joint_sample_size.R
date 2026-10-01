
#' Joint Sample Size for Sensitivity, Specificity, and AUC
#'
#' Finds the \strong{first candidate \code{N} in \code{N_range}} (searched
#' in the order given) at which sensitivity and specificity simultaneously
#' achieve their target precision with joint probability at least
#' \code{target_prob}, subject to the AUC also achieving its target
#' precision. Se and Sp are evaluated by Monte Carlo; AUC precision is
#' computed deterministically via the Hanley-McNeil variance approximation.
#' \strong{Read the \code{Details}: the AUC component is not probabilistic,
#' and the result is \strong{not} the smallest integer \code{N} unless
#' \code{N_range} itself contains every integer.}
#'
#' @details \strong{The estimand.} At each candidate total \code{N}, Se and
#'   Sp are point values (not random draws of a prior, unlike
#'   \code{\link{bam_sample_size}}): sensitivity, specificity and prevalence
#'   are the \strong{expected} (point) values supplied via \code{Se},
#'   \code{Sp} and \code{prev}. The Se and Sp confidence intervals are
#'   \strong{Wilson} intervals (frequentist), not the Beta-Binomial credible
#'   intervals used by \code{bam_sample_size()}; \code{delta_se} and
#'   \code{delta_sp} are their \strong{half}-widths (the target \strong{full}
#'   CI width checked internally is \code{2 * delta}), unlike
#'   \code{bam_sample_size()}'s \code{delta_se}/\code{delta_sp}, which are
#'   full widths. \code{joint_sample_size()} and \code{bam_sample_size()}
#'   therefore measure \strong{related but distinct estimands}, evaluated by
#'   different computational routes (Wilson/Monte Carlo here; Beta-Binomial
#'   credible interval, exact or Monte Carlo, there): \strong{do not treat
#'   one as an independent validation of the other}. The joint probability
#'   returned (\code{joint_prob_se_sp}) is over Se and Sp \strong{only}; see
#'   below for how the AUC criterion enters.
#'
#'   \strong{"First candidate in \code{N_range}", not "smallest \code{N}".}
#'   The search evaluates \code{N_range} in the order given (ascending, for
#'   the default grid) and stops at the first candidate whose joint
#'   probability reaches \code{target_prob}. This is the smallest value
#'   \strong{actually searched}, not necessarily the smallest \strong{integer}
#'   \code{N} with that property, unless \code{N_range} happens to contain
#'   every integer in its span (the default, \code{seq(100, 800, by = 10)},
#'   does not: it is a grid of step 10). For example, at the package's
#'   documented worked example (\code{Se = 0.85}, \code{Sp = 0.90},
#'   \code{prev = 0.20}, \code{delta_se = 0.07}, \code{delta_sp = 0.05},
#'   \code{B = 20000}, \code{seed = 2026}), the default grid returns
#'   \code{n_total = 580} with \code{joint_prob_se_sp = 0.8041}, but a finer,
#'   integer-step search over the same neighbourhood crosses the target
#'   \code{target_prob = 0.80} as early as \code{N = 578}-\code{579}
#'   (verified independently; see the package's audit trail). \code{580} is
#'   correctly described as "the first point of the step-10 grid reaching
#'   the target", \strong{not} as "the minimum sample size". \code{search_type
#'   = "grid_first_candidate"} and \code{N_range_used} in the return value
#'   make this explicit and queryable.
#'
#'   The parameters \code{delta_se}, \code{delta_sp}, and
#'   \code{delta_auc} are \strong{half-widths} of the confidence interval.
#'   The target full CI width checked internally is \code{2 * delta}.
#'
#'   \strong{What \code{joint_prob_se_sp} is, and what it is not.} The
#'   returned probability is the Monte Carlo probability that the Se
#'   \strong{and} Sp confidence intervals both meet their target width. It
#'   is a joint probability over \strong{two} of the three criteria, not
#'   three. The AUC criterion is applied as a \strong{deterministic gate}: at
#'   each candidate \code{N} the Hanley-McNeil variance is evaluated at the
#'   \emph{assumed} \code{AUC} and the candidate is discarded unless the
#'   resulting CI width meets \code{2 * delta_auc}. Because that width is
#'   computed at the assumed parameter rather than simulated, the AUC
#'   component carries no sampling variability: it is a
#'   \strong{median-style criterion}, satisfied in roughly 50\% of real
#'   studies (the ones whose realised AUC standard error lands below its
#'   expectation), \strong{not} an assurance of \code{target_prob}. The
#'   overall probability that all three intervals meet their targets in a
#'   real study is therefore \strong{lower} than the reported
#'   \code{joint_prob_se_sp}. Treat \code{joint_prob_se_sp} as the assurance
#'   for the Se/Sp pair, conditional on an AUC criterion that is met on
#'   average.
#'
#'   \strong{Geometric coherence of \code{AUC} with \code{Se} and \code{Sp}.}
#'   \code{Se}, \code{Sp} and \code{AUC} are not free parameters: any ROC
#'   curve that is concave and passes through the operating point
#'   \eqn{(1 - Sp,\ Se)} has area at least that of the polygon joining
#'   \eqn{(0,0)}, \eqn{(1 - Sp,\ Se)} and \eqn{(1,1)}, namely
#'   \deqn{AUC_{min} = \tfrac{1}{2}(1 - Sp)\,Se + \tfrac{1}{2} Sp\,(Se + 1).}
#'   Supplying an \code{AUC} below \eqn{AUC_{min}} describes a test that
#'   cannot exist, and the function now stops with an error rather than
#'   returning a sample size for it. At the defaults \code{Se = 0.85} and
#'   \code{Sp = 0.90} this bound is \strong{0.875}; the previous default of
#'   \code{AUC = 0.80} was therefore impossible. For reference, the binormal
#'   AUC implied by the same operating point,
#'   \eqn{\Phi[(\Phi^{-1}(Se) + \Phi^{-1}(Sp)) / \sqrt{2}]}, is 0.949.
#'   The default \code{AUC = 0.90} sits between the two: it is geometrically
#'   attainable, and being smaller than the binormal value it yields a
#'   larger Hanley-McNeil variance and hence a conservative sample size.
#'   The AUC gate is always evaluated at the \strong{expected} margins
#'   \eqn{n_d = \lfloor N \cdot prev \rfloor} and \eqn{n_{nd} = N - n_d},
#'   regardless of \code{design}: it is a deterministic criterion (see
#'   above), so there is no sampling distribution of \eqn{n_d} to average
#'   over in the first place. A candidate \code{N} that fails the AUC gate
#'   is \strong{not} assigned a Se/Sp joint probability of 0: it is simply
#'   never evaluated for Se/Sp at all (the Monte Carlo step is skipped for
#'   that \code{N}; see \code{auc_gate} in \code{Return}). "Blocked by the
#'   AUC gate" and "evaluated and found to have joint probability 0" are
#'   different things, and the return value keeps them distinguishable via
#'   \code{auc_gate$auc_pass} (\code{NA} for a candidate skipped for a
#'   different reason -- see \strong{No crossing found} below -- \code{FALSE}
#'   for one that failed the AUC gate, \code{TRUE} for one that passed it and
#'   so had its Se/Sp joint probability actually computed).
#'
#'   \strong{No crossing found.} If no \code{N} in \code{N_range} reaches
#'   \code{target_prob} for Se/Sp (either because the AUC gate blocked every
#'   candidate, or because it was cleared somewhere but the Se/Sp joint
#'   probability never reached \code{target_prob}), this function sets
#'   \code{n_total = NA_integer_}, \code{target_reached = FALSE},
#'   \code{joint_prob_se_sp = NA_real_} (and its alias \code{joint_prob}),
#'   \code{n_diseased = NA_integer_}, \code{n_non_diseased = NA_integer_},
#'   and reports the diagnostics \code{max_joint_prob_evaluated} (the highest
#'   Se/Sp joint probability actually computed anywhere in \code{N_range};
#'   \code{NA_real_} if the AUC gate blocked every candidate, since then no
#'   Se/Sp probability was ever computed at all) and \code{N_at_max_joint_prob}
#'   (the \code{N} at which that maximum was observed; \code{NA_integer_} in
#'   the same AUC-blocked-everywhere case), together with a warning naming
#'   which of the two situations occurred and how to search further (widen
#'   \code{N_range} or relax \code{delta_auc}). \code{n_total} is
#'   \strong{never} \code{max(N_range)} when no candidate reaches the target
#'   -- that would silently misrepresent an unreached target as a validated
#'   sample size, which versions \verb{<= 0.6.5} did (the defect corrected
#'   here, the same class of defect M-01 corrected in
#'   \code{\link{bam_sample_size}}).
#'
#'   \strong{The two sampling designs for Se and Sp differ in what is
#'   random, and therefore in the variance of the joint probability.}
#'
#'   \code{design = "cohort"} (default) --- \emph{prospective cohort /
#'   consecutive series}. Disease status is \strong{random}: each replicate
#'   draws \eqn{n_d \sim Bin(N, prev)} diseased subjects and
#'   \eqn{n_{nd} = N - n_d} non-diseased, then \eqn{X_{Se} \sim Bin(n_d, Se)}
#'   and \eqn{X_{Sp} \sim Bin(n_{nd}, Sp)} conditional on that draw. A
#'   replicate whose margins are degenerate (\eqn{n_d = 0} or
#'   \eqn{n_{nd} = 0}, so that no interval can be formed for that arm) is
#'   counted as a \strong{failure}, not discarded, so the reported joint
#'   probability is unconditional (denominator \code{B}). This is the same
#'   convention, and the same threshold, used by
#'   \code{\link{bam_sample_size}}: an arm of size \strong{one} is not
#'   degenerate. The Wilson interval used here (unlike the Wald interval)
#'   is well-defined and non-degenerate at \eqn{n = 1}: its width does not
#'   collapse to zero the way the Wald width does when \eqn{\hat p} is 0 or
#'   1, so a one-subject arm is scored on its actual (wide) Wilson width
#'   rather than being excluded by convention.
#'
#'   \code{design = "fixed"} --- \emph{fixed disease-status margins} (the
#'   behaviour of package versions <= 0.4.0). \eqn{n_d = \lfloor N \cdot
#'   prev \rfloor} and \eqn{n_{nd} = N - n_d} are treated as fixed by
#'   design, so \eqn{X_{Se} \sim Bin(n_d, Se)} and
#'   \eqn{X_{Sp} \sim Bin(n_{nd}, Sp)} condition on the expected margins
#'   rather than sampling them. This applies only to a design that recruits
#'   the two disease groups separately with pre-specified sizes (e.g. a
#'   case-control accuracy study), and it \strong{overstates the assurance}
#'   for a prospective cohort, where the number of diseased subjects
#'   actually enrolled is itself random. See \code{Note}.
#'
#' @param Se Expected sensitivity. Default 0.85.
#' @param Sp Expected specificity. Default 0.90.
#' @param AUC Expected AUC. Default 0.90. Must be geometrically compatible
#'   with \code{Se} and \code{Sp}: see \code{Details}. Raising \code{Se} or
#'   \code{Sp} raises the minimum attainable AUC, so a high-accuracy
#'   operating point may require \code{AUC} to be supplied explicitly.
#' @param delta_se Half-width target for Se. Default 0.07.
#' @param delta_sp Half-width target for Sp. Default 0.05.
#' @param delta_auc Half-width target for AUC. Default 0.05.
#' @param prev Disease prevalence. Default 0.20.
#' @param design Sampling design for Se and Sp, \code{"cohort"} (default) or
#'   \code{"fixed"}. \code{"cohort"} treats the number of diseased subjects
#'   as random, as in a prospective cohort or consecutive series, and is the
#'   appropriate choice for almost all diagnostic accuracy studies.
#'   \code{"fixed"} conditions on \eqn{n_d = \lfloor N \cdot prev \rfloor}
#'   diseased subjects and reproduces the behaviour of versions <= 0.4.0; it
#'   applies only to a design that recruits the two disease groups
#'   separately with pre-specified sizes, and it is \strong{not} valid for a
#'   prospective cohort. See \code{Details}.
#' @param target_prob Minimum joint probability \strong{for Se and Sp}.
#'   Default 0.80. See \code{Details}.
#' @param N_range Range of total N to search. Default \code{seq(100, 800, by = 10)}.
#' @param B MC replications. Default 5000. A warning is issued when
#'   \code{0 < B < 1000}, since the Monte Carlo error of the reported joint
#'   probability may then be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{n_total}{The first \code{N} in \code{N_range} (searched in the
#'       order given) achieving the joint target for Se and Sp, subject to
#'       the AUC gate. \strong{Not} the smallest integer \code{N} with that
#'       property unless \code{N_range} contains every integer in its span;
#'       see \code{Details}. \code{NA_integer_} when \code{target_reached =
#'       FALSE} (see \code{Details}, "No crossing found") -- never
#'       \code{max(N_range)}.}
#'     \item{target_reached}{Logical. \code{TRUE} if some candidate \code{N}
#'       in \code{N_range} reached \code{target_prob} for Se/Sp (subject to
#'       the AUC gate); \code{FALSE} otherwise, in which case \code{n_total},
#'       \code{n_diseased}, \code{n_non_diseased} and \code{joint_prob_se_sp}
#'       (and its alias \code{joint_prob}) are all \code{NA}, and
#'       \code{max_joint_prob_evaluated} / \code{N_at_max_joint_prob} carry
#'       the diagnostics instead. See \code{Details}, "No crossing found".}
#'     \item{n_diseased}{Number of diseased at \code{n_total}. Under
#'       \code{design = "cohort"} this is \eqn{\lfloor N \cdot prev
#'       \rfloor}, the \strong{expected} count, not a value fixed by
#'       design: the number actually diseased varies from one replicate,
#'       and one real study, to the next. \code{NA_integer_} when
#'       \code{target_reached = FALSE}.}
#'     \item{n_non_diseased}{Number of non-diseased at \code{n_total}. Same
#'       caveat as \code{n_diseased} under \code{design = "cohort"};
#'       \code{NA_integer_} when \code{target_reached = FALSE}.}
#'     \item{design}{The sampling design used for the Se/Sp Monte Carlo
#'       (\code{"cohort"} or \code{"fixed"}). Does not affect the AUC gate.}
#'     \item{joint_prob_se_sp}{Monte Carlo probability that the \strong{Se
#'       and Sp} intervals both meet their target width at \code{n_total}.
#'       This is \strong{not} a three-way joint probability: the AUC
#'       criterion is a deterministic median-style gate (see
#'       \code{Details}). \code{NA_real_} when \code{target_reached = FALSE}
#'       -- see \code{max_joint_prob_evaluated} for the diagnostic value in
#'       that case.}
#'     \item{joint_prob}{\strong{Deprecated} alias of
#'       \code{joint_prob_se_sp}, kept for backward compatibility with
#'       version 0.2.0. Use \code{joint_prob_se_sp}.}
#'     \item{joint_prob_mcse}{Monte Carlo standard error of
#'       \code{joint_prob_se_sp} at \code{n_total}, i.e.
#'       \code{sqrt(joint_prob_se_sp * (1 - joint_prob_se_sp) / B)}.
#'       \code{NA_real_} when \code{target_reached = FALSE} (no \code{n_total}
#'       to evaluate it at).}
#'     \item{search_type}{Always \code{"grid_first_candidate"}: \code{n_total}
#'       is the first candidate in \code{N_range}, searched in the order
#'       given, reaching \code{target_prob} -- see \code{Details}.}
#'     \item{N_range_used}{The candidate \code{N} values actually searched,
#'       i.e. \code{N_range} as supplied (every element is visited by the
#'       search, whether or not its Se/Sp probability was ultimately
#'       computed; see \code{auc_gate}).}
#'     \item{seed}{The \code{seed} argument actually used.}
#'     \item{auc_gate_passed}{Logical: whether any \code{N} in
#'       \code{N_range} satisfied the AUC precision criterion. Kept for
#'       backward compatibility; see \code{auc_gate} for the full detail.}
#'     \item{auc_gate}{A list describing the deterministic Hanley-McNeil AUC
#'       filter applied at every candidate \code{N}: \code{AUC} and
#'       \code{delta_auc} (the parameters used), \code{table} (a
#'       \code{data.frame} with one row per element of \code{N_range_used},
#'       columns \code{N}, \code{n_d_exp}, \code{n_nd_exp}, \code{auc_width}
#'       and \code{auc_pass}; \code{auc_pass} is \code{NA} for a candidate
#'       skipped because its \strong{expected} margin was already 0 -- i.e.
#'       the AUC gate itself was never evaluated for that \code{N} --
#'       \code{FALSE} for a candidate whose AUC gate was evaluated and
#'       failed, so no Se/Sp joint probability was computed for it either,
#'       and \code{TRUE} for a candidate that passed, meaning its Se/Sp
#'       joint probability WAS computed; see \code{Details}), and
#'       \code{first_N_auc_pass} (the first \code{N} in \code{N_range_used}
#'       with \code{auc_pass = TRUE}; \code{NA_integer_} if none).}
#'     \item{AUC_min}{The geometric lower bound on AUC implied by
#'       \code{Se} and \code{Sp}.}
#'     \item{max_joint_prob_evaluated}{The highest Se/Sp joint probability
#'       actually computed anywhere in \code{N_range} (i.e. at a candidate
#'       whose AUC gate passed). Always populated when any Se/Sp probability
#'       was computed, including when \code{target_reached = TRUE} (where it
#'       equals \code{joint_prob_se_sp}, since the search evaluates \code{N}
#'       in the given order and stops at the first crossing).
#'       \code{NA_real_} if the AUC gate blocked every candidate in
#'       \code{N_range} (no Se/Sp probability was ever computed; see
#'       \code{Details}).}
#'     \item{N_at_max_joint_prob}{The \code{N} (one element of
#'       \code{N_range_used}) at which \code{max_joint_prob_evaluated} was
#'       observed. Equals \code{n_total} when \code{target_reached = TRUE}.
#'       \code{NA_integer_} in the same AUC-blocked-everywhere case as
#'       \code{max_joint_prob_evaluated}.}
#'     \item{buderer_N}{Buderer-based total N for comparison.}
#'     \item{B}{The \code{B} argument actually used.}
#'   }
#' @note AUC CI width is computed using the Hanley-McNeil variance
#'   approximation evaluated at the assumed AUC, rather than from simulated
#'   Mann-Whitney statistics. This makes the AUC component deterministic per
#'   \code{N} and therefore a median-style, not an assurance-style,
#'   criterion. See \code{Details}.
#'
#'   \strong{Design matters, a lot, for the Se/Sp joint probability.} Under
#'   the default \code{design = "cohort"}, the number of diseased subjects
#'   is itself random, which adds a source of sampling variability that
#'   \code{design = "fixed"} conditions away. At a realistic adverse
#'   operating point (\code{Se = 0.70}, \code{Sp = 0.80}, \code{prev =
#'   0.20}, \code{delta_se = 0.08}, \code{delta_sp = 0.06}, \code{N =
#'   650}), the fixed-margin joint probability is \strong{0.807} while the
#'   cohort joint probability is \strong{0.716}: the fixed design
#'   overstates the assurance by roughly 9 percentage points. Sizing a
#'   prospective cohort with \code{design = "fixed"} therefore yields a
#'   sample size that is too small for the stated target. Versions <= 0.4.0
#'   offered only the fixed-margin behaviour.
#' @references
#' Hanley JA, McNeil BJ (1982). The meaning and use of the area under a
#' receiver operating characteristic (ROC) curve. \emph{Radiology}
#' 143:29-36.
#'
#' Buderer NMF (1996). Statistical methodology: I. Incorporating the
#' prevalence of disease into the sample size calculation for sensitivity
#' and specificity. \emph{Acad Emerg Med} 3:895-900.
#' \doi{10.1111/j.1553-2712.1996.tb03538.x}
#' @examples
#' result <- joint_sample_size(B = 1000, N_range = seq(100, 700, by = 20))
#' print(result)
#'
#' # An AUC below the geometric minimum is refused:
#' try(joint_sample_size(Se = 0.85, Sp = 0.90, AUC = 0.80))
#'
#' # Legacy fixed-margin design (not valid for a prospective cohort)
#' old <- joint_sample_size(B = 1000, design = "fixed",
#'                          N_range = seq(100, 700, by = 20))
#' @export
joint_sample_size <- function(Se = 0.85,
                              Sp = 0.90,
                              AUC = 0.90,
                              delta_se = 0.07,
                              delta_sp = 0.05,
                              delta_auc = 0.05,
                              prev = 0.20,
                              design = c("cohort", "fixed"),
                              target_prob = 0.80,
                              N_range = seq(100, 800, by = 10),
                              B = 5000,
                              seed = 2026) {
  # --- preserve the caller's RNG state (kind AND seed) ------------------
  # See save_rng_state()/restore_rng_state(): restoring only .Random.seed's
  # VALUE is not enough, because set.seed() called later by unrelated code
  # with no explicit `kind` argument reuses whichever kind is CURRENTLY
  # ACTIVE. The set.seed() call below names its kind explicitly
  # (Mersenne-Twister, R's own default), so the Se/Sp Monte Carlo
  # reproduces the same numbers regardless of the caller's own RNG
  # configuration.
  old_rng_state <- save_rng_state()
  on.exit(restore_rng_state(old_rng_state), add = TRUE)

  # Validate inputs
  design <- match.arg(design)
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1, AUC > 0.5, AUC <= 1)
  stopifnot(delta_se > 0, delta_sp > 0, delta_auc > 0)
  stopifnot(prev > 0, prev < 1)
  stopifnot(target_prob > 0, target_prob < 1)
  stopifnot(B >= 1)
  warn_small_B(B)

  # --- (a): AUC must be geometrically compatible with (Se, Sp) ---------
  # Minimum area of a concave ROC through the operating point (1 - Sp, Se):
  # the polygon (0,0) -> (1-Sp, Se) -> (1,1).
  AUC_min <- 0.5 * (1 - Sp) * Se + 0.5 * Sp * (Se + 1)
  if (AUC < AUC_min - 1e-9) {
    AUC_binormal <- stats::pnorm(
      (stats::qnorm(Se) + stats::qnorm(Sp)) / sqrt(2)
    )
    stop(
      sprintf(
        paste0(
          "AUC = %.3f is geometrically impossible for Se = %.3f and Sp = %.3f.\n",
          "  Any concave ROC curve through the operating point (1 - Sp, Se) = ",
          "(%.3f, %.3f)\n  has area at least AUC_min = %.4f ",
          "[= 0.5*(1-Sp)*Se + 0.5*Sp*(Se+1)].\n",
          "  For reference, the binormal AUC implied by this operating point ",
          "is %.4f.\n",
          "  Supply an AUC in [%.4f, 1], or lower Se / Sp."
        ),
        AUC, Se, Sp, 1 - Sp, Se, AUC_min, AUC_binormal, AUC_min
      ),
      call. = FALSE
    )
  }

  target_se_width <- 2 * delta_se
  target_sp_width <- 2 * delta_sp
  target_auc_width <- 2 * delta_auc

  # N_range_used: every candidate the search visits, in the order supplied
  # (point 1: same order of candidates as before; N_range is never sorted
  # here, matching the pre-existing loop behaviour).
  N_range_used <- N_range
  n_candidates <- length(N_range_used)

  optimal_N <- NA_integer_
  joint_prob_achieved <- NA_real_
  joint_prob_mcse_achieved <- NA_real_
  # (c): distinguish "never computed" from "computed and equal to 0".
  joint_prob <- NA_real_
  auc_gate_passed <- FALSE
  target_reached <- FALSE
  max_joint_prob_evaluated <- NA_real_
  N_at_max_joint_prob <- NA_integer_

  # auc_gate$table: one row per candidate in N_range_used, recording the
  # deterministic AUC gate's inputs/outcome at that N. auc_pass is NA for a
  # candidate whose EXPECTED margin was already 0 (the AUC gate itself was
  # never evaluated there), FALSE for one the AUC gate evaluated and
  # rejected (no Se/Sp joint probability computed for it), TRUE for one
  # that passed (Se/Sp joint probability WAS computed). See @return.
  auc_table_N <- integer(n_candidates)
  auc_table_n_d_exp <- integer(n_candidates)
  auc_table_n_nd_exp <- integer(n_candidates)
  auc_table_width <- rep(NA_real_, n_candidates)
  auc_table_pass <- rep(NA, n_candidates)

  for (idx in seq_len(n_candidates)) {
    N <- N_range_used[idx]
    auc_table_N[idx] <- as.integer(N)

    # Expected margins: used for the AUC gate always, and for the Se/Sp
    # Monte Carlo under design = "fixed" (see below).
    n_d_exp <- floor(N * prev)
    n_nd_exp <- N - n_d_exp
    auc_table_n_d_exp[idx] <- as.integer(n_d_exp)
    auc_table_n_nd_exp[idx] <- as.integer(n_nd_exp)
    # Same degeneracy threshold as the per-replicate check below and as
    # bam_sample_size(): skip only a candidate N whose EXPECTED margin is
    # zero (undefined for both hanley_mcneil_var() and wilson_width()). An
    # expected margin of 1 is not skipped; see @details. auc_table_pass
    # stays NA here: the AUC gate was never even evaluated for this N.
    if (n_d_exp == 0 || n_nd_exp == 0) next

    # --- AUC: DETERMINISTIC gate via Hanley-McNeil (median-style) ---
    # Evaluated at the expected margins regardless of design: this
    # criterion has no sampling distribution to average over (see
    # @details), so there is nothing for "cohort" to change here.
    var_auc <- hanley_mcneil_var(AUC, n_d_exp, n_nd_exp)
    auc_width <- 2 * stats::qnorm(0.975) * sqrt(var_auc)
    auc_pass <- auc_width <= target_auc_width
    auc_table_width[idx] <- auc_width
    auc_table_pass[idx] <- auc_pass

    # If AUC alone fails, skip MC for Se/Sp: a candidate blocked by the AUC
    # gate is NOT assigned a Se/Sp joint probability of 0, it is simply
    # never evaluated for Se/Sp (see @details).
    if (!auc_pass) next
    auc_gate_passed <- TRUE

    # --- Se and Sp: Monte Carlo with Wilson CI ---
    set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

    if (design == "cohort") {
      # Disease status is RANDOM in a prospective cohort: the number of
      # diseased subjects actually enrolled varies from study to study.
      # Conditioning on its expected value (design = "fixed") overstates
      # the assurance; see @details and @note.
      n_d_b <- stats::rbinom(B, N, prev)
      n_nd_b <- N - n_d_b
      # Degenerate iff an arm received NO subjects at all (n = 0), the
      # point at which the Wilson width is literally undefined (division
      # by n). An arm of size 1 is NOT degenerate: the Wilson interval is
      # well-defined there and its width does not collapse to zero (unlike
      # the Wald width), so it is scored on its own (wide) width. This is
      # the same threshold used by bam_sample_size() (see its @details).
      degenerate <- (n_d_b == 0L) | (n_nd_b == 0L)

      # Degenerate replicates count as FAILURES, not exclusions (the
      # denominator stays B); clamp their margins to 1 only so that the
      # vectorised draws below stay well-defined, since `pass` overrides
      # them to FALSE regardless of the simulated width.
      safe_n_d <- pmax(n_d_b, 1L)
      safe_n_nd <- pmax(n_nd_b, 1L)
      x_se <- stats::rbinom(B, safe_n_d, Se)
      x_sp <- stats::rbinom(B, safe_n_nd, Sp)

      se_width <- wilson_width(x_se, safe_n_d)
      sp_width <- wilson_width(x_sp, safe_n_nd)

      pass <- !degenerate &
        (se_width <= target_se_width) & (sp_width <= target_sp_width)
    } else {
      # design == "fixed": disease-status margins fixed by design, valid
      # only when the two groups are recruited separately with
      # pre-specified sizes; NOT valid for a prospective cohort.
      x_se <- stats::rbinom(B, n_d_exp, Se)
      x_sp <- stats::rbinom(B, n_nd_exp, Sp)

      se_width <- wilson_width(x_se, n_d_exp)
      sp_width <- wilson_width(x_sp, n_nd_exp)

      pass <- (se_width <= target_se_width) & (sp_width <= target_sp_width)
    }

    # Joint over Se and Sp ONLY (the AUC gate is deterministic, see @details)
    joint_prob <- mean(pass)

    if (is.na(max_joint_prob_evaluated) || joint_prob > max_joint_prob_evaluated) {
      max_joint_prob_evaluated <- joint_prob
      N_at_max_joint_prob <- as.integer(N)
    }

    if (joint_prob >= target_prob) {
      optimal_N <- as.integer(N)
      joint_prob_achieved <- joint_prob
      joint_prob_mcse_achieved <- sqrt(joint_prob * (1 - joint_prob) / B)
      target_reached <- TRUE
      break
    }
  }

  # auc_gate: always populated, whether or not the target was reached (see
  # @return). first_N_auc_pass is the first candidate (in N_range_used
  # order) whose AUC gate passed, regardless of what its Se/Sp joint
  # probability turned out to be.
  auc_table <- data.frame(
    N = auc_table_N,
    n_d_exp = auc_table_n_d_exp,
    n_nd_exp = auc_table_n_nd_exp,
    auc_width = auc_table_width,
    auc_pass = auc_table_pass
  )
  first_pass_idx <- which(auc_table_pass %in% TRUE)
  first_N_auc_pass <- if (length(first_pass_idx) > 0) {
    auc_table_N[first_pass_idx[1]]
  } else {
    NA_integer_
  }
  auc_gate <- list(
    AUC = AUC,
    delta_auc = delta_auc,
    table = auc_table,
    first_N_auc_pass = first_N_auc_pass
  )

  if (!target_reached) {
    # M-01-style contract (bam_sample_size()): no longer silently returns
    # max(N_range) as if it were a validated solution (the defect this
    # corrects; see @details, "No crossing found"). n_total / n_diseased /
    # n_non_diseased / joint_prob_se_sp are NA; max_joint_prob_evaluated and
    # N_at_max_joint_prob report what the search actually saw.
    if (!auc_gate_passed) {
      # (c): the AUC gate blocked every candidate N, so the Se/Sp Monte
      # Carlo probability was NEVER computed. Reporting 0 here (the old
      # initialisation value) would be indistinguishable from a genuine 0.
      warning("The AUC precision target (delta_auc = ", delta_auc,
              ") was not met at any N in N_range, so no N was ever ",
              "evaluated for Se/Sp. n_total and joint_prob_se_sp are NA ",
              "(target_reached = FALSE; joint_prob_se_sp was never ",
              "computed, not 0). Expand N_range or relax delta_auc.",
              call. = FALSE)
    } else {
      warning("No N in N_range achieved the target joint probability for ",
              "Se and Sp: the highest joint_prob_se_sp seen was ",
              sprintf("%.6f", max_joint_prob_evaluated), " at N = ",
              N_at_max_joint_prob, ", against target_prob = ", target_prob,
              ". n_total and joint_prob_se_sp are NA (target_reached = ",
              "FALSE). Expand N_range; see max_joint_prob_evaluated / ",
              "N_at_max_joint_prob for diagnostics.",
              call. = FALSE)
    }
  }

  # Expected margins at n_total. Under design = "cohort" these are EXPECTED
  # counts, not values fixed by the study design; see @return. NA when the
  # target was not reached: there is no "final N" to report margins for.
  if (target_reached) {
    n_d_final <- as.integer(floor(optimal_N * prev))
    n_nd_final <- as.integer(optimal_N - n_d_final)
  } else {
    n_d_final <- NA_integer_
    n_nd_final <- NA_integer_
  }

  # Buderer comparison
  buderer_N <- buderer_total_N(Se, Sp, delta_se, delta_sp, prev)

  structure(
    list(
      method = "Joint Sample Size for Se + Sp + AUC",
      design = design,
      n_total = optimal_N,
      target_reached = target_reached,
      n_diseased = n_d_final,
      n_non_diseased = n_nd_final,
      joint_prob_se_sp = joint_prob_achieved,
      # Deprecated alias, kept so that 0.2.0 code keeps running.
      joint_prob = joint_prob_achieved,
      joint_prob_mcse = joint_prob_mcse_achieved,
      search_type = "grid_first_candidate",
      N_range_used = as.integer(N_range_used),
      seed = seed,
      auc_gate_passed = auc_gate_passed,
      auc_gate = auc_gate,
      max_joint_prob_evaluated = max_joint_prob_evaluated,
      N_at_max_joint_prob = N_at_max_joint_prob,
      buderer_N = buderer_N,
      Se = Se,
      Sp = Sp,
      AUC = AUC,
      AUC_min = AUC_min,
      prev = prev,
      B = B,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
