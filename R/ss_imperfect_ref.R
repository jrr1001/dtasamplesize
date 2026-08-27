
# --- internal helpers: exact identification of (pi, Se, Sp) under a known,
# conditionally-independent imperfect reference standard -------------------
#
# Notation: T = index test, R = reference standard, D = true disease status
# (unobserved). Given Se_ref = P(R+|D+) and Sp_ref = P(R-|D-) *known* and
# T independent of R given D, the joint distribution of the observed 2x2
# table (T x R) is a smooth, invertible function of theta = (pi, Se, Sp):
#   p11 = P(T+,R+) = pi*Se*Se_ref       + (1-pi)*(1-Sp)*(1-Sp_ref)
#   p10 = P(T+,R-) = pi*Se*(1-Se_ref)   + (1-pi)*(1-Sp)*Sp_ref
#   p01 = P(T-,R+) = pi*(1-Se)*Se_ref   + (1-pi)*Sp*(1-Sp_ref)
#   p00 = P(T-,R-) = pi*(1-Se)*(1-Se_ref) + (1-pi)*Sp*Sp_ref
# so theta is exactly identified from the cell probabilities, with
#   Se_hat = (Sp_ref*p11 - (1-Sp_ref)*p10) / (p11 + p01 - (1-Sp_ref))
#   Sp_hat = (Se_ref*p00 - (1-Se_ref)*p01) / (p00 + p10 - (1-Se_ref))
# (both reduce to the naive proportion when Se_ref = Sp_ref = 1). These
# functions are not exported; they exist only to build the sample-size
# formulas below, and are validated against an independent from-scratch
# Monte Carlo / Fisher-information derivation (see package NEWS / tests).

irs_cell_probs <- function(Se, Sp, prev, Se_ref, Sp_ref) {
  p11 <- prev * Se * Se_ref + (1 - prev) * (1 - Sp) * (1 - Sp_ref)
  p10 <- prev * Se * (1 - Se_ref) + (1 - prev) * (1 - Sp) * Sp_ref
  p01 <- prev * (1 - Se) * Se_ref + (1 - prev) * Sp * (1 - Sp_ref)
  p00 <- prev * (1 - Se) * (1 - Se_ref) + (1 - prev) * Sp * Sp_ref
  c(p11 = p11, p10 = p10, p01 = p01, p00 = p00)
}

# Gradient of Se_hat / Sp_hat w.r.t. (p11, p10, p01, p00), evaluated at the
# planning values of (Se, Sp, prev). J = Se_ref + Sp_ref - 1 is the Youden
# index of the reference standard; both gradients scale as 1/J, which is
# why a weak reference (small J) inflates the corrected-estimand variance.
irs_grad_se <- function(Se, Sp, prev, Se_ref, Sp_ref) {
  J <- Se_ref + Sp_ref - 1
  c(Sp_ref - Se, -(1 - Sp_ref), -Se, 0) / (prev * J)
}
irs_grad_sp <- function(Se, Sp, prev, Se_ref, Sp_ref) {
  J <- Se_ref + Sp_ref - 1
  c(0, -Sp, -(1 - Se_ref), Se_ref - Sp) / ((1 - prev) * J)
}

# N * Var(g'p_hat) by the delta method, for p_hat the multinomial cell
# proportions from N observations: Var(p_hat) = (diag(p) - p p') / N, so
# N * Var(g'p_hat) = g'(diag(p) - p p')g = sum(p*g^2) - (sum(p*g))^2.
irs_delta_nvar <- function(g, p) sum(p * g^2) - (sum(p * g))^2

# Apparent quantities: Se_apparent = P(T+|R+), Sp_apparent = P(T-|R-), and
# the marginal probabilities of testing reference-positive/negative. These
# are the quantities a naive analysis (no misclassification correction)
# actually estimates, and they have an exact closed form -- no delta
# method or inflation factor needed.
irs_apparent <- function(Se, Sp, prev, Se_ref, Sp_ref) {
  p <- irs_cell_probs(Se, Sp, prev, Se_ref, Sp_ref)
  P_pos <- unname(p["p11"] + p["p01"])
  P_neg <- unname(p["p00"] + p["p10"])
  list(Se_apparent = unname(p["p11"]) / P_pos,
       Sp_apparent = unname(p["p00"]) / P_neg,
       P_ref_pos = P_pos,
       P_ref_neg = P_neg,
       p = p)
}

#' Sample Size for a Diagnostic Accuracy Study With an Imperfect Reference
#' Standard
#'
#' Computes the sample size needed to estimate the sensitivity and
#' specificity of an index test to a target precision when the reference
#' standard used to verify disease status is itself imperfect (known
#' \code{Se_ref}, \code{Sp_ref}), under conditional independence of the
#' index test and the reference standard given true disease status. Two
#' different target quantities (\code{estimand}) can be sized for; both are
#' always reported (see Details and Value).
#'
#' @details
#' \strong{Two estimands, two different sample sizes.} An imperfect
#' reference standard forces a choice of what quantity the study is
#' actually powered to estimate precisely:
#' \itemize{
#'   \item \code{estimand = "apparent"} (default): the study targets
#'     precision for the \emph{apparent} sensitivity and specificity,
#'     \eqn{Se_{app} = P(T+\,|\,R+)} and \eqn{Sp_{app} = P(T-\,|\,R-)},
#'     i.e. what a standard 2x2-table analysis against the imperfect
#'     reference actually estimates. Under conditional independence of the
#'     index test \eqn{T} and the reference \eqn{R} given true disease
#'     status \eqn{D}, \eqn{Se_{app}} and \eqn{P(R+)} have an exact closed
#'     form (see below), so sizing for this estimand needs \strong{no
#'     inflation factor at all}: it is the ordinary Buderer calculation
#'     with \eqn{Se_{app}}/\eqn{Sp_{app}} in place of \eqn{Se}/\eqn{Sp} and
#'     \eqn{P(R+)}/\eqn{P(R-)} in place of prevalence.
#'   \item \code{estimand = "corrected"}: the study targets precision for
#'     the index test's \emph{true} sensitivity and specificity, recovered
#'     from the observed 2x2 table via the closed-form, exactly-identified
#'     misclassification correction (below). This is the quantity a
#'     bias-corrected analysis or a latent-class model would report, and
#'     it removes the verification bias -- but its sampling variance is
#'     driven by the whole 2x2 table, not by a ref-positive/negative
#'     subsample, and can be \strong{far} larger than either the naive
#'     Buderer variance or the apparent-estimand variance (see the worked
#'     example below). Both \eqn{N} are always computed and returned
#'     (elements \code{N_apparent}/\code{N_apparent_loss} and
#'     \code{N_corrected}/\code{N_corrected_loss}); \code{estimand} only
#'     selects which one is echoed in the generic \code{N_adjusted} /
#'     \code{n_total} slots, so a caller who reads only those cannot
#'     silently under-power a bias-corrected analysis by relying on the
#'     default.
#' }
#'
#' \strong{Identification.} Given \code{Se_ref} and \code{Sp_ref} known and
#' \eqn{T \perp R \mid D}, the model \eqn{(\pi, Se, Sp)} is exactly
#' identified from the observed cell probabilities
#' \eqn{p_{11}, p_{10}, p_{01}, p_{00}} of the \eqn{T \times R} table:
#' \deqn{\widehat{Se} = \frac{Sp_{ref}\,p_{11} - (1-Sp_{ref})\,p_{10}}
#'   {p_{11} + p_{01} - (1-Sp_{ref})}}
#' \deqn{\widehat{Sp} = \frac{Se_{ref}\,p_{00} - (1-Se_{ref})\,p_{01}}
#'   {p_{00} + p_{10} - (1-Se_{ref})}}
#' Its sampling variance is obtained by the delta method applied to the
#' cell-count multinomial (\code{N * Var(theta_hat) = g' (diag(p) - p p') g}
#' for the appropriate gradient \code{g}); this matches an independent
#' Fisher-information derivation and Monte Carlo simulation to within
#' Monte Carlo error. \eqn{Se_{app} = p_{11}/(p_{11}+p_{01})} and
#' \eqn{P(R+) = p_{11}+p_{01}} follow directly from the same model.
#'
#' \strong{Worked example (package defaults with \code{prev = 0.20}).}
#' At \code{Se = .85, Sp = .90, prev = .20, Se_ref = .90, Sp_ref = .95,
#' d_se = .07, d_sp = .05}: the classical Buderer \eqn{N} (ignoring
#' misclassification entirely) is 500; sizing for the apparent estimand
#' gives \eqn{N \approx} 732; sizing for the corrected estimand gives
#' \eqn{N \approx} 1235. A single scalar "variance inflation factor"
#' cannot represent both of these at once, because they answer different
#' questions -- which is why \code{VIF} (below) is retained only as an
#' informative, prevalence-specific quantity and no longer multiplies any
#' sample size.
#'
#' \strong{The retired \code{VIF}.} \eqn{VIF = 1/J_{ref}^2}, with
#' \eqn{J_{ref} = Se_{ref} + Sp_{ref} - 1} the Youden index of the
#' reference standard, is the classical Rogan-Gladen variance-inflation
#' factor for a \emph{prevalence} estimated from an imperfect screening
#' test (Rogan & Gladen 1978). It is still computed and reported (element
#' \code{VIF}) because it is a correct and citable quantity in its own
#' right, but versions <= 0.5.0 of this function multiplied the Buderer
#' \eqn{Se}/\eqn{Sp} sample size by it, which silently substitutes the
#' prevalence problem's variance ratio for the (different) Se/Sp problem's
#' variance ratio. Measured against the exact formulas above, that
#' substitution is wrong by up to two orders of magnitude for the
#' \code{estimand = "corrected"} sample size, and by a smaller but still
#' material amount for \code{estimand = "apparent"} (which, unlike the old
#' \code{VIF} path, needs no inflation factor to begin with). Do not
#' multiply any sample size by \code{VIF}; it is kept purely as a
#' diagnostic of how informative the reference standard is about
#' prevalence. \code{multiplier_se_corrected}, \code{multiplier_sp_corrected},
#' \code{multiplier_se_apparent} and \code{multiplier_sp_apparent} report
#' the \emph{actual} ratio of the real Se/Sp variance to Buderer's naive
#' variance, for comparison against \code{VIF}.
#'
#' \strong{\code{min_youden} guards near-non-identification, not the size
#' of the multiplier.} Both gradients above scale as \code{1/J_ref}, so a
#' reference standard with \eqn{J_{ref}} near 0 is only weakly informative
#' about \eqn{(Se, Sp)} and both corrected-estimand variances blow up; this
#' is what \code{min_youden} (default 0.5) guards against, exactly as
#' before. What changed is what passing the guard tells you: previously
#' the documentation claimed \code{min_youden = 0.5} "caps the variance
#' inflation factor at 4", which is no longer (and, given the defect this
#' version fixes, was never accurately) a bound on the sample-size
#' multiplier actually realised for either estimand. It is not: a
#' reference standard with \eqn{J_{ref} = 0.65} (comfortably above the
#' default threshold) can still require a \code{"corrected"}-estimand
#' multiplier above 100x at low prevalence, because that multiplier also
#' depends on prevalence and on the true Se/Sp, not on \eqn{J_{ref}} alone.
#' \code{min_youden} is therefore a floor against outright non-identification,
#' not a usability guarantee; check \code{multiplier_se_corrected} /
#' \code{multiplier_sp_corrected} (or just the resulting \code{N_corrected})
#' for that.
#'
#' @param Se Expected sensitivity of index test. Default 0.85.
#' @param Sp Expected specificity of index test. Default 0.90.
#' @param d_se Desired precision for Se (half-width). Default 0.07.
#' @param d_sp Desired precision for Sp (half-width). Default 0.05.
#' @param prev Disease prevalence. Default 0.30.
#' @param Se_ref Sensitivity of reference standard. Default 0.90.
#' @param Sp_ref Specificity of reference standard. Default 0.95.
#' @param loss_rate Expected loss-to-follow-up rate. Default 0.10.
#' @param alpha Significance level. Default 0.05.
#' @param estimand Which target quantity to size the primary \code{N} for:
#'   \code{"apparent"} (default) sizes for \eqn{P(T+\,|\,R+)} /
#'   \eqn{P(T-\,|\,R-)}, the quantity a naive (uncorrected) analysis
#'   against the imperfect reference estimates; \code{"corrected"} sizes
#'   for the true, misclassification-corrected \eqn{Se}/\eqn{Sp}. Both are
#'   always computed and returned regardless of this choice (see Value);
#'   \code{estimand} only selects which pair populates the generic
#'   \code{N_adjusted}/\code{n_total} slots. See Details.
#' @param B MC replications for validation. Default 5000. Set to 0 to
#'   skip MC validation. A warning is issued when \code{0 < B < 1000},
#'   since the Monte Carlo error of the validation probabilities may then
#'   be substantial; silence it with
#'   \code{options(dtasamplesize.warn_small_B = FALSE)}.
#' @param seed Random seed. Default 2026. The RNG state of the calling
#'   session is restored on exit.
#' @param sensitivity_table Logical. If \code{TRUE} (default), compute a
#'   table of \code{N_apparent}/\code{N_corrected} (and the informative
#'   \code{VIF}) for varying \code{Se_ref} and \code{Sp_ref}.
#' @param min_youden Minimum acceptable Youden index of the reference
#'   standard, \eqn{Se_{ref} + Sp_{ref} - 1}. Default 0.5. A reference
#'   standard below this is too weak for either estimand's correction to
#'   be usable, and the function stops with an informative error. This
#'   guards near-non-identification; it does \strong{not} bound how large
#'   the realised sample-size multiplier can be above the threshold. See
#'   \code{Details}. Lower it deliberately (and expect a very large
#'   \code{N}) if you really intend to plan against such a reference.
#' @param max_mc_cells Upper bound on \code{B * N_adjusted}, the number of
#'   cells the Monte Carlo validation allocates. Default \code{2e7}
#'   (roughly a few hundred MB). Exceeding it stops with an informative
#'   error suggesting a smaller \code{B} or \code{B = 0}. This guard also
#'   catches an extreme \code{prev}, which for \code{estimand = "corrected"}
#'   can inflate \code{N_adjusted} by several orders of magnitude even with
#'   a fine reference standard (see Details).
#' @param delta_se Alias for \code{d_se}, matching the \code{delta_*}
#'   naming used elsewhere in the package. If supplied (non-\code{NULL}),
#'   it takes precedence over \code{d_se}. Default \code{NULL} (use
#'   \code{d_se}).
#' @param delta_sp Alias for \code{d_sp}. If supplied (non-\code{NULL}),
#'   it takes precedence over \code{d_sp}. Default \code{NULL} (use
#'   \code{d_sp}).
#' @return Object of class \code{"dtasamplesize"} with additional elements:
#'   \describe{
#'     \item{estimand}{The \code{estimand} argument actually used.}
#'     \item{n_diseased, n_total}{Generic slots (read by
#'       \code{print.dtasamplesize} and by \code{\link{ss_unified}}):
#'       alias of whichever \code{estimand} was selected.}
#'     \item{n_diseased_unadjusted, n_nondiseased_unadjusted, N_unadjusted,
#'       N_buderer}{The classical Buderer numbers, i.e. as if the
#'       reference standard were perfect. \code{N_buderer} is an alias of
#'       \code{N_unadjusted}.}
#'     \item{n_diseased_adjusted, n_nondiseased_adjusted, N_adjusted,
#'       N_adjusted_loss}{Aliases of the \code{estimand}-selected numbers
#'       below (apparent or corrected), kept for backward compatibility
#'       with code that reads these generic names.}
#'     \item{Se_apparent, Sp_apparent, P_ref_pos, P_ref_neg,
#'       n_refpos_apparent, n_refneg_apparent, N_apparent,
#'       N_apparent_loss}{The \code{"apparent"}-estimand quantities:
#'       \eqn{P(T+\,|\,R+)}, \eqn{P(T-\,|\,R-)}, the marginal probabilities
#'       of a reference-positive/negative result, the required
#'       reference-positive/negative subsample sizes, and the resulting
#'       total \code{N} (with and without \code{loss_rate}).}
#'     \item{NVar_se_corrected, NVar_sp_corrected, n_se_corrected,
#'       n_sp_corrected, N_corrected, N_corrected_loss}{The
#'       \code{"corrected"}-estimand quantities: \code{N * Var} of the
#'       misclassification-corrected Se/Sp (the per-arm sample size these
#'       imply, which -- unlike the apparent/Buderer arms -- is already a
#'       total-\code{N} quantity, not a diseased/ref-positive subsample
#'       size), and the resulting total \code{N} (with and without
#'       \code{loss_rate}).}
#'     \item{youden_ref}{Youden index of the reference standard,
#'       \code{Se_ref + Sp_ref - 1}.}
#'     \item{VIF}{The Rogan-Gladen \strong{prevalence} variance-inflation
#'       factor, \code{1 / youden_ref^2}. Informative only -- see Details.
#'       Does not multiply any sample size here.}
#'     \item{multiplier_se_corrected, multiplier_sp_corrected,
#'       multiplier_se_apparent, multiplier_sp_apparent}{The actual ratio
#'       of the real Se/Sp variance (corrected or apparent) to Buderer's
#'       naive variance, for comparison against \code{VIF}.}
#'     \item{sensitivity_table}{Data frame (if requested) with
#'       \code{N_apparent} and \code{N_corrected} (and the informative
#'       \code{VIF}/\code{inflation_factor}) over a grid of \code{Se_ref}
#'       and \code{Sp_ref}.}
#'     \item{mc_validation}{Data frame (if \code{B > 0}) with, for the
#'       unadjusted (Buderer) and adjusted (\code{estimand}-selected)
#'       sample sizes: the probability that the \strong{apparent}
#'       sensitivity CI achieves the target width (\code{P_width_target}),
#'       the mean apparent sensitivity \eqn{P(T+\,|\,R+)}
#'       (\code{se_apparent}) and its \code{bias} against the true Se, and
#'       the mean \strong{corrected} sensitivity estimator
#'       (\code{se_corrected}) and its \code{bias_corrected} -- which
#'       should be close to zero, in contrast to \code{bias}.}
#'     \item{results}{Data frame summarising both estimands side by side
#'       (\code{N}, \code{N_with_loss}, and which one \code{estimand}
#'       selected), printed automatically by \code{print.dtasamplesize}.}
#'   }
#' @references
#' Rogan WJ, Gladen B (1978). Estimating prevalence from the results of a
#' screening test. \emph{Am J Epidemiol} 107:71-76.
#' \doi{10.1093/oxfordjournals.aje.a112510}
#'
#' Buderer NMF (1996). Statistical methodology: I. Incorporating the
#' prevalence of disease into the sample size calculation for sensitivity
#' and specificity. \emph{Acad Emerg Med} 3:895-900.
#' \doi{10.1111/j.1553-2712.1996.tb03538.x}
#' @note \strong{Misclassification bias persists under \code{"apparent"}.}
#'   Sizing for the apparent estimand gives an exact, uninflated CI for
#'   \eqn{P(T+\,|\,R+)} -- but that quantity is \strong{not} the index
#'   test's true sensitivity when the reference standard is imperfect.
#'   \eqn{P(T+\,|\,R+)} converges to a value below the true Se whenever
#'   \code{Sp_ref < 1} (reported as \code{bias} in \code{mc_validation});
#'   precision under \code{"apparent"} does not remove this bias. Recovering
#'   an unbiased estimate of the true Se/Sp requires either sizing for
#'   \code{estimand = "corrected"} and applying the corresponding
#'   misclassification correction at the analysis stage (or an equivalent
#'   latent-class analysis), which needs the (typically much larger)
#'   \code{N_corrected}, or accepting \eqn{P(T+\,|\,R+)} as the reported
#'   estimand and being explicit that it is not the true Se. Both
#'   assume conditional independence between the index test and the
#'   reference standard given true disease status; this may be violated
#'   when both measure the same construct.
#' @examples
#' result <- ss_imperfect_ref(B = 0)
#' print(result)
#' @export
ss_imperfect_ref <- function(Se = 0.85,
                             Sp = 0.90,
                             d_se = 0.07,
                             d_sp = 0.05,
                             prev = 0.30,
                             Se_ref = 0.90,
                             Sp_ref = 0.95,
                             loss_rate = 0.10,
                             alpha = 0.05,
                             estimand = c("apparent", "corrected"),
                             B = 5000,
                             seed = 2026,
                             sensitivity_table = TRUE,
                             min_youden = 0.5,
                             max_mc_cells = 2e7,
                             delta_se = NULL,
                             delta_sp = NULL) {
  # delta_se/delta_sp are aliases of d_se/d_sp, matching the delta_*
  # naming used elsewhere in the package; when supplied they take
  # precedence, but default behaviour (both NULL) is unchanged.
  if (!is.null(delta_se)) d_se <- delta_se
  if (!is.null(delta_sp)) d_sp <- delta_sp

  estimand <- match.arg(estimand)

  # --- preserve the caller's RNG state (kind AND seed) ------------------
  # See save_rng_state()/restore_rng_state(): restoring only .Random.seed's
  # VALUE (as versions <= 0.6.1 did here) is not enough, because
  # set.seed() called later by unrelated code with no explicit `kind`
  # argument reuses whichever kind is CURRENTLY ACTIVE. The set.seed()
  # call below names all three kinds explicitly (Mersenne-Twister /
  # Inversion / Rejection, R's own defaults), so the reported apparent
  # Se/Sp and mc_validation figures -- including the manuscript's
  # published apparent sensitivity of 0.764 (bias -0.086) -- reproduce the
  # same numbers regardless of the caller's own RNG configuration. Before
  # this fix, that 0.764 moved to 0.765 under Wichmann-Hill for otherwise
  # identical arguments and seed.
  old_rng_state <- save_rng_state()
  on.exit(restore_rng_state(old_rng_state), add = TRUE)

  # Validate inputs
  stopifnot(Se > 0, Se < 1, Sp > 0, Sp < 1)
  stopifnot(d_se > 0, d_sp > 0)
  stopifnot(prev > 0, prev < 1)
  stopifnot(Se_ref > 0, Se_ref <= 1, Sp_ref > 0, Sp_ref <= 1)
  stopifnot(Se_ref + Sp_ref > 1)  # must be informative
  stopifnot(loss_rate >= 0, loss_rate < 1)
  stopifnot(min_youden > 0, min_youden <= 1)
  stopifnot(max_mc_cells > 0)
  warn_small_B(B)

  z_alpha <- stats::qnorm(1 - alpha / 2)

  # Youden index of the reference standard and the (informative-only,
  # prevalence-specific) Rogan-Gladen VIF -- see Details for why this no
  # longer multiplies any sample size.
  youden_ref <- Se_ref + Sp_ref - 1
  VIF <- 1 / youden_ref^2

  # --- refuse a reference standard so weak that Se/Sp are barely
  # identified (both gradients below scale as 1/youden_ref) -----------
  if (youden_ref < min_youden) {
    stop(
      sprintf(
        paste0(
          "Reference standard too weak: Youden index = Se_ref + Sp_ref - 1 = ",
          "%.4f (< min_youden = %.4f).\n",
          "  This index is the denominator of both the Rogan-Gladen ",
          "prevalence VIF (1/Youden^2 = %.1f; informational only, see ",
          "'VIF' in the return value) and of the misclassification-",
          "correction gradients used by estimand = \"corrected\" (which ",
          "scale as 1/(prev*Youden) and 1/((1-prev)*Youden)); a reference ",
          "standard this weak leaves the corrected Se/Sp only weakly ",
          "identified.\n",
          "  Note: clearing this guard does NOT bound how large the ",
          "required sample size can get for estimand = \"corrected\" -- ",
          "that also depends on prevalence and on the true Se/Sp, and can ",
          "still run into the hundreds-of-times range at low prevalence ",
          "even with a comfortable Youden index; check ",
          "multiplier_se_corrected / multiplier_sp_corrected. See Details.\n",
          "  Use a better reference standard, or lower min_youden ",
          "deliberately if you really intend to plan against this one ",
          "(and expect a very large N)."
        ),
        youden_ref, min_youden, VIF
      ),
      call. = FALSE
    )
  }

  # --- classical Buderer baseline: as if the reference were perfect ---
  n_unadj_se <- buderer_n(Se, d_se, alpha)
  n_unadj_sp <- buderer_n(Sp, d_sp, alpha)
  N_unadj <- buderer_total_N(Se, Sp, d_se, d_sp, prev, alpha)

  buderer_nvar_se <- Se * (1 - Se) / prev
  buderer_nvar_sp <- Sp * (1 - Sp) / (1 - prev)

  # --- estimand = "apparent": exact closed form, no inflation factor --
  app <- irs_apparent(Se, Sp, prev, Se_ref, Sp_ref)
  Se_apparent <- app$Se_apparent
  Sp_apparent <- app$Sp_apparent
  P_ref_pos <- app$P_ref_pos
  P_ref_neg <- app$P_ref_neg

  n_refpos_apparent <- ceiling(z_alpha^2 * Se_apparent * (1 - Se_apparent) / d_se^2)
  n_refneg_apparent <- ceiling(z_alpha^2 * Sp_apparent * (1 - Sp_apparent) / d_sp^2)
  N_apparent <- ceiling(max(n_refpos_apparent / P_ref_pos,
                             n_refneg_apparent / P_ref_neg))
  N_apparent_loss <- ceiling(N_apparent / (1 - loss_rate))

  multiplier_se_apparent <- (Se_apparent * (1 - Se_apparent) / P_ref_pos) /
    buderer_nvar_se
  multiplier_sp_apparent <- (Sp_apparent * (1 - Sp_apparent) / P_ref_neg) /
    buderer_nvar_sp

  # --- estimand = "corrected": delta-method variance of the exactly-
  # identified misclassification-corrected Se/Sp -----------------------
  g_se <- irs_grad_se(Se, Sp, prev, Se_ref, Sp_ref)
  g_sp <- irs_grad_sp(Se, Sp, prev, Se_ref, Sp_ref)
  NVar_se_corrected <- irs_delta_nvar(g_se, app$p)
  NVar_sp_corrected <- irs_delta_nvar(g_sp, app$p)

  # Note: n_se_corrected / n_sp_corrected are already TOTAL sample sizes
  # (N * Var is a total-N quantity, unlike Buderer's/apparent's diseased-
  # or ref-positive-subsample counts), so N_corrected is their max, with
  # no division by prevalence or P(ref+/-) afterwards.
  n_se_corrected <- ceiling(z_alpha^2 * NVar_se_corrected / d_se^2)
  n_sp_corrected <- ceiling(z_alpha^2 * NVar_sp_corrected / d_sp^2)
  N_corrected <- max(n_se_corrected, n_sp_corrected)
  N_corrected_loss <- ceiling(N_corrected / (1 - loss_rate))

  multiplier_se_corrected <- NVar_se_corrected / buderer_nvar_se
  multiplier_sp_corrected <- NVar_sp_corrected / buderer_nvar_sp

  # --- select the estimand-specific aliases ----------------------------
  if (estimand == "apparent") {
    n_diseased_adj <- n_refpos_apparent
    n_nondiseased_adj <- n_refneg_apparent
    N_adj <- N_apparent
    N_adj_loss <- N_apparent_loss
  } else {
    n_diseased_adj <- n_se_corrected
    n_nondiseased_adj <- n_sp_corrected
    N_adj <- N_corrected
    N_adj_loss <- N_corrected_loss
  }

  # --- bound the memory of the MC validation -------------------------
  # The validation allocates three B x N_adj matrices. Fail loudly, naming
  # the numbers, instead of letting R abort on a multi-gigabyte allocation.
  # For estimand = "corrected" this can trigger at moderate prevalence even
  # with a fine reference standard, since N_adjusted there is not bounded
  # the way the apparent estimand's is (see Details).
  if (B > 0 && B * N_adj > max_mc_cells) {
    stop(
      sprintf(
        paste0(
          "Monte Carlo validation would allocate B * N_adjusted = %.0f * %.0f ",
          "= %.3g cells\n  (> max_mc_cells = %.3g), i.e. several matrices of ",
          "that size.\n",
          "  N_adjusted = %.0f is the required N for estimand = \"%s\" ",
          "(Youden_ref = %.3f, prev = %.3f).\n",
          "  Set B = 0 to skip the validation, lower B, or raise ",
          "max_mc_cells if you have the memory."
        ),
        B, N_adj, B * N_adj, max_mc_cells, N_adj, estimand, youden_ref, prev
      ),
      call. = FALSE
    )
  }

  # Sensitivity table: N_apparent / N_corrected (and the informative VIF)
  # over a grid of Se_ref, Sp_ref, holding Se/Sp/prev fixed at the
  # planning values above.
  sens_table <- NULL
  if (sensitivity_table) {
    se_ref_grid <- seq(0.80, 1.00, 0.05)
    sp_ref_grid <- seq(0.85, 1.00, 0.05)
    grid <- expand.grid(Se_ref = se_ref_grid, Sp_ref = sp_ref_grid)
    grid <- grid[grid$Se_ref + grid$Sp_ref > 1, , drop = FALSE]
    grid$VIF <- 1 / (grid$Se_ref + grid$Sp_ref - 1)^2
    # inflation_factor is a same-value, more descriptive alias of VIF;
    # kept for backward compatibility. Neither multiplies N_apparent or
    # N_corrected below.
    grid$inflation_factor <- grid$VIF

    row_calc <- function(se_ref_i, sp_ref_i) {
      app_i <- irs_apparent(Se, Sp, prev, se_ref_i, sp_ref_i)
      n_rp_i <- ceiling(z_alpha^2 * app_i$Se_apparent * (1 - app_i$Se_apparent) / d_se^2)
      n_rn_i <- ceiling(z_alpha^2 * app_i$Sp_apparent * (1 - app_i$Sp_apparent) / d_sp^2)
      N_app_i <- ceiling(max(n_rp_i / app_i$P_ref_pos, n_rn_i / app_i$P_ref_neg))

      nvar_se_i <- irs_delta_nvar(irs_grad_se(Se, Sp, prev, se_ref_i, sp_ref_i), app_i$p)
      nvar_sp_i <- irs_delta_nvar(irs_grad_sp(Se, Sp, prev, se_ref_i, sp_ref_i), app_i$p)
      n_se_c_i <- ceiling(z_alpha^2 * nvar_se_i / d_se^2)
      n_sp_c_i <- ceiling(z_alpha^2 * nvar_sp_i / d_sp^2)
      N_corr_i <- max(n_se_c_i, n_sp_c_i)

      c(Se_apparent = app_i$Se_apparent, Sp_apparent = app_i$Sp_apparent,
        N_apparent = N_app_i, N_corrected = N_corr_i)
    }
    res_mat <- mapply(row_calc, grid$Se_ref, grid$Sp_ref)
    grid$Se_apparent <- res_mat["Se_apparent", ]
    grid$Sp_apparent <- res_mat["Sp_apparent", ]
    grid$N_apparent <- res_mat["N_apparent", ]
    grid$N_corrected <- res_mat["N_corrected", ]
    grid$N_adj <- if (estimand == "apparent") grid$N_apparent else grid$N_corrected

    sens_table <- grid
  }

  # MC validation
  mc_validation <- NULL
  if (B > 0) {
    set.seed(seed, kind = "Mersenne-Twister",
             normal.kind = "Inversion", sample.kind = "Rejection")

    # Generate true disease status
    true_disease <- stats::rbinom(B * N_adj, 1, prev)

    # For each subject, generate reference test result
    ref_result <- ifelse(
      true_disease == 1,
      stats::rbinom(length(true_disease), 1, Se_ref),
      1L - stats::rbinom(length(true_disease), 1, Sp_ref)
    )

    # Generate index test result
    index_result <- ifelse(
      true_disease == 1,
      stats::rbinom(length(true_disease), 1, Se),
      1L - stats::rbinom(length(true_disease), 1, Sp)
    )

    # Reshape into B datasets of size N_adj
    dim(true_disease) <- c(N_adj, B)
    dim(ref_result) <- c(N_adj, B)
    dim(index_result) <- c(N_adj, B)

    z_val <- z_alpha
    target_width_se <- 2 * d_se

    # For adjusted n: use all N_adj subjects
    # Apparent Se = P(index+ | ref+)
    ref_pos <- colSums(ref_result == 1)
    index_pos_given_ref_pos <- colSums(index_result == 1 & ref_result == 1)
    se_apparent_rep <- index_pos_given_ref_pos / pmax(ref_pos, 1)
    se_se_adj <- sqrt(se_apparent_rep * (1 - se_apparent_rep) / pmax(ref_pos, 1))
    width_adj <- 2 * z_val * se_se_adj
    P_width_adj <- mean(width_adj <= target_width_se)

    # Corrected Se_hat on the full N_adj table (diagnostic only: mean and
    # bias across replicates, not a per-replicate CI width check).
    n10_adj <- colSums(index_result == 1 & ref_result == 0)
    n01_adj <- colSums(index_result == 0 & ref_result == 1)
    p11c_adj <- index_pos_given_ref_pos / N_adj
    p10c_adj <- n10_adj / N_adj
    p01c_adj <- n01_adj / N_adj
    se_hat_adj <- (Sp_ref * p11c_adj - (1 - Sp_ref) * p10c_adj) /
      (p11c_adj + p01c_adj - (1 - Sp_ref))
    se_hat_adj[!is.finite(se_hat_adj)] <- NA
    se_corrected_mean_adj <- mean(se_hat_adj, na.rm = TRUE)

    # For unadjusted n: use only first N_unadj subjects
    n_use <- min(N_unadj, N_adj)
    idx_sub <- index_result[seq_len(n_use), , drop = FALSE]
    ref_sub <- ref_result[seq_len(n_use), , drop = FALSE]
    ref_pos_unadj <- colSums(ref_sub == 1)
    idx_pos_unadj <- colSums(idx_sub == 1 & ref_sub == 1)
    se_apparent_unadj_rep <- idx_pos_unadj / pmax(ref_pos_unadj, 1)
    se_se_unadj <- sqrt(se_apparent_unadj_rep * (1 - se_apparent_unadj_rep) /
                          pmax(ref_pos_unadj, 1))
    width_unadj <- 2 * z_val * se_se_unadj
    P_width_unadj <- mean(width_unadj <= target_width_se)

    n10_unadj <- colSums(idx_sub == 1 & ref_sub == 0)
    n01_unadj <- colSums(idx_sub == 0 & ref_sub == 1)
    p11c_u <- idx_pos_unadj / n_use
    p10c_u <- n10_unadj / n_use
    p01c_u <- n01_unadj / n_use
    se_hat_unadj <- (Sp_ref * p11c_u - (1 - Sp_ref) * p10c_u) /
      (p11c_u + p01c_u - (1 - Sp_ref))
    se_hat_unadj[!is.finite(se_hat_unadj)] <- NA
    se_corrected_mean_unadj <- mean(se_hat_unadj, na.rm = TRUE)

    # Mean of the APPARENT sensitivity P(index+ | ref+). With an imperfect
    # reference this converges to a value below the true Se; the gap is the
    # verification bias that "apparent" sizing does NOT remove (see Note).
    se_apparent_mean_adj <- mean(se_apparent_rep)
    se_apparent_mean_unadj <- mean(se_apparent_unadj_rep)

    mc_validation <- data.frame(
      scenario = c("unadjusted", "adjusted"),
      N = c(N_unadj, N_adj),
      P_width_target = c(P_width_unadj, P_width_adj),
      se_apparent = c(se_apparent_mean_unadj, se_apparent_mean_adj),
      se_true = c(Se, Se),
      bias = c(se_apparent_mean_unadj - Se, se_apparent_mean_adj - Se),
      se_corrected = c(se_corrected_mean_unadj, se_corrected_mean_adj),
      bias_corrected = c(se_corrected_mean_unadj - Se, se_corrected_mean_adj - Se),
      stringsAsFactors = FALSE
    )
  }

  results <- data.frame(
    estimand = c("apparent", "corrected"),
    N = c(N_apparent, N_corrected),
    N_with_loss = c(N_apparent_loss, N_corrected_loss),
    selected = c(estimand == "apparent", estimand == "corrected"),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      method = paste(
        "Sample Size for Imperfect Reference Standard",
        "(exact apparent-Se/Sp closed form and misclassification-",
        "correction delta method; Rogan-Gladen VIF retained as an",
        "informative prevalence quantity only)"
      ),
      estimand = estimand,
      n_diseased = n_diseased_adj,
      n_total = N_adj_loss,

      # classical Buderer baseline (reference treated as perfect)
      n_diseased_unadjusted = n_unadj_se,
      n_nondiseased_unadjusted = n_unadj_sp,
      N_unadjusted = N_unadj,
      N_buderer = N_unadj,

      # estimand-selected aliases (backward-compatible generic names)
      n_diseased_adjusted = n_diseased_adj,
      n_nondiseased_adjusted = n_nondiseased_adj,
      N_adjusted = N_adj,
      N_adjusted_loss = N_adj_loss,

      # estimand = "apparent": exact closed form
      Se_apparent = Se_apparent,
      Sp_apparent = Sp_apparent,
      P_ref_pos = P_ref_pos,
      P_ref_neg = P_ref_neg,
      n_refpos_apparent = n_refpos_apparent,
      n_refneg_apparent = n_refneg_apparent,
      N_apparent = N_apparent,
      N_apparent_loss = N_apparent_loss,

      # estimand = "corrected": delta-method variance of the exactly-
      # identified misclassification correction
      NVar_se_corrected = NVar_se_corrected,
      NVar_sp_corrected = NVar_sp_corrected,
      n_se_corrected = n_se_corrected,
      n_sp_corrected = n_sp_corrected,
      N_corrected = N_corrected,
      N_corrected_loss = N_corrected_loss,

      # informative only -- see Details/Note; does not size anything here
      youden_ref = youden_ref,
      VIF = VIF,
      Se_ref = Se_ref,
      Sp_ref = Sp_ref,
      multiplier_se_corrected = multiplier_se_corrected,
      multiplier_sp_corrected = multiplier_sp_corrected,
      multiplier_se_apparent = multiplier_se_apparent,
      multiplier_sp_apparent = multiplier_sp_apparent,

      sensitivity_table = sens_table,
      mc_validation = mc_validation,
      results = results,
      call = match.call()
    ),
    class = "dtasamplesize"
  )
}
