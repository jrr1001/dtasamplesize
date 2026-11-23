#' @title Design and Analysis Priors from Pilot Data (Power-Prior Update)
#' @description Updates Beta priors for sensitivity (Se), specificity (Sp)
#'   and, in a cohort design, prevalence with the counts of a previous
#'   stage (a pilot or earlier study), using a power prior with a
#'   \strong{fixed weight} \code{a0}. This is the fixed-weight
#'   power-prior option described in section 4.1 of Wilson et al. (2022);
#'   with \code{a0 = 1} it reproduces their Beta-Binomial updating
#'   (their VAP example, section 5.1). The result is an object of class
#'   \code{"dta_bam_priors"} holding the \strong{design priors} and the
#'   \strong{analysis priors}, to be passed to
#'   \code{\link{bam_sample_size}} (argument \code{priors}) or
#'   unpacked by hand.
#'
#' @param tp,fn,tn,fp Non-negative whole numbers: true positives, false
#'   negatives, true negatives and false positives observed in the
#'   previous stage.
#' @param base_prior_se,base_prior_sp,base_prior_prev Base (initial) Beta
#'   priors, each a length-2 numeric vector \code{c(a, b)} with
#'   \code{a, b > 0}, for sensitivity, specificity and prevalence, before
#'   incorporating the pilot.
#' @param a0 Power-prior weight in \eqn{(0, 1]}: the fraction of the pilot
#'   counts that is added to the base prior. Default \code{1} (the pilot
#'   counts enter in full; the updating of Wilson et al. 2022).
#' @param design \code{"cohort"} (default) or \code{"case_control"}: the
#'   \strong{sampling design of the pilot} (the previous stage), not of the
#'   study to be planned. The main study sized by
#'   \code{\link{bam_sample_size}} is always a cohort (the total \code{N} is
#'   fixed and the number of diseased subjects is random), so a pilot with
#'   \code{design = "case_control"} followed by a cohort main study is a valid
#'   combination, and the \code{"dta_bam_priors"} object can be passed as
#'   \code{priors} without any warning. In a
#'   cohort design the prevalence prior is also updated with the pilot
#'   (diseased = \code{tp + fn}, non-diseased = \code{tn + fp}). In a
#'   case-control design the group sizes are set by the design and carry no
#'   information on prevalence, so the prevalence prior is \strong{not}
#'   updated (as in Williamson et al. 2023): it stays equal to the base
#'   prior.
#' @param use_in \code{"both"} (default) or \code{"design"}. With
#'   \code{"both"} the updated prior is used as the design prior and as
#'   the analysis prior (as in the VAP example of Wilson et al., section
#'   5.1, and in Williamson et al.). With \code{"design"} the updated
#'   prior is the design prior and the \strong{base} prior is kept as
#'   the analysis prior.
#'
#' @details Updating (with \eqn{x_0} counts from the pilot and base prior
#'   Beta(\eqn{a}, \eqn{b})):
#'   \itemize{
#'     \item Se: Beta(\eqn{a + a_0\,tp}, \eqn{b + a_0\,fn});
#'     \item Sp: Beta(\eqn{a + a_0\,tn}, \eqn{b + a_0\,fp});
#'     \item prevalence (\code{design = "cohort"} only):
#'       Beta(\eqn{a + a_0\,(tp + fn)}, \eqn{b + a_0\,(tn + fp)}).
#'   }
#'
#'   \strong{Choosing a0.} The choice of \code{a0} rests on how similar the
#'   pilot population is to the population of the planned study, and on a
#'   sensitivity analysis over \code{a0}; it is \emph{not} based on a
#'   prior-data conflict check (\code{\link{bam_prior_data_conflict}}),
#'   which does not say how heterogeneous the pilot and the study are nor
#'   which source is at fault.
#'
#'   \strong{Not implemented.} The procedure under weaker assumptions in
#'   Supplement B.2 of Wilson et al. (2022) (a hierarchical prior fitted
#'   by MCMC) is not a power prior and is not implemented here; it is
#'   cited as an alternative.
#'
#'   \strong{Warning on manual use.} If you pass the updated prior by hand
#'   as \code{prior_se} (or \code{prior_sp}) to
#'   \code{\link{bam_sample_size}} and leave \code{analysis_prior_se}
#'   (\code{analysis_prior_sp}) \code{NULL}, that updated prior is used
#'   \emph{also in the analysis}. To keep the base prior as the analysis
#'   prior, use \code{use_in = "design"} and pass the whole object with
#'   \code{priors =}, or pass \code{analysis_prior_se} and
#'   \code{analysis_prior_sp} explicitly.
#'
#' @return An object of class \code{"dta_bam_priors"}, a list with
#'   components:
#'   \describe{
#'     \item{\code{prior_se}, \code{prior_sp}, \code{prior_prev}}{Design
#'       priors, each \code{c(a, b)}.}
#'     \item{\code{analysis_prior_se}, \code{analysis_prior_sp}}{Analysis
#'       priors (never \code{NULL}): equal to the design priors when
#'       \code{use_in = "both"}, equal to the base priors when
#'       \code{use_in = "design"}. There is no analysis prior for
#'       prevalence.}
#'     \item{\code{base}}{List with the base priors \code{prior_se},
#'       \code{prior_sp}, \code{prior_prev}.}
#'     \item{\code{pilot}}{List with the counts \code{tp}, \code{fn},
#'       \code{tn}, \code{fp}.}
#'     \item{\code{a0}, \code{design}, \code{use_in}}{The arguments used.}
#'   }
#'
#' @references
#' Wilson KJ, Williamson SF, Allen AJ, Williams CJ, Hellyer TP, Lendrem
#' BC (2022). Bayesian sample size determination for diagnostic accuracy
#' studies. Statistics in Medicine 41:2908-2922.
#' \doi{10.1002/sim.9393}
#'
#' Williamson SF, Williams CJ, Lendrem BC, Wilson KJ (2023). Sample size
#' determination for point-of-care COVID-19 diagnostic tests: a Bayesian
#' approach. \emph{Diagn Progn Res} 7:17. \doi{10.1186/s41512-023-00153-1}
#'
#' @seealso \code{\link{bam_sample_size}}, \code{\link{bam_prior_data_conflict}}
#'
#' @examples
#' # VAP example of Wilson et al. (2022), section 5.1 (a0 = 1, cohort)
#' p <- bam_prior_from_pilot(
#'   tp = 16, fn = 1, tn = 20, fp = 35,
#'   base_prior_se = c(9.9, 1.1), base_prior_sp = c(1, 1),
#'   base_prior_prev = c(12, 43))
#' p$prior_se    # c(25.9, 2.1)
#' p$prior_prev  # c(29, 98)
#' p$prior_sp    # c(21, 36)
#' p
#'
#' # Half weight to the pilot, updated prior only in the design
#' bam_prior_from_pilot(16, 1, 20, 35, c(9.9, 1.1), c(1, 1), c(12, 43),
#'                      a0 = 0.5, use_in = "design")
#' @export
bam_prior_from_pilot <- function(tp, fn, tn, fp,
                                 base_prior_se, base_prior_sp,
                                 base_prior_prev,
                                 a0 = 1,
                                 design = c("cohort", "case_control"),
                                 use_in = c("both", "design")) {
  design <- match.arg(design)
  use_in <- match.arg(use_in)

  .chk_count(tp, "tp"); .chk_count(fn, "fn")
  .chk_count(tn, "tn"); .chk_count(fp, "fp")
  .chk_beta_prior(base_prior_se, "base_prior_se")
  .chk_beta_prior(base_prior_sp, "base_prior_sp")
  .chk_beta_prior(base_prior_prev, "base_prior_prev")
  if (!is.numeric(a0) || length(a0) != 1L || !is.finite(a0) ||
      a0 <= 0 || a0 > 1)
    stop("`a0` must be a single number in (0, 1].", call. = FALSE)

  base_prior_se   <- as.numeric(base_prior_se)
  base_prior_sp   <- as.numeric(base_prior_sp)
  base_prior_prev <- as.numeric(base_prior_prev)

  new_se <- c(base_prior_se[1] + a0 * tp, base_prior_se[2] + a0 * fn)
  new_sp <- c(base_prior_sp[1] + a0 * tn, base_prior_sp[2] + a0 * fp)
  new_prev <- if (design == "cohort") {
    c(base_prior_prev[1] + a0 * (tp + fn),
      base_prior_prev[2] + a0 * (tn + fp))
  } else {
    base_prior_prev
  }

  if (use_in == "both") {
    an_se <- new_se
    an_sp <- new_sp
  } else {
    an_se <- base_prior_se
    an_sp <- base_prior_sp
  }

  structure(
    list(
      prior_se = new_se, prior_sp = new_sp, prior_prev = new_prev,
      analysis_prior_se = an_se, analysis_prior_sp = an_sp,
      base = list(prior_se = base_prior_se, prior_sp = base_prior_sp,
                  prior_prev = base_prior_prev),
      pilot = list(tp = tp, fn = fn, tn = tn, fp = fp),
      a0 = a0, design = design, use_in = use_in),
    class = "dta_bam_priors")
}

#' Print Method for dta_bam_priors Objects
#'
#' @param x An object of class \code{"dta_bam_priors"}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dta_bam_priors <- function(x, ...) {
  fmt <- function(ab) {
    sprintf("Beta(%s, %s)  mean = %.4f  a + b = %s",
            format(ab[1], digits = 10), format(ab[2], digits = 10),
            ab[1] / sum(ab), format(sum(ab), digits = 10))
  }
  cat("\nBeta priors from previous-stage (pilot) data -- power prior with",
      "fixed weight a0\n\n")
  cat("  design:", x$design, "| a0 =", format(x$a0), "| use_in =", x$use_in,
      "\n")
  cat("  pilot counts: tp =", x$pilot$tp, " fn =", x$pilot$fn,
      " tn =", x$pilot$tn, " fp =", x$pilot$fp, "\n\n")
  cat("  Base priors (before the pilot)\n")
  cat("    Se   ", fmt(x$base$prior_se), "\n")
  cat("    Sp   ", fmt(x$base$prior_sp), "\n")
  cat("    prev ", fmt(x$base$prior_prev), "\n\n")
  cat("  Design priors (updated with the pilot)\n")
  cat("    Se   ", fmt(x$prior_se), "\n")
  cat("    Sp   ", fmt(x$prior_sp), "\n")
  cat("    prev ", fmt(x$prior_prev),
      if (x$design == "case_control") "(not updated: case_control)" else "",
      "\n\n")
  cat("  Analysis priors",
      if (x$use_in == "both") "(= design priors)" else "(= base priors)",
      "\n")
  cat("    Se   ", fmt(x$analysis_prior_se), "\n")
  cat("    Sp   ", fmt(x$analysis_prior_sp), "\n")
  cat("    (no analysis prior for prevalence)\n\n")
  invisible(x)
}

# --- internal validators (shared with bam_prior_data_conflict) -------------

.chk_count <- function(v, nm) {
  if (!is.numeric(v) || length(v) != 1L || is.na(v) || !is.finite(v) ||
      v < 0 || abs(v - round(v)) > 1e-8)
    stop("`", nm, "` must be a single non-negative whole number.",
         call. = FALSE)
  invisible(TRUE)
}

.chk_beta_prior <- function(v, nm) {
  if (!is.numeric(v) || length(v) != 2L || anyNA(v) || any(!is.finite(v)) ||
      any(v <= 0))
    stop("`", nm, "` must be a numeric vector c(a, b) with a > 0 and b > 0.",
         call. = FALSE)
  invisible(TRUE)
}
