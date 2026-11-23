#' @title Prior-Data Conflict Check (Exact Beta-Binomial Tail Probabilities)
#' @description For observed 2x2 counts, computes the \strong{exact} tail
#'   probabilities of the Beta-Binomial prior predictive distribution of
#'   (i) the number of diseased (prevalence prior), (ii) the number of true
#'   positives among the diseased (Se prior) and (iii) the number of true
#'   negatives among the non-diseased (Sp prior). No simulation is used
#'   (log-scale \code{lchoose}/\code{lbeta} and direct summation).
#'
#' @param tp,fn,tn,fp Non-negative whole numbers: the observed counts.
#' @param prior_se,prior_sp,prior_prev Design priors, each \code{c(a, b)}
#'   (\code{a, b > 0}) or \code{NULL} to skip that quantity.
#' @param analysis_prior_se,analysis_prior_sp Analysis priors for Se and Sp.
#'   \code{NULL} (default) means \dQuote{equal to the design prior}. There is
#'   no analysis prior for prevalence.
#' @param priors Optionally an object of class \code{"dta_bam_priors"}
#'   (from \code{\link{bam_prior_from_pilot}}); then the other prior
#'   arguments must not be given (error otherwise). If \code{design} is
#'   not supplied it is taken from the object, that is, it is the design of
#'   the pilot. To check the data of a study planned with
#'   \code{\link{bam_sample_size}()}, which is a cohort, give
#'   \code{design = "cohort"} explicitly when the pilot was a case-control
#'   study.
#' @param design \code{"cohort"} (default) or \code{"case_control"}. In a
#'   case-control design the number of diseased is fixed by the design, so
#'   the prevalence row does not exist and a note says so. When
#'   \code{priors} is given and \code{design} is not, the design of the
#'   object (the design of the pilot) is used; for the data of a study
#'   planned with \code{\link{bam_sample_size}()} (a cohort) indicate
#'   \code{design = "cohort"}.
#' @param alpha Flagging threshold in (0, 1), default \code{0.05}.
#'
#' @details \strong{Quantities and references.} With \eqn{X \sim
#'   \mathrm{BetaBinomial}(n, a, b)}:
#'   \itemize{
#'     \item \code{prevalence}: \eqn{x = tp + fn} diseased out of
#'       \eqn{n = tp + fn + tn + fp}, against the prevalence prior (only
#'       in a cohort design). This prior acts in the design.
#'     \item Se: \eqn{x = tp} out of \eqn{n = tp + fn}.
#'     \item Sp: \eqn{x = tn} out of \eqn{n = tn + fp}.
#'   }
#'   For Se and Sp two references are reported: the predictive of the
#'   \strong{analysis prior} (\code{reference = "analysis"}: the inferential
#'   reference, as in sections 4.3 and 5.3 of Wilson et al. 2022) and that
#'   of the \strong{design prior} (\code{reference = "design"}: calibration
#'   of the assurance). When both priors coincide a single row is returned
#'   (\code{reference = "analysis"}, \code{same_as_design = TRUE}) and the
#'   output says so.
#'
#'   \strong{Tail probabilities and flagging rule.} \code{p_lower} =
#'   \eqn{P(X \le x_{obs})}, \code{p_upper} = \eqn{P(X \ge x_{obs})},
#'   \code{p_two_sided} = \eqn{\min(1, 2 \min(p_{lower}, p_{upper}))}. Each
#'   tail is flagged when \code{p <= alpha} and the two-sided value is
#'   flagged when \code{p_two_sided <= alpha} (a tolerance of 1e-12 absorbs
#'   floating-point error). Wilson et al. define the upper tail
#'   \eqn{\Pr(n \ge n_{obs})} in section 4.3 but report percentiles in
#'   section 5.3 that correspond to \eqn{P(X \le x_{obs})}; both tails are
#'   returned here.
#'
#'   \strong{Flat priors.} With Beta(1, 1) the predictive is uniform on
#'   \eqn{0,\dots,n} and the smallest possible tail is \eqn{1/(n+1)}: with
#'   the one-tail rule nothing can be flagged if \eqn{n \le 18} (for 10
#'   observations the smallest tail is 1/11 = 0.091), with the two-sided
#'   rule nothing if \eqn{n \le 38}, and for larger \eqn{n} only extreme
#'   results are flagged (30 of 30 gives a tail of 1/31 = 0.032), which
#'   reflects the flatness of the prior rather than a shift of the true
#'   value. The column \code{min_attainable_p} is the
#'   smallest two-sided p attainable for that \eqn{n} and prior, and
#'   \code{min_attainable_tail} the smallest attainable single tail. A
#'   warning is issued when \code{min_attainable_tail} exceeds \code{alpha}
#'   (nothing can be flagged) or when only \code{min_attainable_p} does
#'   (the two-sided rule cannot flag, a single tail can).
#'   \strong{Do not use the check with a flat initial prior.}
#'
#'   \strong{Intended use.} The reference use is \emph{after the study},
#'   with the study data. A secondary use at the planning stage is to
#'   compare the initial prior with the pilot, to review the initial prior.
#'   The check does \strong{not} inform the choice of the pilot weight
#'   \code{a0} in \code{\link{bam_prior_from_pilot}}: heterogeneity between
#'   pilot and study is not observed when planning, and a conflict does not
#'   say which source fails.
#'
#' @return A data frame of class \code{c("dta_bam_conflict", "data.frame")}
#'   with one row per quantity and reference and columns \code{quantity},
#'   \code{reference} (\code{"analysis"}, \code{"design"},
#'   \code{"prevalence"}), \code{same_as_design}, \code{x_obs}, \code{n},
#'   \code{a}, \code{b}, \code{p_lower}, \code{p_upper},
#'   \code{p_two_sided}, \code{flag_lower}, \code{flag_upper},
#'   \code{flag_two_sided}, \code{min_attainable_p} and
#'   \code{min_attainable_tail} (the smallest attainable one-sided tail).
#'   Attributes \code{alpha}, \code{design} and \code{notes}.
#'
#' @references
#' Wilson KJ, Williamson SF, Allen AJ, Williams CJ, Hellyer TP, Lendrem
#' BC (2022). Bayesian sample size determination for diagnostic accuracy
#' studies. Statistics in Medicine 41:2908-2922.
#' \doi{10.1002/sim.9393}
#'
#' @seealso \code{\link{bam_prior_from_pilot}}, \code{\link{bam_sample_size}}
#'
#' @examples
#' # VAP example (Wilson et al. 2022): 53 diseased of 150; 51 TP of 53
#' bam_prior_data_conflict(
#'   tp = 51, fn = 2, tn = 42, fp = 55,
#'   prior_se = c(25.9, 2.1), prior_sp = c(21, 36),
#'   prior_prev = c(29, 98))
#' @export
bam_prior_data_conflict <- function(tp, fn, tn, fp,
                                    prior_se = NULL, prior_sp = NULL,
                                    prior_prev = NULL,
                                    analysis_prior_se = NULL,
                                    analysis_prior_sp = NULL,
                                    priors = NULL,
                                    design = c("cohort", "case_control"),
                                    alpha = 0.05) {
  design_given <- !missing(design)
  design <- match.arg(design)
  .chk_count(tp, "tp"); .chk_count(fn, "fn")
  .chk_count(tn, "tn"); .chk_count(fp, "fp")
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1)
    stop("`alpha` must be a single number in (0, 1).", call. = FALSE)

  if (!is.null(priors)) {
    if (!inherits(priors, "dta_bam_priors"))
      stop("`priors` must be an object of class \"dta_bam_priors\" ",
           "(see bam_prior_from_pilot()).", call. = FALSE)
    if (!is.null(prior_se) || !is.null(prior_sp) || !is.null(prior_prev) ||
        !is.null(analysis_prior_se) || !is.null(analysis_prior_sp))
      stop("Do not combine `priors` with prior_se, prior_sp, prior_prev, ",
           "analysis_prior_se or analysis_prior_sp.", call. = FALSE)
    prior_se <- priors$prior_se
    prior_sp <- priors$prior_sp
    prior_prev <- priors$prior_prev
    analysis_prior_se <- priors$analysis_prior_se
    analysis_prior_sp <- priors$analysis_prior_sp
    if (!design_given) design <- priors$design
  }
  for (nm in c("prior_se", "prior_sp", "prior_prev",
               "analysis_prior_se", "analysis_prior_sp")) {
    v <- get(nm)
    if (!is.null(v)) .chk_beta_prior(v, nm)
  }
  if (is.null(prior_se) && is.null(prior_sp) && is.null(prior_prev) &&
      is.null(analysis_prior_se) && is.null(analysis_prior_sp))
    stop("Supply at least one prior (prior_se, prior_sp, prior_prev, ",
         "analysis_prior_se, analysis_prior_sp) or `priors`.", call. = FALSE)

  n_dis <- tp + fn
  n_non <- tn + fp
  n_all <- n_dis + n_non
  notes <- character(0)
  rows <- list()

  add_row <- function(quantity, reference, same, x, n, ab) {
    pm <- .bb_pmf(0:n, n, ab[1], ab[2])
    pl <- min(1, sum(pm[seq_len(x + 1L)]))
    pu <- min(1, sum(pm[(x + 1L):(n + 1L)]))
    p2 <- min(1, 2 * min(pl, pu))
    tol <- 1e-12
    min_tail <- min(pm[1L], pm[n + 1L])
    data.frame(
      quantity = quantity, reference = reference, same_as_design = same,
      x_obs = x, n = n, a = ab[1], b = ab[2],
      p_lower = pl, p_upper = pu, p_two_sided = p2,
      flag_lower = pl <= alpha + tol, flag_upper = pu <= alpha + tol,
      flag_two_sided = p2 <= alpha + tol,
      min_attainable_p = min(1, 2 * min_tail),
      min_attainable_tail = min_tail,
      stringsAsFactors = FALSE)
  }

  arm_rows <- function(quantity, x, n, d_prior, a_prior, arm_label) {
    if (n < 1) {
      notes <<- c(notes, paste0(arm_label, ": no observations (n = 0); ",
                                "row omitted."))
      return(NULL)
    }
    if (is.null(d_prior) && is.null(a_prior)) return(NULL)
    if (is.null(d_prior))
      return(list(add_row(quantity, "analysis", NA, x, n, a_prior)))
    if (is.null(a_prior) ||
        identical(as.numeric(a_prior), as.numeric(d_prior)))
      return(list(add_row(quantity, "analysis", TRUE, x, n, d_prior)))
    list(add_row(quantity, "analysis", FALSE, x, n, a_prior),
         add_row(quantity, "design", FALSE, x, n, d_prior))
  }

  # Prevalence (cohort only; acts in the design)
  if (design == "cohort") {
    if (!is.null(prior_prev) && n_all >= 1)
      rows <- c(rows, list(add_row("prevalence (diseased)", "prevalence",
                                   NA, n_dis, n_all, prior_prev)))
  } else if (!is.null(prior_prev)) {
    notes <- c(notes, paste0(
      "prevalence row omitted: in a case_control design the number of ",
      "diseased is fixed by the design and carries no information on ",
      "prevalence."))
  } else {
    notes <- c(notes, paste0(
      "prevalence row does not exist in a case_control design (the number ",
      "of diseased is fixed by the design)."))
  }
  rows <- c(rows, arm_rows("sensitivity (true positives)", tp, n_dis,
                           prior_se, analysis_prior_se, "sensitivity"))
  rows <- c(rows, arm_rows("specificity (true negatives)", tn, n_non,
                           prior_sp, analysis_prior_sp, "specificity"))

  if (length(rows) == 0L)
    stop("No check could be computed with the supplied priors and counts.",
         call. = FALSE)
  res <- do.call(rbind, rows)
  rownames(res) <- NULL

  # Warnings (single combined warning)
  msgs <- character(0)
  flat <- vapply(seq_len(nrow(res)),
                 function(i) isTRUE(all.equal(c(res$a[i], res$b[i]), c(1, 1))),
                 logical(1))
  if (any(flat))
    msgs <- c(msgs, paste0(
      "flat Beta(1, 1) prior in: ",
      paste(unique(paste0(res$quantity, " [", res$reference, "]")[flat]),
            collapse = "; "),
      ". The predictive is uniform and the smallest possible tail is ",
      "1/(n + 1); the check is not informative with a flat prior and only ",
      "extreme results can be flagged."))
  lab <- paste0(res$quantity, " [", res$reference, "], n = ", res$n)
  noflag_any <- res$min_attainable_tail > alpha + 1e-12
  noflag_two <- !noflag_any & res$min_attainable_p > alpha + 1e-12
  if (any(noflag_any))
    msgs <- c(msgs, paste0(
      "the check cannot flag anything, in either tail or two-sided ",
      "(min_attainable_tail > alpha = ", alpha, ") for: ",
      paste(unique(lab[noflag_any]), collapse = "; "), "."))
  if (any(noflag_two))
    msgs <- c(msgs, paste0(
      "the two-sided check cannot flag anything (min_attainable_p > alpha = ",
      alpha, "; a single tail can still be flagged) for: ",
      paste(unique(lab[noflag_two]), collapse = "; "), "."))
  if (length(msgs))
    warning(paste(msgs, collapse = " "), call. = FALSE)

  attr(res, "alpha") <- alpha
  attr(res, "design") <- design
  attr(res, "notes") <- notes
  class(res) <- c("dta_bam_conflict", "data.frame")
  res
}

#' Print Method for dta_bam_conflict Objects
#'
#' @param x An object of class \code{"dta_bam_conflict"}.
#' @param digits Significant digits for the probabilities (default 4).
#' @param ... Ignored.
#' @return Invisibly returns \code{x}.
#' @export
print.dta_bam_conflict <- function(x, digits = 4, ...) {
  alpha <- attr(x, "alpha")
  cat("\nPrior-data conflict check: exact Beta-Binomial prior predictive\n")
  cat("  design:", attr(x, "design"), "| alpha =", alpha, "\n")
  cat("  p_lower = P(X <= x_obs); p_upper = P(X >= x_obs); ",
      "p_two_sided = min(1, 2 * min(p_lower, p_upper))\n", sep = "")
  cat("  Flagging rule: each tail flagged if p <= alpha; two-sided flagged ",
      "if p_two_sided <= alpha.\n", sep = "")
  cat("  reference = 'analysis': inferential reference (analysis prior); ",
      "'design': calibration of the assurance;\n  'prevalence': prevalence ",
      "prior (acts in the design).\n\n", sep = "")
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  attr(d, "alpha") <- NULL; attr(d, "design") <- NULL
  attr(d, "notes") <- NULL
  show <- d[, c("quantity", "reference", "x_obs", "n", "a", "b",
                "p_lower", "p_upper", "p_two_sided", "flag_lower",
                "flag_upper", "flag_two_sided", "min_attainable_p"),
            drop = FALSE]
  for (cn in c("p_lower", "p_upper", "p_two_sided", "min_attainable_p"))
    show[[cn]] <- signif(show[[cn]], digits)
  print.data.frame(show, row.names = FALSE)
  same <- !is.na(d$same_as_design) & d$same_as_design
  if (any(same))
    cat("\n  Analysis and design priors coincide for: ",
        paste(d$quantity[same], collapse = "; "), " (one row each).\n",
        sep = "")
  for (nt in attr(x, "notes")) cat("  Note:", nt, "\n")
  cat("\n")
  invisible(x)
}

# Beta-Binomial pmf on the log scale (exact; no simulation)
.bb_pmf <- function(x, n, a, b) {
  exp(lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b))
}
