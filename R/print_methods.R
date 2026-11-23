#' Print Method for dtasamplesize Objects
#'
#' @param x An object of class \code{"dtasamplesize"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.dtasamplesize <- function(x, ...) {
  cat("\n", x$method, "\n\n")
  .print_bam_conventions(x)
  if (!is.null(x$status) && x$status != "converged") {
    cat("  status:", x$status,
        "-- NO N SATISFIES target_assurance; N_total/n_diseased are NA.\n")
    cat("  (see the warning issued by this call for diagnostic detail)\n")
  }
  # Is this a joint_sample_size() object? (R/joint_sample_size.R, since v0.6.6):
  # identified by its unique search_type value, since method is
  # a free-text header shared by every "dtasamplesize" object.
  is_joint <- !is.null(x$search_type) &&
    identical(x$search_type, "grid_first_candidate")

  # bam_sample_size() and joint_sample_size() (same contract, since v0.6.6):
  # target_reached = FALSE means no candidate N reached the target
  # probability, so N_total/n_total/joint_assurance (or
  # joint_prob_se_sp) are NA and must never be printed as if they were a
  # recommended design -- print the diagnostic fields instead.
  # With assurance_for = "se"/"sp" the field `joint_assurance` holds the
  # MARGINAL assurance of that arm: label it as such (joint mode: unchanged).
  ja_label <- switch(
    if (is.null(x$assurance_for)) "joint" else x$assurance_for,
    se = "marginal Se assurance",
    sp = "marginal Sp assurance",
    "joint_assurance"
  )
  marginal <- !identical(ja_label, "joint_assurance")
  if (!is.null(x$target_reached) && !isTRUE(x$target_reached)) {
    if (is_joint) {
      cat("  Target NOT reached -- no sample size returned.\n")
      cat("  (n_total / joint_prob_se_sp are NA; see the warning issued by",
          "this call)\n")
    } else {
      cat("  Target assurance NOT reached -- no sample size returned.\n")
      cat("  (N_total / n_total /", paste0(ja_label, " are NA; see the warning"),
          "issued by this call)\n")
    }
    if (!is.null(x$max_assurance_evaluated))
      cat(if (marginal) paste0("  max_assurance_evaluated (", ja_label, "):")
          else "  max_assurance_evaluated:", x$max_assurance_evaluated, "\n")
    if (!is.null(x$N_at_max_assurance))
      cat("  N_at_max_assurance:", x$N_at_max_assurance, "\n")
    if (!is.null(x$max_joint_prob_evaluated))
      cat("  max_joint_prob_evaluated:", x$max_joint_prob_evaluated, "\n")
    if (!is.null(x$N_at_max_joint_prob))
      cat("  N_at_max_joint_prob:", x$N_at_max_joint_prob, "\n")
    if (is_joint && !is.null(x$auc_gate) &&
          !isTRUE(x$auc_gate_passed))
      cat("  AUC gate: never passed anywhere in N_range (delta_auc too",
          "tight, or N_range too small); see auc_gate.\n")
    cat("  n_diseased:", x$n_diseased, "\n")
    if (!is.null(x$n_non_diseased))
      cat("  n_non_diseased:", x$n_non_diseased, "\n")
    if (!is.null(x$results)) {
      cat("\n")
      print(x$results)
    }
    return(invisible(x))
  }
  cat("  n_diseased:", x$n_diseased, "\n")
  if (!is.null(x$n_non_diseased))
    cat("  n_non_diseased:", x$n_non_diseased, "\n")
  cat("  N_total:", x$n_total, "\n")
  if (!is.null(x$joint_assurance)) {
    cat(paste0("  ", ja_label, ":"), x$joint_assurance, "\n")
  }
  if (is_joint) {
    if (!is.null(x$joint_prob_se_sp))
      cat("  joint_prob_se_sp:", x$joint_prob_se_sp, "\n")
    if (!is.null(x$N_range_used) && length(x$N_range_used) >= 2) {
      step <- diff(x$N_range_used)[1]
      cat("  n_total is the first candidate in N_range (grid step ",
          step, ") reaching target_prob; see ?joint_sample_size.\n",
          sep = "")
    } else {
      cat("  n_total is the first candidate in N_range reaching",
          "target_prob; see ?joint_sample_size.\n")
    }
    if (!is.null(x$joint_prob_mcse))
      cat("  joint_prob_mcse:", x$joint_prob_mcse, "\n")
    if (!is.null(x$auc_gate) && !is.null(x$auc_gate$first_N_auc_pass))
      cat("  AUC gate (Hanley-McNeil): first candidate in N_range passing =",
          x$auc_gate$first_N_auc_pass, "\n")
  }
  if (!is.null(x$N_buderer))
    cat("  N_buderer:", x$N_buderer, "\n")
  if (!is.null(x$N_imperfect))
    cat("  N_imperfect:", x$N_imperfect, "\n")
  if (!is.null(x$seed))
    cat("  seed:", x$seed, "\n")
  if (!is.null(x$B))
    cat("  B:", x$B, "\n")
  if (!is.null(x$results)) {
    cat("\n")
    print(x$results)
  }
  invisible(x)
}

# Conventions of a bam_sample_size() result (v0.7.0): width criterion, level,
# assurance objective and the analysis priors. Silent for every other
# "dtasamplesize" object (they have no `width_type` field).
.print_bam_conventions <- function(x) {
  if (is.null(x$width_type)) return(invisible(NULL))
  fmt_beta <- function(p) paste0("Beta(", format(p[1]), ", ", format(p[2]), ")")
  cat("  width criterion: ", .bam_width_label(x$width_type, x$alpha_ci),
      "; delta_se, delta_sp compared exactly (width <= delta)\n", sep = "")
  obj <- switch(
    if (is.null(x$assurance_for)) "joint" else x$assurance_for,
    joint = "joint (Se and Sp widths both met in the same study)",
    se = "marginal Se (only the Se width decides; prior_sp and delta_sp are not used for N_total)",
    sp = "marginal Sp (only the Sp width decides; prior_se and delta_se are not used for N_total)"
  )
  cat("  assurance objective: ", obj, "\n", sep = "")
  if (!is.null(x$analysis_prior_se)) {
    # Compare VALUES, not the "given" flag: an analysis prior that was passed
    # explicitly (or comes from a `priors` object with use_in = "both") but
    # equals the design prior is "same as design".
    same <- vapply(c("se", "sp"), function(arm) {
      identical(as.numeric(x[[paste0("analysis_prior_", arm)]]),
                as.numeric(x[[paste0("prior_", arm)]]))
    }, logical(1))
    if (!all(same)) {
      cat("  analysis priors differ from the design priors:\n")
      for (arm in c("se", "sp")) {
        ap <- x[[paste0("analysis_prior_", arm)]]
        dp <- x[[paste0("prior_", arm)]]
        cat("    ", toupper(substring(arm, 1, 1)), substring(arm, 2),
            ": analysis prior ", fmt_beta(ap), " (a + b = ", format(sum(ap)),
            ") vs design ", fmt_beta(dp), " (a + b = ", format(sum(dp)),
            ")", if (same[[arm]]) " [same as design]", "\n", sep = "")
      }
    } else {
      cat("  analysis priors: same as design priors (Se a + b = ",
          format(sum(x$analysis_prior_se)), ", Sp a + b = ",
          format(sum(x$analysis_prior_sp)), ")\n", sep = "")
    }
  }
  invisible(NULL)
}
