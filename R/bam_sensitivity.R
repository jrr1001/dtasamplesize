# Batch scenarios and local (Hellinger) prior sensitivity for
# bam_sample_size(method = "exact").  (v0.7.0)

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Scenario-level arguments accepted by bam_sensitivity().
.bam_sens_scenario_args <- c(
  "prior_se", "prior_sp", "prior_prev", "analysis_prior_se",
  "analysis_prior_sp", "priors", "delta_se", "delta_sp", "width_type",
  "assurance_for", "alpha_ci", "target_assurance", "N_max"
)

# Arguments that define the priors (design and analysis).
.bam_sens_prior_args <- c(
  "prior_se", "prior_sp", "prior_prev", "analysis_prior_se",
  "analysis_prior_sp", "priors"
)

# Internal defaults for the legacy per-arm Monte Carlo diagnostics of
# bam_sample_size(). Under method = "exact" the headline N_total and
# joint_assurance do not depend on B, n_range or seed, so these are set to
# the cheapest values that do not trigger the small-B advisory; the caller can
# override them through `...`.
.bam_sens_cheap_defaults <- list(B = 1000L, n_range = 20L)

# Warnings that only restate what the returned table already reports
# (no crossing found; legacy per-arm search exhausted).
.bam_sens_muffle <- function(w) {
  msg <- conditionMessage(w)
  grepl("^No n in n_range achieved target assurance", msg) ||
    grepl("^No integer N in 2:", msg) ||
    grepl("^No N in N_range achieved the target", msg)
}

# The warning of .bam_warn_prior_only(): the analysis prior of a used arm
# already meets its target width before any data are observed.
.bam_sens_is_prior_only <- function(w) {
  grepl("^The analysis prior for", conditionMessage(w))
}

# Run an expression, muffling the warnings above. If `sink` (an environment)
# is given, the prior-only warning of .bam_warn_prior_only() is also muffled
# and its message is recorded in `sink$msgs`, so that the caller can report
# it ONCE per scenario or contour point instead of once per internal call
# (the N - 1 and n_ref companion calls repeat it). The record travels in the
# value returned by each job, so it works the same in sequential and in
# parallel runs (the warnings of PSOCK workers never reach the main session).
.bam_sens_quiet <- function(expr, sink = NULL) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (.bam_sens_muffle(w)) invokeRestart("muffleWarning")
      if (!is.null(sink) && .bam_sens_is_prior_only(w)) {
        sink$msgs <- c(sink$msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    }
  )
}

# A sink for .bam_sens_quiet().
.bam_sens_sink <- function() {
  e <- new.env(parent = emptyenv())
  e$msgs <- character()
  e
}

# The (unique) prior-only messages recorded in a sink.
.bam_sens_sunk <- function(sink) unique(sink$msgs)

# lapply() or a PSOCK-cluster parLapply(); the cluster workers load the
# package with the same library path (or, for a development tree loaded with
# pkgload::load_all(), load that tree).
.bam_sens_lapply <- function(X, FUN, parallel = FALSE, cores = 2L) {
  if (!parallel || length(X) <= 1L) return(lapply(X, FUN))
  cores <- min(as.integer(cores), length(X))
  libs <- .libPaths()
  ns <- asNamespace("dtasamplesize")
  dev_path <- if (exists(".__DEVTOOLS__", envir = ns, inherits = FALSE)) {
    getNamespaceInfo(ns, "path")
  } else {
    NULL
  }
  init <- function(libs, dev_path) {
    .libPaths(libs)
    if (!is.null(dev_path)) {
      getExportedValue("pkgload", "load_all")(dev_path, quiet = TRUE)
    } else {
      loadNamespace("dtasamplesize")
    }
    invisible(TRUE)
  }
  # Must not carry the package namespace as its enclosure: it is shipped to
  # workers that have not loaded the package yet.
  environment(init) <- globalenv()
  cl <- parallel::makeCluster(cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, init, libs, dev_path)
  parallel::parLapply(cl, X, FUN)
}

.bam_sens_check_parallel <- function(parallel, cores) {
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    stop("`parallel` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
      cores < 1 || cores != round(cores)) {
    stop("`cores` must be a single positive whole number.", call. = FALSE)
  }
  invisible(TRUE)
}

# Names the user may pass through `...` to bam_sample_size().
.bam_sens_check_dots <- function(dots, reserved, where) {
  if (length(dots) == 0L) return(invisible(TRUE))
  nms <- names(dots)
  if (is.null(nms) || any(!nzchar(nms))) {
    stop("All arguments passed through `...` of ", where,
         " must be named.", call. = FALSE)
  }
  ok <- setdiff(names(formals(bam_sample_size)), "method")
  bad <- setdiff(nms, ok)
  if (length(bad)) {
    if ("method" %in% bad) {
      stop(where, " always uses method = \"exact\"; do not pass `method`.",
           call. = FALSE)
    }
    stop("Unknown argument(s) passed through `...` of ", where, ": ",
         paste(bad, collapse = ", "),
         ". Valid names are the arguments of bam_sample_size().",
         call. = FALSE)
  }
  if ("N_range" %in% nms && "N_range" %in% reserved) {
    stop(where, " scans every integer N from 2 and also evaluates N_total - 1 ",
         "(or n_ref), so `N_range` cannot be passed through `...`.",
         call. = FALSE)
  }
  clash <- intersect(nms, reserved)
  if (length(clash)) {
    stop("Argument(s) not allowed through `...` of ", where, ": ",
         paste(clash, collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

# --- data frame of scenarios -> list of scenarios -----------------------------

.bam_sens_pair_cols <- c("prior_se", "prior_sp", "prior_prev",
                         "analysis_prior_se", "analysis_prior_sp")
.bam_sens_scalar_cols <- c("delta_se", "delta_sp", "width_type",
                           "assurance_for", "alpha_ci", "target_assurance",
                           "N_max")

.bam_sens_df_to_list <- function(df) {
  nms <- names(df)
  pair_cols <- c(paste0(rep(.bam_sens_pair_cols, each = 2L),
                        c("_a", "_b")))
  known <- c("scenario", "name", pair_cols, .bam_sens_scalar_cols)
  bad <- setdiff(nms, known)
  if (length(bad)) {
    stop("Unknown column(s) in the scenario data frame: ",
         paste(bad, collapse = ", "), ". Allowed: scenario (or name), ",
         paste(pair_cols, collapse = ", "), ", ",
         paste(.bam_sens_scalar_cols, collapse = ", "), ".", call. = FALSE)
  }
  for (p in .bam_sens_pair_cols) {
    has_a <- paste0(p, "_a") %in% nms
    has_b <- paste0(p, "_b") %in% nms
    if (xor(has_a, has_b)) {
      stop("Columns `", p, "_a` and `", p, "_b` must be given together ",
           "(both or none).", call. = FALSE)
    }
  }
  if (all(c("scenario", "name") %in% nms)) {
    stop("Use either a `scenario` or a `name` column, not both.",
         call. = FALSE)
  }
  lab_col <- intersect(c("scenario", "name"), nms)
  n <- nrow(df)
  labels <- if (length(lab_col)) {
    as.character(df[[lab_col]])
  } else if (is.character(attr(df, "row.names")) ||
             .row_names_info(df) > 0L) {
    rownames(df)
  } else {
    NULL
  }
  out <- vector("list", n)
  for (i in seq_len(n)) {
    sc <- list()
    for (p in .bam_sens_pair_cols) {
      ca <- paste0(p, "_a")
      cb <- paste0(p, "_b")
      if (!(ca %in% nms)) next
      va <- df[[ca]][i]
      vb <- df[[cb]][i]
      if (is.na(va) && is.na(vb)) next
      if (xor(is.na(va), is.na(vb))) {
        stop("Row ", i, ": `", ca, "` and `", cb, "` must both be given or ",
             "both be NA.", call. = FALSE)
      }
      sc[[p]] <- c(as.numeric(va), as.numeric(vb))
    }
    for (s in .bam_sens_scalar_cols) {
      if (!(s %in% nms)) next
      v <- df[[s]][i]
      if (is.na(v)) next
      sc[[s]] <- if (is.factor(v)) as.character(v) else v
    }
    out[[i]] <- sc
  }
  if (!is.null(labels)) names(out) <- labels
  out
}

.bam_sens_normalise_scenarios <- function(scenarios) {
  if (is.data.frame(scenarios)) {
    if (nrow(scenarios) < 1L) {
      stop("`scenarios` must contain at least one scenario.", call. = FALSE)
    }
    scenarios <- .bam_sens_df_to_list(scenarios)
  } else if (!is.list(scenarios)) {
    stop("`scenarios` must be a list of scenarios or a data frame.",
         call. = FALSE)
  }
  if (length(scenarios) < 1L) {
    stop("`scenarios` must contain at least one scenario.", call. = FALSE)
  }
  ok <- vapply(scenarios, function(s) is.list(s) && !inherits(s, "dta_bam_priors"),
               logical(1))
  if (!all(ok)) {
    stop("Each scenario must be a list of arguments (a list of lists, or a ",
         "data frame with one row per scenario).", call. = FALSE)
  }
  nms <- names(scenarios)
  if (is.null(nms)) {
    nms <- paste0("S", seq_along(scenarios))
  } else if (any(is.na(nms) | !nzchar(nms))) {
    stop("Either all scenarios or none must have names (empty names found).",
         call. = FALSE)
  }
  if (anyDuplicated(nms)) {
    stop("Scenario names must be unique; duplicated: ",
         paste(unique(nms[duplicated(nms)]), collapse = ", "), ".",
         call. = FALSE)
  }
  names(scenarios) <- nms
  for (nm in nms) {
    sc <- scenarios[[nm]]
    if (length(sc) && (is.null(names(sc)) || any(!nzchar(names(sc))))) {
      stop("Scenario '", nm, "': all entries must be named.", call. = FALSE)
    }
    bad <- setdiff(names(sc), .bam_sens_scenario_args)
    if (length(bad)) {
      stop("Scenario '", nm, "': unknown argument(s): ",
           paste(bad, collapse = ", "), ". A scenario may set: ",
           paste(.bam_sens_scenario_args, collapse = ", "), ".",
           call. = FALSE)
    }
    if (anyDuplicated(names(sc))) {
      stop("Scenario '", nm, "': repeated argument name(s).", call. = FALSE)
    }
    if ("priors" %in% names(sc) &&
        any(setdiff(.bam_sens_prior_args, "priors") %in% names(sc))) {
      stop("Scenario '", nm, "': `priors` supplies prior_se, prior_sp, ",
           "prior_prev, analysis_prior_se and analysis_prior_sp; do not also ",
           "set them in the same scenario.", call. = FALSE)
    }
  }
  scenarios
}

# One bam_sample_size(method = "exact") call plus its N - 1 companion.
.bam_sens_one <- function(job) {
  args <- job$args
  sink <- .bam_sens_sink()
  res <- tryCatch(
    .bam_sens_quiet(do.call(bam_sample_size,
                            c(args, list(method = "exact"))), sink),
    error = function(e) e
  )
  if (inherits(res, "error")) return(list(error = conditionMessage(res)))
  A_nm1 <- NA_real_
  if (isTRUE(res$target_reached) && res$N_total >= 3L) {
    args_m1 <- args
    args_m1$N_range <- res$N_total - 1L
    r1 <- tryCatch(
      .bam_sens_quiet(do.call(bam_sample_size,
                              c(args_m1, list(method = "exact"))), sink),
      error = function(e) e
    )
    if (inherits(r1, "error")) return(list(error = conditionMessage(r1)))
    A_nm1 <- r1$max_assurance_evaluated
  }
  list(
    N_total = if (isTRUE(res$target_reached)) as.integer(res$N_total)
              else NA_integer_,
    A_N = if (isTRUE(res$target_reached)) res$joint_assurance else NA_real_,
    A_N_minus_1 = A_nm1,
    target_reached = isTRUE(res$target_reached),
    max_assurance_evaluated = res$max_assurance_evaluated,
    N_at_max_assurance = as.integer(res$N_at_max_assurance),
    width_type = res$width_type, assurance_for = res$assurance_for,
    alpha_ci = res$alpha_ci, target_assurance = res$target_assurance,
    delta_se = res$delta_se, delta_sp = res$delta_sp,
    prior_se = res$prior_se, prior_sp = res$prior_sp,
    prior_prev = res$prior_prev,
    analysis_prior_se = res$analysis_prior_se,
    analysis_prior_sp = res$analysis_prior_sp,
    prior_only = .bam_sens_sunk(sink)
  )
}

# The single warning of bam_sensitivity() for the scenarios whose analysis
# prior already meets the width target with no data.
.bam_sens_warn_prior_only <- function(names_affected) {
  warning(
    "In ", length(names_affected), " scenario(s) (",
    paste0("'", names_affected, "'", collapse = ", "), ") the analysis prior ",
    "of a used arm already has a prior-only interval width within its target ",
    "delta: the width criterion is met before any data are observed, so ",
    "success is dominated by the analysis prior. Use a weaker or flat ",
    "analysis prior (analysis_prior_se, analysis_prior_sp) if this is not ",
    "intended. See bam_sample_size() for the detail of each call.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# bam_sensitivity()
# ---------------------------------------------------------------------------

#' Batch Scenarios for the Exact Bayesian Assurance Method
#'
#' Runs \code{\link{bam_sample_size}(method = "exact")} over a set of
#' scenarios (design and analysis priors, widths, criterion, objective,
#' credible level, target assurance) and returns one row per scenario with the
#' smallest total sample size \eqn{N} that meets the target, the assurance
#' \eqn{A(N)} reached there and the assurance \eqn{A(N - 1)} one subject
#' earlier. It is the batch counterpart of a single call and the tool behind
#' a prior or scenario sensitivity table.
#'
#' @section Conventions (inherited from \code{bam_sample_size()}):
#' Each scenario is evaluated under the exact Beta-Binomial calculation, so
#' every convention of \code{\link{bam_sample_size}} applies and is reported in
#' the output columns: the width criterion (\code{width_type}: total width of
#' the equal-tailed credible interval, \code{"full_width"}, or median minus
#' one-sided lower bound, \code{"lower_half_width"}), the credible level
#' (\code{alpha_ci}), the objective (\code{assurance_for}: joint over Se and
#' Sp, or marginal for Se or Sp alone) and the target assurance
#' (\code{target_assurance}). \eqn{A(N)} is the assurance under the
#' \emph{design} priors (which generate the data) with the interval computed
#' under the \emph{analysis} priors; when \code{analysis_prior_se} or
#' \code{analysis_prior_sp} is not given the analysis prior equals the design
#' prior.
#'
#' @section Scenario formats:
#' \describe{
#'   \item{List}{A list of scenarios, each a list of arguments. The list may
#'     be named (the names become the \code{scenario} column) or unnamed
#'     (\code{S1}, \code{S2}, ...); a mix is an error. A scenario may set any
#'     of \code{prior_se}, \code{prior_sp}, \code{prior_prev},
#'     \code{analysis_prior_se}, \code{analysis_prior_sp} (each
#'     \code{c(a, b)}), \code{priors} (an object from
#'     \code{\link{bam_prior_from_pilot}}, instead of the five prior
#'     arguments), \code{delta_se}, \code{delta_sp}, \code{width_type},
#'     \code{assurance_for}, \code{alpha_ci}, \code{target_assurance} and
#'     \code{N_max}.}
#'   \item{Data frame}{One row per scenario. Beta priors are given as two
#'     columns, \code{prior_se_a} and \code{prior_se_b}, and likewise
#'     \code{prior_sp_a}/\code{prior_sp_b}, \code{prior_prev_a}/
#'     \code{prior_prev_b}, \code{analysis_prior_se_a}/\code{_b} and
#'     \code{analysis_prior_sp_a}/\code{_b} (both columns of a pair must be
#'     present; in a row they are both \code{NA} or both given). Scalar
#'     columns: \code{delta_se}, \code{delta_sp}, \code{width_type},
#'     \code{assurance_for}, \code{alpha_ci}, \code{target_assurance},
#'     \code{N_max}. The scenario label is the \code{scenario} (or
#'     \code{name}) column if present, otherwise character row names,
#'     otherwise \code{S1}, \code{S2}, .... An \code{NA} cell means "not set
#'     by this scenario". A \code{priors} object can only be given in the list
#'     format. The data frame and the equivalent list give identical results.}
#' }
#' Whatever a scenario does not set is taken from \code{...} (arguments common
#' to all scenarios) and, if absent there too, from the default of
#' \code{\link{bam_sample_size}}. A scenario value overrides the common one.
#' If a scenario sets a \code{priors} object, the common prior arguments are
#' ignored for that scenario; if the common arguments contain a \code{priors}
#' object and a scenario sets individual prior arguments, an error is raised
#' by \code{bam_sample_size()} (a \code{priors} object supplies all five
#' priors).
#'
#' @section Internal defaults for the legacy diagnostics:
#' \code{bam_sample_size()} also computes legacy per-arm Monte Carlo fields
#' that this function does not use and that do not influence the exact
#' headline result (\code{N_total}, \code{joint_assurance}). To keep each call
#' cheap, \code{B = 1000} and \code{n_range = 20} are used unless you pass
#' \code{B} or \code{n_range} through \code{...}. The warnings about an
#' exhausted search (no crossing within \code{N_max}) are not re-emitted: they
#' are reported in the \code{target_reached} column.
#'
#' @section Warning about the analysis prior:
#' If the effective analysis prior of an arm used by the objective already has
#' a prior-only interval (no data) within its target width, the width
#' criterion is met before any data are observed and success is dominated by
#' the analysis prior (see \code{\link{bam_sample_size}}). The individual
#' calls do not emit that warning (the companion call at \eqn{N - 1} would
#' repeat it): the function emits a single warning that names the affected
#' scenarios, in the same way with \code{parallel = TRUE} and
#' \code{parallel = FALSE}. Use \code{\link{bam_sample_size}} on a scenario
#' for the detail.
#'
#' @param scenarios A list of scenarios or a data frame with one row per
#'   scenario (see \emph{Scenario formats}).
#' @param ... Arguments common to all scenarios, passed to
#'   \code{\link{bam_sample_size}} (for example \code{prior_prev},
#'   \code{delta_se}, \code{N_max}). \code{method} is always \code{"exact"} and
#'   cannot be passed; neither can \code{N_range}: every scenario scans all
#'   integers from 2 and the assurance at \eqn{N - 1} is computed separately,
#'   which a user-supplied \code{N_range} would invalidate (an error is
#'   raised).
#' @param parallel Logical. If \code{TRUE}, scenarios are distributed over a
#'   PSOCK cluster of \code{cores} processes (works on Windows). Workers load
#'   the package with the same library path as the calling session. The result
#'   does not depend on this argument. Default \code{FALSE}.
#' @param cores Number of worker processes when \code{parallel = TRUE}
#'   (capped at the number of scenarios). Default 2.
#'
#' @return A data frame with one row per scenario and columns:
#'   \describe{
#'     \item{scenario}{Scenario label.}
#'     \item{width_type, assurance_for, alpha_ci, target_assurance}{The
#'       conventions used (see above).}
#'     \item{delta_se, delta_sp}{Width targets used (total width under
#'       \code{"full_width"}, median minus one-sided bound under
#'       \code{"lower_half_width"}).}
#'     \item{prior_se_a, prior_se_b, prior_sp_a, prior_sp_b, prior_prev_a,
#'       prior_prev_b}{Design priors used.}
#'     \item{analysis_prior_se_a, analysis_prior_se_b, analysis_prior_sp_a,
#'       analysis_prior_sp_b}{Effective analysis priors.}
#'     \item{N_total}{Smallest total \eqn{N} (integer scan from 2) whose
#'       assurance reaches \code{target_assurance}; \code{NA} if none up to
#'       \code{N_max}.}
#'     \item{A_N}{Assurance \eqn{A(N)} at \code{N_total} (\code{NA} if not
#'       reached).}
#'     \item{A_N_minus_1}{Assurance at \code{N_total - 1}, obtained with
#'       \code{bam_sample_size(N_range = N_total - 1)}; it is below the target
#'       by construction. \code{NA} if not reached or \code{N_total} is 2.}
#'     \item{target_reached}{Logical.}
#'     \item{max_assurance_evaluated, N_at_max_assurance}{Highest assurance
#'       seen by the search and where; informative mainly when the target is
#'       not reached (then \code{N_total}, \code{A_N} and \code{A_N_minus_1}
#'       are \code{NA}).}
#'   }
#'
#' @seealso \code{\link{bam_sample_size}}, \code{\link{bam_sensitivity_local}}
#' @examples
#' # Three scenarios, with widths wide enough for the example to be fast and for
#' # no analysis prior to meet its width before any data are observed
#' scen <- list(
#'   informative = list(prior_se = c(17, 3), prior_sp = c(18, 2)),
#'   weaker      = list(prior_se = c(8.5, 1.5), prior_sp = c(9, 1)),
#'   se_only     = list(prior_se = c(17, 3), prior_sp = c(18, 2),
#'                      assurance_for = "se")
#' )
#' bam_sensitivity(scen, prior_prev = c(4, 16), delta_se = 0.22,
#'                 delta_sp = 0.22)
#'
#' # The same, as a data frame
#' df <- data.frame(
#'   scenario = c("informative", "weaker"),
#'   prior_se_a = c(17, 8.5), prior_se_b = c(3, 1.5),
#'   prior_sp_a = c(18, 9), prior_sp_b = c(2, 1)
#' )
#' bam_sensitivity(df, prior_prev = c(4, 16), delta_se = 0.22,
#'                 delta_sp = 0.22)
#' @export
bam_sensitivity <- function(scenarios, ..., parallel = FALSE, cores = 2L) {
  .bam_sens_check_parallel(parallel, cores)
  dots <- list(...)
  .bam_sens_check_dots(dots, reserved = "N_range", where = "bam_sensitivity()")
  scenarios <- .bam_sens_normalise_scenarios(scenarios)
  nms <- names(scenarios)

  common <- c(.bam_sens_cheap_defaults[setdiff(names(.bam_sens_cheap_defaults),
                                               names(dots))], dots)
  jobs <- lapply(nms, function(nm) {
    sc <- scenarios[[nm]]
    args <- common
    if ("priors" %in% names(sc)) {
      args[intersect(.bam_sens_prior_args, names(args))] <- NULL
    } else if (any(setdiff(.bam_sens_prior_args, "priors") %in% names(sc)) &&
               "priors" %in% names(args)) {
      stop("Scenario '", nm, "' sets individual prior arguments but the ",
           "common arguments contain a `priors` object; a `priors` object ",
           "supplies all five priors. Put the priors in the scenarios or ",
           "drop the common `priors`.", call. = FALSE)
    }
    args[names(sc)] <- sc
    list(name = nm, args = args)
  })

  res <- .bam_sens_lapply(jobs, .bam_sens_one, parallel = parallel,
                          cores = cores)

  err <- vapply(res, function(r) !is.null(r$error), logical(1))
  # One summary warning for the scenarios whose analysis prior already meets
  # the width with no data (the internal calls do not emit it). It is emitted
  # before the error check so that it is not lost when a scenario fails.
  affected <- vapply(res, function(r) length(r$prior_only) > 0L, logical(1))
  if (any(affected)) .bam_sens_warn_prior_only(nms[affected])
  if (any(err)) {
    i <- which(err)[1L]
    stop("Scenario '", nms[i], "': ", res[[i]]$error, call. = FALSE)
  }

  col <- function(f) vapply(res, f, numeric(1))
  out <- data.frame(
    scenario = nms,
    width_type = vapply(res, function(r) r$width_type, character(1)),
    assurance_for = vapply(res, function(r) r$assurance_for, character(1)),
    alpha_ci = col(function(r) r$alpha_ci),
    target_assurance = col(function(r) r$target_assurance),
    delta_se = col(function(r) r$delta_se),
    delta_sp = col(function(r) r$delta_sp),
    prior_se_a = col(function(r) r$prior_se[1]),
    prior_se_b = col(function(r) r$prior_se[2]),
    prior_sp_a = col(function(r) r$prior_sp[1]),
    prior_sp_b = col(function(r) r$prior_sp[2]),
    prior_prev_a = col(function(r) r$prior_prev[1]),
    prior_prev_b = col(function(r) r$prior_prev[2]),
    analysis_prior_se_a = col(function(r) r$analysis_prior_se[1]),
    analysis_prior_se_b = col(function(r) r$analysis_prior_se[2]),
    analysis_prior_sp_a = col(function(r) r$analysis_prior_sp[1]),
    analysis_prior_sp_b = col(function(r) r$analysis_prior_sp[2]),
    N_total = vapply(res, function(r) r$N_total, integer(1)),
    A_N = col(function(r) r$A_N),
    A_N_minus_1 = col(function(r) r$A_N_minus_1),
    target_reached = vapply(res, function(r) r$target_reached, logical(1)),
    max_assurance_evaluated = col(function(r) r$max_assurance_evaluated),
    N_at_max_assurance = vapply(res, function(r) r$N_at_max_assurance,
                                integer(1)),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Hellinger contour of a Beta prior
# ---------------------------------------------------------------------------

# H^2 = 1 - BC, with the Bhattacharyya coefficient
# BC = B((a + a0)/2, (b + b0)/2) / sqrt(B(a, b) B(a0, b0)),
# computed on the log scale and with expm1() for accuracy near 0.
.bam_hellinger_d2 <- function(a, b, a0, b0) {
  lbc <- lbeta((a + a0) / 2, (b + b0) / 2) -
    0.5 * (lbeta(a, b) + lbeta(a0, b0))
  d2 <- -expm1(lbc)
  pmax(d2, 0)
}

.bam_dist <- function(a, b, a0, b0, distance) {
  d2 <- .bam_hellinger_d2(a, b, a0, b0)
  if (distance == "H") sqrt(d2) else d2
}

# Radius r(phi) such that dist((a0 + r cos, b0 + r sin), (a0, b0)) = epsilon.
# Returns NA (with the reason as an attribute) when no root is found.
.bam_contour_radius <- function(phi, a0, b0, epsilon, distance) {
  cs <- cos(phi)
  sn <- sin(phi)
  f <- function(r) .bam_dist(a0 + r * cs, b0 + r * sn, a0, b0, distance) - epsilon

  # a > 0 and b > 0 bound r in the directions with cos < 0 or sin < 0
  r_cap <- Inf
  if (cs < 0) r_cap <- min(r_cap, a0 / (-cs))
  if (sn < 0) r_cap <- min(r_cap, b0 / (-sn))
  bounded <- is.finite(r_cap)
  r_top <- if (bounded) r_cap * (1 - 1e-9) else 1e12 * (a0 + b0)

  # bracket by doubling: f(0) = -epsilon < 0
  lo <- 0
  hi <- min(sqrt(epsilon) * (a0 + b0), r_top)
  if (!(hi > 0)) hi <- r_top
  repeat {
    fh <- f(hi)
    if (is.finite(fh) && fh > 0) break
    if (hi >= r_top) return(NA_real_)
    lo <- hi
    hi <- min(2 * hi, r_top)
  }
  stats::uniroot(f, lower = lo, upper = hi, tol = 1e-15 * (1 + hi),
                 maxiter = 1000L)$root
}

# Relative tolerance of the distance reached by a contour point: the point is
# accepted when |d - epsilon| <= .bam_contour_tol * epsilon.
.bam_contour_tol <- 1e-6

# Estimated rounding error of the distance at a contour point, relative to the
# distance itself. The distance is 1 - BC = -expm1(log BC), with
# log BC = lbeta(am, bm) - (lbeta(a, b) + lbeta(a0, b0)) / 2, am = (a + a0) / 2,
# bm = (b + b0) / 2. The three terms are added (and subtracted) without
# cancellation of their errors. Each lbeta() carries an absolute rounding
# error of a few .Machine$double.eps that does NOT vanish when |lbeta| is
# small: for arguments below 10, lbeta() is evaluated as the logarithm of a
# quotient of gamma functions, so its absolute error is a few eps even when
# the result is close to 0 (for example Beta(1, 1)). Hence an absolute floor
# of 4 per term, and the absolute error of log BC, and of H^2 = 1 - BC when BC
# is close to 1, is estimated by
# err = 4 * eps * ((4 + |l1|) + (4 + |l2|) / 2 + (4 + |l3|) / 2),
# with l1 = lbeta(am, bm), l2 = lbeta(a, b), l3 = lbeta(a0, b0). In large
# priors the floor is negligible. The factor 4 is a safety factor. The
# estimate was checked against a reference without cancellation (Taylor series
# of log BC in the differences of the parameters): 12 priors x 5 values of
# epsilon x 72 directions, ratio actual error / bound at most 0.50. It is an
# estimate, not a guarantee. Relative to H^2 the bound is err / H^2; the
# distance "H" = sqrt(H^2) has a relative error of half of that, err / (2 H^2).
.bam_contour_precision <- function(a, b, a0, b0, d, distance) {
  err <- 4 * .Machine$double.eps *
    ((4 + abs(lbeta((a + a0) / 2, (b + b0) / 2))) +
       (4 + abs(lbeta(a, b))) / 2 + (4 + abs(lbeta(a0, b0))) / 2)
  d2 <- if (distance == "H") d^2 else d
  rel <- err / d2
  if (distance == "H") rel <- rel / 2
  rel
}

.bam_check_epsilon <- function(epsilon) {
  if (missing(epsilon) || is.null(epsilon)) {
    stop("`epsilon` has no default: give the distance to the initial prior ",
         "(Wilson et al. 2022 use 0.00354 with distance = \"H2\").",
         call. = FALSE)
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon <= 0 || epsilon >= 1) {
    stop("`epsilon` must be a single number in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

.bam_check_grid <- function(n_phi, phi) {
  if (!is.null(phi)) {
    if (!is.numeric(phi) || length(phi) < 1L || anyNA(phi) ||
        any(!is.finite(phi))) {
      stop("`phi` must be a numeric vector of finite angles (radians).",
           call. = FALSE)
    }
    return(as.numeric(phi))
  }
  if (!is.numeric(n_phi) || length(n_phi) != 1L || !is.finite(n_phi) ||
      n_phi != round(n_phi) || n_phi < 3) {
    stop("`n_phi` must be a single whole number >= 3 (number of directions).",
         call. = FALSE)
  }
  n_phi <- as.integer(n_phi)
  -pi + 2 * pi * (0:(n_phi - 1L)) / n_phi
}

#' Contour of Beta Priors at a Given Hellinger-Type Distance
#'
#' Returns the points \eqn{(a, b)} of the contour
#' \eqn{\{(a, b): d(\mathrm{Beta}(a, b), \mathrm{Beta}(a_0, b_0)) = \epsilon\}}
#' around an initial Beta prior, as used for the local prior sensitivity
#' analysis of Wilson et al. (2022, section 4.2). The points are parametrised
#' in polar coordinates in the hyperparameter plane,
#' \eqn{a = a_0 + r\cos\varphi}, \eqn{b = b_0 + r\sin\varphi}, and for each
#' direction \eqn{\varphi} the radius \eqn{r(\varphi)} is solved numerically
#' (\code{\link[stats]{uniroot}}) so that the distance equals \code{epsilon}
#' to numerical precision: the result is a \emph{contour}, not a ball.
#' The precision is checked at every point (see \emph{Precision}).
#'
#' @section Distance:
#' \code{distance = "H2"} (default) is \eqn{1 - BC}, where
#' \deqn{BC = \frac{B((a + a_0)/2, (b + b_0)/2)}{\sqrt{B(a, b)\,B(a_0, b_0)}}}
#' is the Bhattacharyya coefficient; this is the expression printed by Wilson
#' et al. (2022) under the name "Hellinger distance" and used in their code
#' (it is the \emph{squared} Hellinger distance). \code{distance = "H"} is
#' \eqn{\sqrt{1 - BC}}, the Hellinger distance proper. The two differ only by
#' the scale of \code{epsilon}: the \code{"H"} contour at \code{epsilon} is the
#' \code{"H2"} contour at \code{epsilon^2}.
#'
#' @section About the value 0.00354:
#' Wilson et al. (2022) use \code{epsilon = 0.00354} with the \code{"H2"}
#' expression and describe it as equivalent to a shift of 0.1 in the mean of a
#' standard normal random variable. For \eqn{N(0, 1)} against \eqn{N(\delta,
#' 1)}, \eqn{1 - BC = 1 - \exp(-\delta^2/8)}; with \eqn{\delta = 0.1} this is
#' 0.00125 (and \eqn{\sqrt{1 - BC} = 0.0353}), not 0.00354. \code{epsilon}
#' therefore has no default here; choose it explicitly and state which
#' distance it refers to. With \code{"H2"} and \code{epsilon = 0.00354} the
#' contour of Se \eqn{\mathrm{Beta}(9.9, 1.1)} spans \eqn{a \in [7.63, 12.80]}
#' and \eqn{b \in [0.892, 1.362]}, which is the contour of their Figure 1.
#'
#' @section Angular grid and failures:
#' By default the directions are \code{n_phi = 72} equally spaced angles in
#' \eqn{[-\pi, \pi)}, \eqn{\varphi_k = -\pi + 2\pi k/72}, \eqn{k = 0, \dots,
#' 71} (reproducible, no randomness). Give another number with \code{n_phi}
#' or an explicit vector of angles (radians) with \code{phi}, which overrides
#' \code{n_phi}. In the directions with \eqn{\cos\varphi < 0} or
#' \eqn{\sin\varphi < 0} the radius is bounded by \eqn{a > 0} and
#' \eqn{b > 0}; the root is searched inside that interval (up to
#' \eqn{r_{max}(1 - 10^{-9})}). The distance tends to 1 when \eqn{a} or
#' \eqn{b} tends to 0 and when \eqn{r \to \infty}, so a root exists for any
#' \code{epsilon} in (0, 1) except those numerically indistinguishable from 1;
#' if no root is found in a direction the point is \code{NA} and a warning
#' names the directions.
#'
#' @section Precision:
#' The radius is solved to a much finer tolerance than the one that matters,
#' so in ordinary cases (for example the contour of Se
#' \eqn{\mathrm{Beta}(9.9, 1.1)} at \code{epsilon = 0.00354}) the distance at
#' the returned point equals \code{epsilon} to about \eqn{10^{-12}} in
#' relative terms. Precision is lost when the initial prior is very
#' concentrated and \code{epsilon} very small: the Bhattacharyya coefficient
#' is then 1 to within rounding error, the distance computed as \eqn{1 - BC}
#' is dominated by cancellation (and, with \code{distance = "H"}, the square
#' root of a quantity dominated by rounding error), and it can
#' even be exactly 0.
#'
#' Two checks are applied at every returned point, and a point is flagged
#' (imprecise) if either fails.
#' \enumerate{
#'   \item \emph{Estimated rounding error.} The rounding error of
#'     \eqn{\log BC = \log B((a + a_0)/2, (b + b_0)/2) - (\log B(a, b) +
#'     \log B(a_0, b_0))/2} is estimated as
#'     \eqn{4\,\mathrm{eps}\,((4 + |\log B((a + a_0)/2, (b + b_0)/2)|) +
#'     (4 + |\log B(a, b)|)/2 + (4 + |\log B(a_0, b_0)|)/2)}, with eps the
#'     machine epsilon (\code{.Machine$double.eps}); this is the absolute
#'     error of \eqn{1 - BC} when \eqn{BC} is close to 1. The constant 4 in
#'     each term is an absolute floor: for arguments below 10 the logarithm of
#'     the beta function is computed as the logarithm of a quotient of gamma
#'     functions, with an absolute error of a few eps even when its value is
#'     close to 0 (for example Beta(1, 1)); in large priors the floor is
#'     negligible. Divided by \eqn{H^2 = 1 - BC} of the point (by
#'     \eqn{2H^2} with \code{distance = "H"}) it gives the relative error
#'     bound, returned in \code{attr(, "precision_bound")} (one element per
#'     direction; \code{NA} where there is no point). The point is flagged
#'     when this bound exceeds \eqn{10^{-6}}. The bound is a deterministic
#'     function of the prior, \code{epsilon} and the direction through the
#'     point, not of how the rounding happens to fall. It is a conservative
#'     estimate of the rounding error (it contains a safety factor of 4),
#'     \strong{not a guarantee}: it flags points that are in fact accurate and
#'     cannot prove that an unflagged point is exact. It was checked against a
#'     reference without cancellation (Taylor series of \eqn{\log BC} in the
#'     differences of the parameters) for 12 priors, 5 values of
#'     \code{epsilon} and 72 directions: the ratio of the actual error to the
#'     bound did not exceed 0.50.
#'   \item \emph{Recomputation.} The distance is recomputed at the returned
#'     point, which is flagged if it is 0, not finite or differs from
#'     \code{epsilon} by more than a relative \eqn{10^{-6}}. This check uses
#'     the same formula that loses precision, so it alone cannot be relied on
#'     (it can accept points whose error is larger than the tolerance because
#'     the rounding happens to cancel).
#' }
#' Examples that are flagged at every direction: \code{prior = c(1e5, 1e4)},
#' \code{epsilon = 1e-6}, \code{distance = "H"} (the recomputed distance is 0),
#' and \code{prior = c(1000, 100)}, \code{epsilon = 1e-8},
#' \code{distance = "H2"} (the bound is about \eqn{6 \times 10^{-5}}). When
#' any point is flagged a warning says that the numerical precision is not
#' sufficient for that prior and that \code{epsilon}, and suggests a larger
#' \code{epsilon} or the other distance. The points are returned unchanged,
#' flagged in the logical vector \code{attr(, "imprecise")} (one element per
#' direction).
#'
#' @param prior Numeric vector \code{c(a0, b0)}, the initial Beta prior
#'   (both positive).
#' @param epsilon Distance to the initial prior, in (0, 1). No default.
#' @param n_phi Number of equally spaced directions in \eqn{[-\pi, \pi)}
#'   (whole number >= 3). Default 72.
#' @param phi Optional numeric vector of angles in radians; overrides
#'   \code{n_phi}.
#' @param distance \code{"H2"} (default; \eqn{1 - BC}) or \code{"H"}
#'   (\eqn{\sqrt{1 - BC}}).
#'
#' @return A data frame with columns \code{phi}, \code{r}, \code{a}, \code{b}
#'   and \code{distance} (the distance recomputed at the returned point),
#'   with attributes \code{prior}, \code{epsilon}, \code{distance_type},
#'   \code{imprecise} (logical, one element per direction: \code{TRUE} where
#'   the estimated rounding error is too large, or the recomputed distance
#'   differs from \code{epsilon} by more than a relative \eqn{10^{-6}}, or is
#'   0 or not finite; see \emph{Precision}) and \code{precision_bound}
#'   (numeric, one element per direction: the estimated relative rounding
#'   error of the distance at the point, \code{NA} where there is no point;
#'   a conservative estimate, not a guarantee).
#'
#' @references
#' Wilson KJ et al. (2022). Bayesian sample size determination for diagnostic
#' accuracy studies. \emph{Stat Med} 41:2908-2922. \doi{10.1002/sim.9393}
#' @seealso \code{\link{bam_sensitivity_local}}
#' @examples
#' ct <- bam_prior_contour(c(9.9, 1.1), epsilon = 0.00354, n_phi = 36)
#' range(ct$a)
#' range(ct$b)
#' @export
bam_prior_contour <- function(prior, epsilon, n_phi = 72L, phi = NULL,
                              distance = c("H2", "H")) {
  distance <- match.arg(distance)
  .chk_beta_prior(prior, "prior")
  .bam_check_epsilon(if (missing(epsilon)) NULL else epsilon)
  grid <- .bam_check_grid(n_phi, phi)
  a0 <- as.numeric(prior[1])
  b0 <- as.numeric(prior[2])

  r <- vapply(grid, .bam_contour_radius, numeric(1), a0 = a0, b0 = b0,
              epsilon = epsilon, distance = distance)
  a <- a0 + r * cos(grid)
  b <- b0 + r * sin(grid)
  dd <- rep(NA_real_, length(grid))
  okp <- !is.na(r)
  dd[okp] <- .bam_dist(a[okp], b[okp], a0, b0, distance)
  if (any(!okp)) {
    warning("No contour point found for ", sum(!okp), " of ", length(grid),
            " direction(s) (phi = ",
            paste(format(grid[!okp], digits = 4), collapse = ", "),
            "): the distance did not reach epsilon = ", epsilon,
            " within the admissible radius; epsilon may be too close to 1.",
            call. = FALSE)
  }
  # Precision of the points found: the distance recomputed at the returned
  # point must equal epsilon within a relative tolerance. With very
  # concentrated priors and a very small epsilon the Bhattacharyya
  # coefficient is 1 to within rounding error, the distance is lost to
  # cancellation (it can be exactly 0 or non-finite), and the radius found is
  # not the radius of the contour.
  # Besides that recomputation, which is itself subject to the rounding it
  # tries to detect, the rounding error of the distance is estimated at every
  # point (.bam_contour_precision()) and the point is flagged when that
  # estimate, relative to the distance, exceeds the tolerance. The estimate
  # does not depend on which directions happen to round well, so the flags
  # are deterministic.
  precision_bound <- rep(NA_real_, length(grid))
  precision_bound[okp] <- .bam_contour_precision(a[okp], b[okp], a0, b0,
                                                 dd[okp], distance)
  imprecise <- rep(FALSE, length(grid))
  imprecise[okp] <- !is.finite(dd[okp]) | dd[okp] == 0 |
    abs(dd[okp] - epsilon) > .bam_contour_tol * epsilon |
    is.na(precision_bound[okp]) | precision_bound[okp] > .bam_contour_tol
  if (any(imprecise)) {
    warning("The numerical precision is not sufficient for ", sum(imprecise),
            " of ", length(grid), " direction(s) (phi = ",
            paste(format(grid[imprecise], digits = 4), collapse = ", "),
            "): the estimated rounding error of the distance, relative to ",
            "epsilon = ", format(epsilon), ", exceeds ",
            format(.bam_contour_tol), ", or the distance recomputed at the ",
            "returned point differs from epsilon by more than that relative ",
            "tolerance (or is 0 or not finite), for the prior c(", a0, ", ",
            b0, ") and distance = \"", distance,
            "\". Try a larger epsilon or the other distance. The points are ",
            "returned but flagged in attr(, \"imprecise\").",
            call. = FALSE)
  }
  out <- data.frame(phi = grid, r = r, a = a, b = b, distance = dd)
  attr(out, "prior") <- c(a0, b0)
  attr(out, "epsilon") <- epsilon
  attr(out, "distance_type") <- distance
  attr(out, "imprecise") <- imprecise
  attr(out, "precision_bound") <- precision_bound
  out
}

# ---------------------------------------------------------------------------
# bam_sensitivity_local()
# ---------------------------------------------------------------------------

# One contour point: re-apply the pilot (if any) to the perturbed initial
# prior, then A at n_ref and n*.
.bam_local_point <- function(job) {
  pri <- list(se = job$base_se, sp = job$base_sp, prev = job$base_prev)
  pri[[job$which]] <- job$init
  if (!is.null(job$pilot)) {
    po <- do.call(bam_prior_from_pilot, c(
      job$pilot,
      list(base_prior_se = pri$se, base_prior_sp = pri$sp,
           base_prior_prev = pri$prev, a0 = job$a0, design = job$design,
           use_in = job$use_in)
    ))
    pargs <- list(priors = po)
    used <- switch(job$which, se = po$prior_se, sp = po$prior_sp,
                   prev = po$prior_prev)
  } else {
    pargs <- list(prior_se = pri$se, prior_sp = pri$sp,
                  prior_prev = pri$prev)
    used <- job$init
  }
  sink <- .bam_sens_sink()
  run <- function(extra) {
    .bam_sens_quiet(do.call(bam_sample_size,
                            c(pargs, job$args, extra, list(method = "exact"))),
                    sink)
  }
  r_n <- tryCatch(run(list()), error = function(e) e)
  if (inherits(r_n, "error")) return(list(error = conditionMessage(r_n)))
  r_a <- if (is.null(job$n_ref)) {
    NULL
  } else {
    tryCatch(run(list(N_range = job$n_ref)), error = function(e) e)
  }
  if (inherits(r_a, "error")) return(list(error = conditionMessage(r_a)))
  list(
    a_post = used[1], b_post = used[2],
    n_star = if (isTRUE(r_n$target_reached)) as.integer(r_n$N_total)
             else NA_integer_,
    target_reached = isTRUE(r_n$target_reached),
    A_nref = if (is.null(r_a)) NA_real_ else r_a$max_assurance_evaluated,
    A_at_nstar = if (isTRUE(r_n$target_reached)) r_n$joint_assurance
                 else NA_real_,
    width_type = r_n$width_type, assurance_for = r_n$assurance_for,
    alpha_ci = r_n$alpha_ci, target_assurance = r_n$target_assurance,
    delta_se = r_n$delta_se, delta_sp = r_n$delta_sp,
    prior_only = .bam_sens_sunk(sink)
  )
}

#' Local Sensitivity to the Initial Prior by Hellinger Distance
#'
#' Local prior sensitivity analysis in the style of Wilson et al. (2022,
#' section 4.2), for \code{\link{bam_sample_size}(method = "exact")}. One
#' parameter at a time (Se, Sp or prevalence) the \emph{initial} Beta prior
#' (the one before the pilot data) is moved along the contour of points at
#' distance \code{epsilon} from it (see \code{\link{bam_prior_contour}}); for
#' each point the pilot is applied again with
#' \code{\link{bam_prior_from_pilot}} (if \code{pilot} is given), and the
#' function computes the assurance at a reference sample size \code{n_ref}
#' and the smallest sample size \eqn{n^*} that meets the target. The
#' summary gives the minimum and maximum of both over the contour and the
#' initial priors that produce them.
#'
#' @section Conventions:
#' All conventions of \code{bam_sample_size()} apply and are reported in the
#' result: width criterion (\code{width_type}), objective
#' (\code{assurance_for}, joint or marginal), credible level
#' (\code{alpha_ci}) and target (\code{target_assurance}), passed through
#' \code{...}. The sample size is the \emph{total} \eqn{N}. With a pilot, the
#' perturbed initial prior is updated with the pilot counts (power prior of
#' weight \code{a0}, \code{design}, \code{use_in}; see
#' \code{bam_prior_from_pilot()}), and the updated prior is used as design
#' prior and, with \code{use_in = "both"}, as analysis prior. Without a pilot
#' the perturbed prior is used as given as design prior, and the analysis
#' prior equals the design prior. The distance is the one chosen with
#' \code{distance} (default \code{"H2"}, \eqn{1 - BC}); see
#' \code{bam_prior_contour()}, including the remark on the value 0.00354.
#'
#' @section Differences from the published procedure:
#' The radius is solved exactly for every direction (a contour, not the
#' nearest point of a coarse grid), the assurance and \eqn{n^*} are exact
#' (no optimiser, no Monte Carlo), and the minimum and maximum of \eqn{n^*}
#' are taken over the whole contour. The success rule of
#' \code{bam_sample_size()} (exact enumeration of the Beta-Binomial outcomes)
#' is not the rule of the code of Wilson et al., so numerical results may
#' differ from theirs.
#'
#' @section Discarded directions:
#' A direction is left out of the minima and maxima (and of \code{table}) when
#' no contour point is found (no root) or when the point found does not reach
#' \code{epsilon} with enough numerical precision (see
#' \code{\link{bam_prior_contour}}, \emph{Precision}). Their numbers are
#' reported in \code{summary} (\code{n_dir_no_root}, \code{n_dir_imprecise})
#' and in the printed output. If no direction is left, an error is raised.
#'
#' @section Warning about the analysis prior:
#' The warning of \code{\link{bam_sample_size}} about an analysis prior that
#' already meets its target width with no data is not emitted by each internal
#' call (there are one or two per contour point, and the one at the base
#' prior). A single warning gives the number of contour points affected, and
#' whether the base prior is, identically with \code{parallel = TRUE} and
#' \code{parallel = FALSE}.
#'
#' @param base_prior_se,base_prior_sp,base_prior_prev Initial Beta priors
#'   \code{c(a, b)} (before the pilot) for Se, Sp and prevalence.
#' @param which Parameter whose initial prior is perturbed, one of
#'   \code{"se"} (default), \code{"sp"}, \code{"prev"}. One at a time.
#' @param epsilon Distance of the perturbed priors to the initial prior, in
#'   (0, 1). No default.
#' @param n_phi,phi,distance As in \code{\link{bam_prior_contour}}
#'   (72 equally spaced directions in \eqn{[-\pi, \pi)} by default).
#' @param pilot \code{NULL}, or a list \code{list(tp, fn, tn, fp)} of pilot
#'   counts to apply to each perturbed initial prior.
#' @param a0,design,use_in Passed to \code{\link{bam_prior_from_pilot}}
#'   (ignored when \code{pilot = NULL}).
#' @param n_ref Reference total sample size at which the assurance is
#'   computed. \code{NULL} (default) uses the \eqn{n^*} of the base prior.
#' @param parallel,cores As in \code{\link{bam_sensitivity}}.
#' @param ... Further arguments of \code{\link{bam_sample_size}} (widths
#'   \code{delta_se}, \code{delta_sp}, \code{width_type}, \code{assurance_for},
#'   \code{alpha_ci}, \code{target_assurance}, \code{N_max}, ...). The prior
#'   arguments, \code{N_range} and \code{method} cannot be passed. As in
#'   \code{\link{bam_sensitivity}}, \code{B = 1000} and \code{n_range = 20}
#'   are used for the unused legacy diagnostics unless given.
#'
#' @return An object of class \code{"dta_bam_sens_local"}, a list with
#'   \describe{
#'     \item{table}{One row per contour point: \code{phi}, \code{r},
#'       \code{a_init}, \code{b_init} (perturbed initial prior),
#'       \code{distance}, \code{a_post}, \code{b_post} (the design prior of
#'       the perturbed parameter after the pilot), \code{A_nref} (assurance
#'       at \code{n_ref}), \code{n_star} (\code{NA} if the target is not
#'       reached up to \code{N_max}), \code{target_reached}.}
#'     \item{summary}{One row: \code{A_min}, \code{A_max} with the initial
#'       priors that give them (\code{A_min_a_init}, \code{A_min_b_init},
#'       \code{A_max_a_init}, \code{A_max_b_init}), and \code{n_star_min},
#'       \code{n_star_max} with their initial priors, plus \code{n_star_NA}
#'       (points where the target was not reached), \code{n_dir_no_root}
#'       (directions without a contour point) and \code{n_dir_imprecise}
#'       (directions discarded for insufficient numerical precision). Ties go
#'       to the first point of the grid.}
#'     \item{base}{Results at the unperturbed initial prior: \code{n_star},
#'       \code{A_nref}, \code{A_at_nstar} and the \code{priors} used.}
#'     \item{n_ref, which, epsilon, distance, conventions}{Inputs and the
#'       conventions actually used.}
#'   }
#'
#' @references
#' Wilson KJ et al. (2022). Bayesian sample size determination for diagnostic
#' accuracy studies. \emph{Stat Med} 41:2908-2922. \doi{10.1002/sim.9393}
#' @seealso \code{\link{bam_prior_contour}}, \code{\link{bam_prior_from_pilot}},
#'   \code{\link{bam_sensitivity}}
#' @examples
#' # VAP-style example (Wilson et al. 2022), only 6 directions to be fast;
#' # use the default 72 directions in a real analysis.
#' res <- suppressWarnings(bam_sensitivity_local(
#'   base_prior_se = c(9.9, 1.1), base_prior_sp = c(2, 2),
#'   base_prior_prev = c(12, 43), which = "se", epsilon = 0.00354,
#'   n_phi = 6, pilot = list(tp = 16, fn = 1, tn = 20, fp = 35),
#'   delta_se = 0.16, assurance_for = "se", n_ref = 106
#' ))
#' res$summary
#' res
#' @export
bam_sensitivity_local <- function(base_prior_se, base_prior_sp,
                                  base_prior_prev,
                                  which = c("se", "sp", "prev"),
                                  epsilon, n_phi = 72L, phi = NULL,
                                  distance = c("H2", "H"), pilot = NULL,
                                  a0 = 1,
                                  design = c("cohort", "case_control"),
                                  use_in = c("both", "design"),
                                  n_ref = NULL, parallel = FALSE,
                                  cores = 2L, ...) {
  if (missing(base_prior_se) || missing(base_prior_sp) ||
      missing(base_prior_prev)) {
    stop("`base_prior_se`, `base_prior_sp` and `base_prior_prev` are ",
         "required (initial Beta priors c(a, b), before the pilot).",
         call. = FALSE)
  }
  which <- match.arg(which)
  distance <- match.arg(distance)
  design <- match.arg(design)
  use_in <- match.arg(use_in)
  .chk_beta_prior(base_prior_se, "base_prior_se")
  .chk_beta_prior(base_prior_sp, "base_prior_sp")
  .chk_beta_prior(base_prior_prev, "base_prior_prev")
  .bam_check_epsilon(if (missing(epsilon)) NULL else epsilon)
  .bam_sens_check_parallel(parallel, cores)
  grid <- .bam_check_grid(n_phi, phi)

  if (!is.null(pilot)) {
    if (!is.list(pilot) ||
        !setequal(names(pilot), c("tp", "fn", "tn", "fp"))) {
      stop("`pilot` must be NULL or a list with elements tp, fn, tn, fp.",
           call. = FALSE)
    }
    pilot <- pilot[c("tp", "fn", "tn", "fp")]
  }
  if (!is.null(n_ref)) {
    if (!is.numeric(n_ref) || length(n_ref) != 1L || !is.finite(n_ref) ||
        n_ref != round(n_ref) || n_ref < 2) {
      stop("`n_ref` must be NULL or a single whole number >= 2.",
           call. = FALSE)
    }
    n_ref <- as.integer(n_ref)
  }

  dots <- list(...)
  .bam_sens_check_dots(
    dots, reserved = c(.bam_sens_prior_args, "N_range"),
    where = "bam_sensitivity_local()"
  )
  args <- c(.bam_sens_cheap_defaults[setdiff(names(.bam_sens_cheap_defaults),
                                             names(dots))], dots)

  base_priors <- list(se = as.numeric(base_prior_se),
                      sp = as.numeric(base_prior_sp),
                      prev = as.numeric(base_prior_prev))
  mk_job <- function(init, n_ref_use) {
    list(init = init, which = which, base_se = base_priors$se,
         base_sp = base_priors$sp, base_prev = base_priors$prev,
         pilot = pilot, a0 = a0, design = design, use_in = use_in,
         n_ref = n_ref_use, args = args)
  }

  # Unperturbed base prior
  base_res <- .bam_local_point(mk_job(base_priors[[which]], n_ref))
  if (!is.null(base_res$error)) {
    stop("Base prior: ", base_res$error, call. = FALSE)
  }
  base_prior_only <- length(base_res$prior_only) > 0L
  if (is.null(n_ref)) {
    if (!isTRUE(base_res$target_reached)) {
      stop("The base prior does not reach the target up to N_max, so there ",
           "is no default `n_ref`; give `n_ref` or increase `N_max`.",
           call. = FALSE)
    }
    n_ref <- base_res$n_star
    base_res <- .bam_local_point(mk_job(base_priors[[which]], n_ref))
    if (!is.null(base_res$error)) {
      stop("Base prior: ", base_res$error, call. = FALSE)
    }
    base_prior_only <- base_prior_only || length(base_res$prior_only) > 0L
  }

  ct <- bam_prior_contour(base_priors[[which]], epsilon = epsilon, phi = grid,
                          distance = distance)
  # Directions discarded from the minimum and maximum: no root found
  # (NA point) or insufficient numerical precision (see bam_prior_contour()).
  imprecise <- attr(ct, "imprecise")
  n_dir_no_root <- sum(is.na(ct$a))
  n_dir_imprecise <- sum(!is.na(ct$a) & imprecise)
  ct <- ct[!is.na(ct$a) & !imprecise, , drop = FALSE]
  if (nrow(ct) == 0L) {
    stop("The contour has no valid point (", n_dir_no_root,
         " direction(s) without a root, ", n_dir_imprecise,
         " with insufficient numerical precision).", call. = FALSE)
  }

  jobs <- lapply(seq_len(nrow(ct)), function(i) {
    mk_job(c(ct$a[i], ct$b[i]), n_ref)
  })
  res <- .bam_sens_lapply(jobs, .bam_local_point, parallel = parallel,
                          cores = cores)
  err <- vapply(res, function(r) !is.null(r$error), logical(1))
  # One summary warning (the internal calls do not emit the prior-only
  # warning of bam_sample_size()): how many contour points, and whether the
  # base prior, have an analysis prior that already meets the width with no
  # data. It is emitted before the error check so that it is not lost when a
  # contour point fails.
  n_prior_only <- sum(vapply(res, function(r) length(r$prior_only) > 0L,
                             logical(1)))
  if (n_prior_only > 0L || base_prior_only) {
    warning(
      "The analysis prior of a used arm already has a prior-only interval ",
      "width within its target delta (the width criterion is met before any ",
      "data are observed, so success is dominated by the analysis prior) in ",
      n_prior_only, " of ", length(res), " contour point(s)",
      if (base_prior_only) " and in the base prior" else "",
      ". Use a weaker or flat analysis prior if this is not intended. See ",
      "bam_sample_size() for the detail of each call.", call. = FALSE
    )
  }
  if (any(err)) {
    stop("Contour point ", which(err)[1L], ": ", res[[which(err)[1L]]]$error,
         call. = FALSE)
  }

  tb <- data.frame(
    phi = ct$phi, r = ct$r, a_init = ct$a, b_init = ct$b,
    distance = ct$distance,
    a_post = vapply(res, function(r) r$a_post, numeric(1)),
    b_post = vapply(res, function(r) r$b_post, numeric(1)),
    A_nref = vapply(res, function(r) r$A_nref, numeric(1)),
    n_star = vapply(res, function(r) r$n_star, integer(1)),
    target_reached = vapply(res, function(r) r$target_reached, logical(1)),
    stringsAsFactors = FALSE
  )

  i_amin <- which.min(tb$A_nref)
  i_amax <- which.max(tb$A_nref)
  has_n <- any(!is.na(tb$n_star))
  i_nmin <- if (has_n) which.min(tb$n_star) else NA_integer_
  i_nmax <- if (has_n) which.max(tb$n_star) else NA_integer_
  pick <- function(i, col) if (is.na(i)) NA_real_ else tb[[col]][i]
  sm <- data.frame(
    A_min = tb$A_nref[i_amin],
    A_min_a_init = tb$a_init[i_amin], A_min_b_init = tb$b_init[i_amin],
    A_max = tb$A_nref[i_amax],
    A_max_a_init = tb$a_init[i_amax], A_max_b_init = tb$b_init[i_amax],
    n_star_min = if (has_n) tb$n_star[i_nmin] else NA_integer_,
    n_star_min_a_init = pick(i_nmin, "a_init"),
    n_star_min_b_init = pick(i_nmin, "b_init"),
    n_star_max = if (has_n) tb$n_star[i_nmax] else NA_integer_,
    n_star_max_a_init = pick(i_nmax, "a_init"),
    n_star_max_b_init = pick(i_nmax, "b_init"),
    n_star_NA = sum(is.na(tb$n_star)),
    n_dir_no_root = as.integer(n_dir_no_root),
    n_dir_imprecise = as.integer(n_dir_imprecise)
  )

  base_priors_used <- if (!is.null(pilot)) {
    do.call(bam_prior_from_pilot, c(pilot, list(
      base_prior_se = base_priors$se, base_prior_sp = base_priors$sp,
      base_prior_prev = base_priors$prev, a0 = a0, design = design,
      use_in = use_in)))
  } else {
    NULL
  }

  structure(
    list(
      table = tb, summary = sm,
      base = list(n_star = base_res$n_star, A_nref = base_res$A_nref,
                  A_at_nstar = base_res$A_at_nstar,
                  initial_priors = base_priors, priors = base_priors_used),
      n_ref = n_ref, which = which, epsilon = epsilon, distance = distance,
      grid = list(n_phi = if (is.null(phi)) length(grid) else NULL,
                  phi = grid),
      pilot = pilot, a0 = if (is.null(pilot)) NULL else a0,
      design = if (is.null(pilot)) NULL else design,
      use_in = if (is.null(pilot)) NULL else use_in,
      conventions = list(
        width_type = base_res$width_type,
        assurance_for = base_res$assurance_for,
        alpha_ci = base_res$alpha_ci,
        target_assurance = base_res$target_assurance,
        delta_se = base_res$delta_se, delta_sp = base_res$delta_sp,
        sample_size = "total N"
      )
    ),
    class = "dta_bam_sens_local"
  )
}

#' @export
print.dta_bam_sens_local <- function(x, ...) {
  cv <- x$conventions
  nm <- c(se = "Se", sp = "Sp", prev = "prevalence")[[x$which]]
  dn <- if (x$distance == "H2") "H2 = 1 - Bhattacharyya coefficient"
        else "H = sqrt(1 - Bhattacharyya coefficient)"
  cat("\nLocal prior sensitivity (Hellinger contour), exact BAM\n\n")
  cat("  Perturbed initial prior:", nm, " (one parameter at a time)\n")
  cat("  Distance:", dn, "| epsilon =", format(x$epsilon), "|",
      nrow(x$table), "contour points\n")
  if (is.null(x$pilot)) {
    cat("  Pilot: none (perturbed prior used as is; analysis prior = design",
        "prior)\n")
  } else {
    p <- x$pilot
    cat("  Pilot re-applied: tp =", p$tp, " fn =", p$fn, " tn =", p$tn,
        " fp =", p$fp, "| a0 =", format(x$a0), "|", x$design, "| use_in =",
        x$use_in, "\n")
  }
  cat("  Conventions: ", cv$width_type, " (delta_se = ", format(cv$delta_se),
      ", delta_sp = ", format(cv$delta_sp), "), assurance_for = ",
      cv$assurance_for, ", credible level = ", format(cv$alpha_ci),
      ", target assurance = ", format(cv$target_assurance), ", ",
      cv$sample_size, "\n", sep = "")
  cat("  Reference sample size n_ref =", x$n_ref, "\n\n")
  s <- x$summary
  cat(sprintf("  Base prior:  A(n_ref) = %.4f   n* = %s\n",
              x$base$A_nref, format(x$base$n_star)))
  cat(sprintf("  A(n_ref):    min = %.4f  at (a, b) = (%.4f, %.4f)\n",
              s$A_min, s$A_min_a_init, s$A_min_b_init))
  cat(sprintf("               max = %.4f  at (a, b) = (%.4f, %.4f)\n",
              s$A_max, s$A_max_a_init, s$A_max_b_init))
  cat(sprintf("  n*:          min = %s  at (a, b) = (%.4f, %.4f)\n",
              format(s$n_star_min), s$n_star_min_a_init, s$n_star_min_b_init))
  cat(sprintf("               max = %s  at (a, b) = (%.4f, %.4f)\n",
              format(s$n_star_max), s$n_star_max_a_init, s$n_star_max_b_init))
  if (s$n_star_NA > 0) {
    cat("  n* not reached up to N_max at", s$n_star_NA, "point(s)\n")
  }
  n_disc <- s$n_dir_no_root + s$n_dir_imprecise
  if (n_disc > 0) {
    cat("  Directions discarded from the minima and maxima: ", n_disc,
        " (", s$n_dir_no_root, " without a root, ", s$n_dir_imprecise,
        " with insufficient numerical precision)\n", sep = "")
  }
  cat("\n")
  invisible(x)
}
