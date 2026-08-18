#' Sensitivity Analysis Across Prior Scenarios
#'
#' @title Sample Size Sensitivity Analysis Over Prior Scenarios
#' @description Runs \code{\link{ss_unified}} once for each of a set of
#'   prior scenarios (e.g. optimistic/pessimistic assumptions about Se, Sp,
#'   and prevalence) and collects the resulting sample sizes into a single
#'   data frame. Replaces the hand-written \code{lapply()} loop over
#'   scenarios that this kind of comparison otherwise requires.
#'
#' @param scenarios Either (a) a \strong{named list} of scenarios, where
#'   each element is itself a list with components \code{prior_se},
#'   \code{prior_sp}, and \code{prior_prev} -- each a length-2 numeric
#'   vector \code{c(shape1, shape2)} giving a Beta prior, passed straight to
#'   \code{\link{ss_unified}} -- or (b) a \strong{data frame} with one row
#'   per scenario and the columns \code{scenario}, \code{prior_se_shape1},
#'   \code{prior_se_shape2}, \code{prior_sp_shape1}, \code{prior_sp_shape2},
#'   \code{prior_prev_shape1}, and \code{prior_prev_shape2}.
#' @param ... Additional arguments passed unchanged to
#'   \code{\link{ss_unified}} for every scenario (e.g. \code{delta_se},
#'   \code{delta_sp}, \code{delta_auc}, \code{target_assurance},
#'   \code{N_range}, \code{B}, \code{seed}). Must not include
#'   \code{prior_se}, \code{prior_sp}, or \code{prior_prev}: those come from
#'   \code{scenarios}.
#' @param verbose Logical. If \code{TRUE}, reports progress as each
#'   scenario is run. Default \code{FALSE}.
#' @return A data frame with one row per scenario, in the order supplied,
#'   and columns:
#'   \describe{
#'     \item{scenario}{Scenario name (list names, or the \code{scenario}
#'       column when \code{scenarios} is a data frame).}
#'     \item{N_total}{Total enrolled sample size from \code{ss_unified()}
#'       for that scenario (loss-adjusted).}
#'     \item{joint_assurance}{Achieved joint assurance at that sample
#'       size.}
#'     \item{N_buderer}{Classical Buderer total N for the same scenario,
#'       for comparison. \code{NA} if the installed \code{ss_unified()}
#'       does not report it.}
#'     \item{n_diseased}{Required number of diseased subjects.}
#'   }
#' @examples
#' \donttest{
#' scenarios <- list(
#'   optimistic = list(prior_se = c(18, 2), prior_sp = c(19, 1),
#'                      prior_prev = c(5, 15)),
#'   pessimistic = list(prior_se = c(14, 6), prior_sp = c(16, 4),
#'                       prior_prev = c(5, 15))
#' )
#' suppressWarnings(sensitivity_analysis(
#'   scenarios,
#'   delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
#'   N_range = seq(300, 900, by = 100), B = 200
#' ))
#' }
#' @seealso \code{\link{ss_unified}}
#' @export
sensitivity_analysis <- function(scenarios, ..., verbose = FALSE) {
  scenario_list <- as_scenario_list(scenarios)
  scenario_names <- names(scenario_list)
  n_scenarios <- length(scenario_list)

  extra_args <- list(...)
  conflicting <- intersect(names(extra_args), c("prior_se", "prior_sp", "prior_prev"))
  if (length(conflicting)) {
    stop(
      "sensitivity_analysis(): ", paste(sQuote(conflicting), collapse = ", "),
      " must be supplied through 'scenarios', not through '...'.",
      call. = FALSE
    )
  }

  N_total <- numeric(n_scenarios)
  joint_assurance <- numeric(n_scenarios)
  N_buderer <- numeric(n_scenarios)
  n_diseased <- numeric(n_scenarios)

  for (i in seq_len(n_scenarios)) {
    sc <- scenario_list[[i]]
    nm <- scenario_names[i]

    if (isTRUE(verbose)) {
      message(sprintf(
        "sensitivity_analysis(): running scenario '%s' (%d/%d)...",
        nm, i, n_scenarios
      ))
    }

    result <- tryCatch(
      do.call(ss_unified, c(
        list(prior_se = sc$prior_se, prior_sp = sc$prior_sp,
             prior_prev = sc$prior_prev),
        extra_args
      )),
      error = function(e) {
        stop(
          sprintf(
            paste0(
              "sensitivity_analysis(): scenario '%s' (position %d of %d) ",
              "failed in ss_unified(): %s"
            ),
            nm, i, n_scenarios, conditionMessage(e)
          ),
          call. = FALSE
        )
      }
    )

    N_total[i] <- result$n_total
    joint_assurance[i] <- result$joint_assurance
    N_buderer[i] <- scenario_field_or_na(result, "N_buderer")
    n_diseased[i] <- result$n_diseased
  }

  data.frame(
    scenario = scenario_names,
    N_total = N_total,
    joint_assurance = joint_assurance,
    N_buderer = N_buderer,
    n_diseased = n_diseased,
    stringsAsFactors = FALSE
  )
}

# --- internal helpers ------------------------------------------------------

# Extracts a named element from a dtasamplesize result, returning NA_real_
# (rather than erroring) when the element is absent, so sensitivity_analysis()
# keeps working against older ss_unified() versions that do not report it.
scenario_field_or_na <- function(x, name) {
  val <- x[[name]]
  if (is.null(val)) NA_real_ else as.numeric(val)
}

# Normalizes the 'scenarios' argument (named list or data frame) into a
# named list of validated scenario lists, each with prior_se/prior_sp/
# prior_prev.
as_scenario_list <- function(scenarios) {
  if (missing(scenarios) || is.null(scenarios)) {
    stop("sensitivity_analysis(): 'scenarios' must be supplied.", call. = FALSE)
  }

  if (is.data.frame(scenarios)) {
    return(scenario_df_to_list(scenarios))
  }

  if (!is.list(scenarios)) {
    stop(
      "sensitivity_analysis(): 'scenarios' must be a named list of prior ",
      "scenarios or a data frame; got an object of class '",
      paste(class(scenarios), collapse = "/"), "'.",
      call. = FALSE
    )
  }

  n <- length(scenarios)
  if (n < 1) {
    stop(
      "sensitivity_analysis(): 'scenarios' must have at least one element.",
      call. = FALSE
    )
  }

  nms <- names(scenarios)
  if (is.null(nms) || any(!nzchar(nms))) {
    bad <- if (is.null(nms)) seq_len(n) else which(!nzchar(nms))
    stop(
      "sensitivity_analysis(): every element of 'scenarios' must be named; ",
      "unnamed element(s) at position(s): ", paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(nms)) {
    dup <- unique(nms[duplicated(nms)])
    stop(
      "sensitivity_analysis(): scenario names must be unique; duplicated ",
      "name(s): ", paste(dup, collapse = ", "), ".",
      call. = FALSE
    )
  }

  for (i in seq_len(n)) {
    scenarios[[i]] <- validate_scenario(scenarios[[i]], nms[i])
  }

  scenarios
}

# Converts the tabular scenario layout (one row per scenario, shape1/shape2
# columns) into the same named-list-of-lists shape used internally.
scenario_df_to_list <- function(df) {
  required_cols <- c("scenario", "prior_se_shape1", "prior_se_shape2",
                      "prior_sp_shape1", "prior_sp_shape2",
                      "prior_prev_shape1", "prior_prev_shape2")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols)) {
    stop(
      "sensitivity_analysis(): 'scenarios' data frame is missing required ",
      "column(s): ", paste(missing_cols, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (nrow(df) < 1) {
    stop(
      "sensitivity_analysis(): 'scenarios' data frame has no rows.",
      call. = FALSE
    )
  }

  nms <- as.character(df$scenario)
  if (any(!nzchar(nms))) {
    stop(
      "sensitivity_analysis(): 'scenarios' data frame has empty 'scenario' ",
      "name(s) in row(s): ", paste(which(!nzchar(nms)), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(nms)) {
    dup <- unique(nms[duplicated(nms)])
    stop(
      "sensitivity_analysis(): 'scenarios' data frame has duplicated ",
      "'scenario' name(s): ", paste(dup, collapse = ", "), ".",
      call. = FALSE
    )
  }

  out <- vector("list", nrow(df))
  names(out) <- nms
  for (i in seq_len(nrow(df))) {
    sc <- list(
      prior_se = c(df$prior_se_shape1[i], df$prior_se_shape2[i]),
      prior_sp = c(df$prior_sp_shape1[i], df$prior_sp_shape2[i]),
      prior_prev = c(df$prior_prev_shape1[i], df$prior_prev_shape2[i])
    )
    out[[i]] <- validate_scenario(sc, nms[i])
  }
  out
}

# Validates a single scenario's prior_se/prior_sp/prior_prev elements,
# naming both the offending scenario and the offending element in the error.
validate_scenario <- function(sc, nm) {
  if (!is.list(sc)) {
    stop(
      sprintf(
        paste0(
          "sensitivity_analysis(): scenario '%s' must be a list with ",
          "elements 'prior_se', 'prior_sp', and 'prior_prev'."
        ),
        nm
      ),
      call. = FALSE
    )
  }

  required <- c("prior_se", "prior_sp", "prior_prev")
  missing_el <- setdiff(required, names(sc))
  if (length(missing_el)) {
    stop(
      sprintf(
        "sensitivity_analysis(): scenario '%s' is missing required element(s): %s.",
        nm, paste(missing_el, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (el in required) {
    v <- sc[[el]]
    if (!is.numeric(v) || length(v) != 2 || anyNA(v) || any(!is.finite(v)) ||
          any(v <= 0)) {
      stop(
        sprintf(
          paste0(
            "sensitivity_analysis(): scenario '%s' element '%s' must be a ",
            "length-2 positive numeric vector c(shape1, shape2); got %s."
          ),
          nm, el, paste(deparse(v), collapse = " ")
        ),
        call. = FALSE
      )
    }
  }

  sc[required]
}
