#' Plot the Assurance Curve from a Unified Sample Size Search
#'
#' Plots achieved joint assurance against total sample size N, using the
#' \code{grid_results} element of an object returned by
#' \code{\link{ss_unified}}. A dashed horizontal line marks the target
#' assurance, and a dotted vertical line marks the optimal N, labeled with
#' its value.
#'
#' @param x An object of class \code{"dtasamplesize"} returned by
#'   \code{\link{ss_unified}}, called with \code{full_grid = TRUE} so that
#'   \code{grid_results} covers the whole \code{N_range} rather than
#'   stopping at the first N that reaches the target.
#' @param title Optional plot title. Default \code{NULL} (no title): the
#'   package leaves the caption to the caller, as is customary for
#'   journal figures.
#' @param ... Currently unused; present for method consistency.
#' @return A \code{ggplot} object, which can be further customized with
#'   additional \pkg{ggplot2} layers.
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   result <- suppressWarnings(ss_unified(
#'     delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
#'     N_range = seq(300, 900, by = 100), B = 200, full_grid = TRUE
#'   ))
#'   plot_assurance_curve(result)
#' }
#' }
#' @seealso \code{\link{ss_unified}}
#' @export
plot_assurance_curve <- function(x, title = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is needed for this function. Install it with install.packages('ggplot2').", call. = FALSE)

  if (!inherits(x, "dtasamplesize")) {
    stop(
      "plot_assurance_curve(): 'x' must be an object of class ",
      "'dtasamplesize', as returned by ss_unified().",
      call. = FALSE
    )
  }

  grid <- x$grid_results
  if (is.null(grid) || nrow(grid) < 2) {
    stop(
      "plot_assurance_curve(): 'x$grid_results' is missing or has fewer ",
      "than 2 rows. Re-run ss_unified() with 'full_grid = TRUE' to obtain ",
      "the assurance curve over the full N_range.",
      call. = FALSE
    )
  }
  if (!all(c("N", "assurance") %in% names(grid))) {
    stop(
      "plot_assurance_curve(): 'x$grid_results' must have columns 'N' and ",
      "'assurance'.",
      call. = FALSE
    )
  }

  target <- unified_target_assurance(x)
  optimal_N <- unified_optimal_N(x, target)

  p <- ggplot2::ggplot(
    grid, ggplot2::aes(x = .data[["N"]], y = .data[["assurance"]])
  ) +
    ggplot2::geom_line(color = "grey30", linewidth = 0.7) +
    ggplot2::geom_point(color = "grey30", size = 1.6) +
    ggplot2::labs(
      x = "Total sample size (N)",
      y = "Joint assurance",
      title = title
    ) +
    ggplot2::theme_minimal()

  if (is.finite(target)) {
    p <- p + ggplot2::geom_hline(
      yintercept = target, linetype = "dashed", color = "grey50"
    )
  }

  if (is.finite(optimal_N)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = optimal_N, linetype = "dotted", color = "black"
      ) +
      ggplot2::annotate(
        "text",
        x = optimal_N,
        y = min(grid$assurance, na.rm = TRUE),
        label = paste0("N = ", optimal_N),
        hjust = -0.1, vjust = -0.4, color = "black", size = 3.2
      )
  }

  p
}

#' Plot the Variance Inflation Heatmap for an Imperfect Reference Standard
#'
#' Plots the \code{sensitivity_table} element of an object returned by
#' \code{\link{ss_imperfect_ref}} as a heatmap of the variance inflation
#' factor over reference-standard sensitivity and specificity.
#'
#' @param x An object of class \code{"dtasamplesize"} returned by
#'   \code{\link{ss_imperfect_ref}}, called with
#'   \code{sensitivity_table = TRUE}.
#' @param title Optional plot title. Default \code{NULL} (no title).
#' @param midpoint Midpoint of the diverging fill scale. Default
#'   \code{NULL}, which uses the median inflation factor in the table.
#' @param ... Currently unused; present for method consistency.
#' @return A \code{ggplot} object, which can be further customized with
#'   additional \pkg{ggplot2} layers.
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   result <- ss_imperfect_ref(B = 0, sensitivity_table = TRUE)
#'   plot_inflation_heatmap(result)
#' }
#' }
#' @seealso \code{\link{ss_imperfect_ref}}
#' @export
plot_inflation_heatmap <- function(x, title = NULL, midpoint = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is needed for this function. Install it with install.packages('ggplot2').", call. = FALSE)

  if (!inherits(x, "dtasamplesize")) {
    stop(
      "plot_inflation_heatmap(): 'x' must be an object of class ",
      "'dtasamplesize', as returned by ss_imperfect_ref().",
      call. = FALSE
    )
  }

  tbl <- x$sensitivity_table
  if (is.null(tbl)) {
    stop(
      "plot_inflation_heatmap(): 'x$sensitivity_table' is missing. ",
      "Re-run ss_imperfect_ref() with 'sensitivity_table = TRUE'.",
      call. = FALSE
    )
  }
  if (!all(c("Se_ref", "Sp_ref") %in% names(tbl))) {
    stop(
      "plot_inflation_heatmap(): 'x$sensitivity_table' must have columns ",
      "'Se_ref' and 'Sp_ref'.",
      call. = FALSE
    )
  }
  fill_col <- if ("inflation_factor" %in% names(tbl)) "inflation_factor" else "VIF"
  if (!(fill_col %in% names(tbl))) {
    stop(
      "plot_inflation_heatmap(): 'x$sensitivity_table' must have a 'VIF' ",
      "or 'inflation_factor' column.",
      call. = FALSE
    )
  }

  if (is.null(midpoint)) midpoint <- stats::median(tbl[[fill_col]], na.rm = TRUE)

  ggplot2::ggplot(
    tbl,
    ggplot2::aes(
      x = .data[["Se_ref"]], y = .data[["Sp_ref"]], fill = .data[[fill_col]]
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = "#f0f0f0", mid = "#969696", high = "#252525",
      midpoint = midpoint, name = "Inflation\nfactor"
    ) +
    ggplot2::labs(
      x = "Reference standard sensitivity (Se_ref)",
      y = "Reference standard specificity (Sp_ref)",
      title = title
    ) +
    ggplot2::theme_minimal()
}

#' Plot a Method Comparison from a Unified Sample Size Search
#'
#' Plots the \code{comparison} element of an object returned by
#' \code{\link{ss_unified}} as a bar chart of required total sample size N
#' by method, with the N value labeled above each bar.
#'
#' @param x An object of class \code{"dtasamplesize"} returned by
#'   \code{\link{ss_unified}}.
#' @param title Optional plot title. Default \code{NULL} (no title).
#' @param ... Currently unused; present for method consistency.
#' @return A \code{ggplot} object, which can be further customized with
#'   additional \pkg{ggplot2} layers.
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   result <- suppressWarnings(ss_unified(
#'     delta_se = 0.07, delta_sp = 0.05, delta_auc = 0,
#'     N_range = seq(300, 900, by = 100), B = 200
#'   ))
#'   plot_method_comparison(result)
#' }
#' }
#' @seealso \code{\link{ss_unified}}
#' @export
plot_method_comparison <- function(x, title = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is needed for this function. Install it with install.packages('ggplot2').", call. = FALSE)

  if (!inherits(x, "dtasamplesize")) {
    stop(
      "plot_method_comparison(): 'x' must be an object of class ",
      "'dtasamplesize', as returned by ss_unified().",
      call. = FALSE
    )
  }

  comp <- x$comparison
  if (is.null(comp)) {
    stop(
      "plot_method_comparison(): 'x$comparison' is missing. This element ",
      "is produced by ss_unified().",
      call. = FALSE
    )
  }
  if (!all(c("method", "N") %in% names(comp))) {
    stop(
      "plot_method_comparison(): 'x$comparison' must have columns 'method' ",
      "and 'N'.",
      call. = FALSE
    )
  }

  ggplot2::ggplot(
    comp, ggplot2::aes(x = .data[["method"]], y = .data[["N"]])
  ) +
    ggplot2::geom_col(fill = "grey40", width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data[["N"]]),
      vjust = -0.4, color = "black", size = 3.5
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = "Method", y = "Required total sample size (N)", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
}

# --- internal helpers ------------------------------------------------------

# Recovers the target assurance used by an ss_unified() call: the
# x$target_assurance field (present since the field was added to
# ss_unified()'s return value), falling back to the formal default of
# ss_unified() itself for older result objects that predate the field.
unified_target_assurance <- function(x) {
  if (!is.null(x$target_assurance)) {
    val <- suppressWarnings(as.numeric(x$target_assurance))
    if (!is.na(val)) return(val)
  }

  default <- formals(ss_unified)[["target_assurance"]]
  if (!is.null(default)) {
    val <- suppressWarnings(as.numeric(default))
    if (!is.na(val)) return(val)
  }

  NA_real_
}

# Best-effort recovery of the optimal N: prefers the N_effective field
# already produced by ss_unified(), falling back to the first grid_results
# row whose assurance meets the target.
unified_optimal_N <- function(x, target) {
  if (!is.null(x$N_effective)) {
    val <- suppressWarnings(as.numeric(x$N_effective))
    if (!is.na(val)) return(val)
  }

  grid <- x$grid_results
  if (!is.null(grid) && is.finite(target)) {
    idx <- which(grid$assurance >= target)
    if (length(idx)) return(grid$N[idx[1]])
  }

  NA_real_
}
