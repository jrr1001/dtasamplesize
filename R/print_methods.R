#' Print Method for dtasamplesize Objects
#'
#' @param x An object of class \code{"dtasamplesize"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.dtasamplesize <- function(x, ...) {
  cat("\n", x$method, "\n\n")
  if (!is.null(x$status) && x$status != "converged") {
    cat("  status:", x$status,
        "-- NO N SATISFIES target_assurance; N_total/n_diseased are NA.\n")
    cat("  (see the warning issued by this call for diagnostic detail)\n")
  }
  cat("  n_diseased:", x$n_diseased, "\n")
  if (!is.null(x$n_non_diseased))
    cat("  n_non_diseased:", x$n_non_diseased, "\n")
  cat("  N_total:", x$n_total, "\n")
  if (!is.null(x$joint_assurance))
    cat("  joint_assurance:", x$joint_assurance, "\n")
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
