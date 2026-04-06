#' Print Method for dtasamplesize Objects
#'
#' @param x An object of class \code{"dtasamplesize"}.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}.
#' @export
print.dtasamplesize <- function(x, ...) {
  cat("\n", x$method, "\n\n")
  cat("  n_diseased:", x$n_diseased, "\n")
  if (!is.null(x$n_non_diseased))
    cat("  n_non_diseased:", x$n_non_diseased, "\n")
  cat("  N_total:", x$n_total, "\n")
  if (!is.null(x$results)) {
    cat("\n")
    print(x$results)
  }
  invisible(x)
}
