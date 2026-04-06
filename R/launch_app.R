#' Launch the dtasamplesize Shiny App
#'
#' Opens an interactive Shiny application for exploring sample size
#' calculations for diagnostic test accuracy studies.
#'
#' @return Runs the Shiny app (does not return a value).
#' @export
launch_app <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with: ",
         "install.packages('shiny')")
  }
  appDir <- system.file("shiny", package = "dtasamplesize")
  if (appDir == "") {
    stop("Could not find Shiny app directory. ",
         "Try re-installing 'dtasamplesize'.")
  }
  shiny::runApp(appDir, display.mode = "normal")
}
