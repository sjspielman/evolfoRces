#' Run the evolfoRces Shiny application locally
#'
#' @export
run <- function() {
  appDir <- system.file("evolfoRces", package = "evolfoRces")
  shiny::runApp(appDir, display.mode = "normal")
}
