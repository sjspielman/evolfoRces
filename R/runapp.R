#' Run the evolfoRces Shiny application locally
#'
#' @export
run_evolfoRces <- function() {
  appDir <- system.file("evolfoRces", package = "evolfoRces")
  shiny::runApp(appDir, display.mode = "normal")
}


#' Run the evolfoRces Shiny application locally, allowing for people to forget about cap "R"
#'
#' @export
run_evolforces <- function() {
  appDir <- system.file("evolfoRces", package = "evolfoRces")
  shiny::runApp(appDir, display.mode = "normal")
}
