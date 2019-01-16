packages <- c("tidyverse", "shinythemes", "viridis", "cowplot", "shiny", "colourpicker")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}



library(shiny)
runApp(".")