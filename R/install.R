#' Install necessary packages
#' 

library(devtools)

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)

install_github(repo = "DTUAqua/DATRAS/DATRAS")
install_github(repo = "casperwberg/surveyIndex/surveyIndex")
