#' Install necessary packages
#' 

library(devtools)

Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)

install_github(repo = "DTUAqua/DATRAS/DATRAS")
# remotes::install_github("DTUAqua/DATRAS/DATRAS") # for MRAN


install_github(repo = "casperwberg/surveyIndex/surveyIndex")
# remotes::install_github("casperwberg/surveyIndex/surveyIndex") # for MRAN