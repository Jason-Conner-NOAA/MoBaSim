#' Make IDWs for MoBaSim Simple OM
#' Author: Jason Conner (jason.conner@noaa.gov)
#' Date: 11/5/2020
#' 
#' Description:
#' Takes Survey data and a spatial grid (sf) as inputs and outputs simple IDW.


make_IDW <- function(data, var_name, year, input_grid, idp = 6, nmax = 8) {
  ## Should probably add tests for CRS, extent, var_name exists
  sp_yr <- dplyr::filter(data, YEAR==year)

  
  sp_idw <- gstat::idw(formula = eval(parse(text=paste0("sp_yr$",var_name))) ~ 1, 
                       locations = sp_yr, 
                       idp = idp, 
                       newdata = input_grid) %>%
    dplyr::mutate(YEAR = year) %>%
    cbind(input_grid)
  
  
  return(sp_idw)
  
}
