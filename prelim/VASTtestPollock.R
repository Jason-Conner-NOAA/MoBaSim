
# Load packages
library(sumfish)
library(TMB)               
library(VAST)
setwd(file.path(getwd(),"PollockSim"))

# load data set
# see `?load_example` for list of stocks with example data 
# that are installed automatically with `FishStatsUtils`. 
species <- 21740
racebase <- sumfish::getRacebase(2000:2019,'EBS_SHELF')

Data <- sumHaul(racebase) %>%# readRDS("C:/R/VAST2019/haulSum82_19.RDS") %>% ## Hardcoded from saved haulSum
   dplyr::filter(SPECIES_CODE==species)

Data_Geostat <-  transmute(Data,
                           Catch_KG = wCPUE*100,
                           Year = YEAR,
                           Vessel = "missing",
                           AreaSwept_km2 = 1,
                           Lat = START_LATITUDE,
                           Lon = START_LONGITUDE,
                           Pass = 0
)

strata.limits <- data.frame(STRATA = as.factor('All_areas'))

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x=50, Region="User", purpose="index", fine_scale=T, ObsModel= c(2,1),
                          strata.limits=strata.limits, bias.correct=F )

EBSgrid <- read.csv(file="F:/R/VAST2019/EBSThorsonGrid.csv")
input_grid=cbind(Lat=EBSgrid$Lat,Lon=EBSgrid$Lon,Area_km2=EBSgrid$Shape_Area/1000000) 

# Run model
fit = fit_model( "settings"=settings, "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'], 
                 "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=Data_Geostat[,'Vessel'],
                 "input_grid"=input_grid, max_cells=2000, test_fit=F)

# Plot results
plot( fit )
