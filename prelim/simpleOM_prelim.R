library(tidyverse)
library(sf)
library(sp)
library(FNN)
library(gstat) 
#library(raster)


# Load Data ---------------------------------------------------------------
load(here::here("data","EBSbundle_1_2.rdata"))
source(here::here("R","make_IDW.R"))
#coldPool <- read.csv(here::here("data","cpa_areas2019.csv"))

crsString <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
speciesCode <- 21720

predictGrid <- EBSbundle$EBSfullPredict %>%
  mutate(LAT_m= LAT, LON_m=LONG) %>%
  st_as_sf(coords=c("LONG","LAT"), crs=crsString) 

# Make spatial objects of survey data
speciesData <- filter(EBSbundle$EBSspp, SPECIES_CODE==speciesCode) %>%
  st_as_sf(coords = c("LON_DEGREES","LAT_DEGREES"),
           crs = st_crs("+proj=longlat +datum=NAD83 +no_defs"),
           remove = F)  %>%
  st_transform(crs = crsString) 
  
# plot(st_geometry(predictGrid))
# plot(st_geometry(speciesData), add=T)


IDW_all <- list()

years <- unique(speciesData$YEAR)

for (yr in seq_along(years)) {
  ## IDW for CPUE
    sp_idw <- make_IDW(data=speciesData, var_name="wCPUE", year=years[yr], input_grid=predictGrid) %>%
      dplyr::select(predict_id, YEAR, predictCPUE = var1.pred, BOTTOM_DEPTH, GEAR_TEMPERATURE = paste0("GEAR_TEMPERATURE",years[yr]), Area_m)
    
    # IDW_cpue[[yr]] <- sp_idw
  
  ## Calculate local variance - then IDW 1 standard deviation for input into rnorm() during simulation
    sp_yr <- dplyr::filter(speciesData, YEAR== years[yr]) %>%
      st_drop_geometry()
    nn_4 <- FNN::get.knn(dplyr::select(sp_yr,LONG,LAT),k=4)
  
    sp_yr_sd <- data.frame()
    
    for (i in seq_along(sp_yr$CATCHJOIN)) {
      
      nn_4_sd = sqrt(var(c(sp_yr[i,"wCPUE"],sp_yr[nn_4$nn.index[i,],"wCPUE"])))   # multiply by 2 for 2 sd
      sp_yr_sd <- bind_rows(sp_yr_sd, data.frame(nn4sd = nn_4_sd))
    }
    
    sp_yr_sd <- cbind(dplyr::filter(speciesData, YEAR== years[yr]), sp_yr_sd)
    
    sp_sd_idw <- make_IDW(data=sp_yr_sd, var_name="nn4sd", year=years[yr], input_grid=predictGrid) %>%
      st_drop_geometry() %>%
      dplyr::select(predict_id, predictSD = var1.pred)
    
    # IDW_sd[[yr]] <- sp_sd_idw
    
  ## Calculate IDW for presence 
    presenceData <- dplyr::mutate(speciesData, presence = ifelse(wCPUE==0,0,1))
    pres_idw <- make_IDW(data=presenceData, var_name="presence", year=years[yr], input_grid=predictGrid) %>%
      st_drop_geometry() %>%
      dplyr::select(predict_id, predictPresence = var1.pred)
    
    # IDW_presence[[yr]] <- pres_idw

  ## Stack IDWs
    sp_all_idw <- merge(sp_idw, sp_sd_idw, by="predict_id") %>%
      merge(pres_idw, by="predict_id") %>%
      dplyr::filter(predict_id %in% EBSbundle$EBSpredict$predict_id)
    IDW_all[[yr]] <- sp_all_idw
}

names(IDW_all) <- years

saveRDS(IDW_all, file = here::here("EOM",paste0("IDW_all_",speciesCode,".RDS")))

# IDW Testing -------------------------------------------------------------
t = 5 # Test year index
test_cpue <- st_join(speciesData[speciesData$YEAR==years[t],],IDW_all[[t]], join = st_nearest_feature ) %>%
  mutate(cpueDiff = wCPUE - predictCPUE)
hist(test_cpue$cpueDiff)

qqplot(speciesData[speciesData$YEAR==years[t],]$wCPUE, IDW_all[[t]]["BOTTOM_DEPTH">0,]$predictCPUE)
abline(0,1)



# Create a set of simulated populations -----------------------------------

sim_dist <- lapply(years, 
                   FUN= function(y) {
                     yearSet <- IDW_all[[as.character(y)]]
                     sim <- yearSet %>%
                       dplyr::mutate(sim_presence = ifelse(predictPresence >= runif(1),1,0), 
                                     sim_calc = rnorm(nrow(yearSet), mean = predictCPUE, sd = predictSD),
                                     simCPUE = sim_presence * ifelse(sim_calc > 0, sim_calc, 0 )
                                     )
                            
                   }
)
names(sim_dist) <- years

saveRDS(sim_dist, file = here::here("EOM",paste0("SimulatedDistribution_",speciesCode,"_",Sys.Date(),".RDS")))

# Simulate Random Surveys --------------------------------------------------------

survey_set = lapply(years, 
                    FUN= function(y) {
                      yearSet <- sim_dist[[as.character(y)]]
                      sim <- sample(yearSet$simCPUE,350, replace = F)
                      
                    }
)

names(survey_set) <- years


# Test simulated surveys --------------------------------------------------

meanComp <- sapply(years, FUN= function(y) {
  trueMean <- mean(sim_dist[[as.character(y)]]$simCPUE)
  simMean <- mean(survey_set[[as.character(y)]])
  diff <- trueMean - simMean
}
)

hist(meanComp)








  