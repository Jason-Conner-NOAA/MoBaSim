#' Simulate Empirical Observation Model 
#' 
#' * Generates IDWs for probability of presence, CPUE, and the local variance (n=5) to the prediction grid
#' * Generates simulated species distributions by sampling IDWs for 1) presence, 2) if present, sample from the normal distribution
#'   with mean = predicted CPUE and sd = local standard deviation, 3) if sampled CPUE is negative, set it to 0

library(tidyverse)
library(sf)
library(sp)
library(FNN)
library(gstat) 


# Set Species code --------------------------------------------------------
  # 21740=pollock, 21720=cod, 10210=yellowfin, 10110=arrowtooth
  speciesCode <- 21720
  
  
# Load Data ---------------------------------------------------------------
  load(here::here("data","EBSbundle_1_2.rdata"))
  source(here::here("R","make_IDW.R"))
  
  grid_index <- readRDS(here::here("data","grid_index.rds"))
  crsString <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
  
  predictGrid <- EBSbundle$EBSfullPredict %>%
    inner_join(grid_index, by = "predict_id") %>% 
    mutate(LAT_m= LAT, LON_m=LONG) %>%
    st_as_sf(coords=c("LONG","LAT"), crs=crsString) 
  
  # Make spatial objects of survey data
  speciesData <- filter(EBSbundle$EBSspp, SPECIES_CODE==speciesCode) %>%
    st_as_sf(coords = c("LON_DEGREES","LAT_DEGREES"),
             crs = st_crs("+proj=longlat +datum=NAD83 +no_defs"),
             remove = F)  %>%
    st_transform(crs = crsString) 
    
  IDW_all <- list()
  years <- unique(speciesData$YEAR)


# Generate IDW predictions ------------------------------------------------
  for (yr in seq_along(years)) {
    ## IDW for CPUE
      sp_idw <- make_IDW(data=speciesData, 
                         var_name="wCPUE", 
                         year=years[yr], 
                         input_grid=predictGrid, 
                         idp = 6, 
                         nmax = 8) %>%
        dplyr::select(predict_id, grid_id, subgrid_id, YEAR, predictCPUE = var1.pred, BOTTOM_DEPTH, GEAR_TEMPERATURE = paste0("GEAR_TEMPERATURE",years[yr]), Area_m)
  
    
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
      
      sp_sd_idw <- make_IDW(data=sp_yr_sd, 
                            var_name="nn4sd", 
                            year=years[yr], 
                            input_grid=predictGrid, 
                            idp = 6, 
                            nmax = 8) %>%
        st_drop_geometry() %>%
        dplyr::select(predict_id, predictSD = var1.pred)
  
      
    ## Calculate IDW for presence 
      presenceData <- dplyr::mutate(speciesData, presence = ifelse(wCPUE==0,0,1))
      pres_idw <- make_IDW(data=presenceData, 
                           var_name="presence", 
                           year=years[yr], 
                           input_grid=predictGrid, 
                           idp = 6, 
                           nmax = 8) %>%
        st_drop_geometry() %>%
        dplyr::select(predict_id, predictPresence = var1.pred)
  
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
sim_dist <- lapply(IDW_all, 
                   FUN= function(i) {
                       sim <- i %>%
                         dplyr::mutate(sim_presence = ifelse(predictPresence >= runif(1),1,0), 
                                     sim_calc = rnorm(nrow(i), mean = predictCPUE, sd = predictSD),
                                     simCPUE = sim_presence * ifelse(sim_calc > 0, sim_calc, 0 )
                       )
                     
                   }
)
names(sim_dist) <- years

saveRDS(sim_dist, file = here::here("EOM",paste0("SimulatedDistribution_",speciesCode,"_",Sys.Date(),".RDS")))


# Simulate Random Surveys --------------------------------------------------------

survey_set = lapply(sim_dist, 
                    FUN= function(i) {
                      sim <- sample(i$simCPUE,350, replace = F)
                      
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








  