library(tidyverse)
library(sf)
library(sp)
library(FNN)
library(gstat) 
library(raster)


# Load Data ---------------------------------------------------------------
load(here::here("data","EBSbundle.rdata"))
#coldPool <- read.csv(here::here("data","cpa_areas2019.csv"))

crsString <- "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
speciesCode <- 10110

predictGrid <- EBSbundle$EBSfullPredict %>%
  mutate(LAT_m= LAT, LON_m=LONG) %>%
  st_as_sf(coords=c("LONG","LAT"), crs=crsString) 

# Make spatial objects of survey data
speciesData <- filter(EBSbundle$EBSspp, SPECIES_CODE==speciesCode) %>%
  st_as_sf(coords = c("LON_DEGREES","LAT_DEGREES"),
           crs = st_crs("+proj=longlat +datum=NAD83 +no_defs"),
           remove = F)  %>%
  st_transform(crs = crsString) 
  
# plot(predictGrid)
# plot(st_geometry(atf), add=T)

IDW_cpue <- list()
IDW_sd <- list()
IDW_presence <- list()
IDW_all <- data.frame()

years <- unique(speciesData$YEAR)

for (yr in seq_along(years)) {
  ## IDW for CPUE
    sp_idw <- make_IDW(data=speciesData, var_name="wCPUE", year=years[yr], input_grid=predictGrid)
    
    IDW_cpue[[yr]] <- sp_idw
  
  ## Calculate local variance - then IDW 2 standard deviations
    sp_yr <- dplyr::filter(speciesData, YEAR== years[yr]) %>%
      st_drop_geometry()
    nn_4 <- FNN::get.knn(dplyr::select(sp_yr,LONG,LAT),k=4)
  
    sp_yr_sd <- data.frame()
    
    for (i in seq_along(sp_yr$CATCHJOIN)) {
      
      nn_4_sd = 2*sqrt(var(c(sp_yr[i,"wCPUE"],sp_yr[nn_4$nn.index[i,],"wCPUE"])))
      sp_yr_sd <- bind_rows(sp_yr_sd, data.frame(nn4sd = nn_4_sd))
    }
    
    sp_yr_sd <- cbind(dplyr::filter(speciesData, YEAR== years[yr]), sp_yr_sd)
    
    sp_sd_idw <- make_IDW(data=sp_yr_sd, var_name="nn4sd", year=years[yr], input_grid=predictGrid)
    
    IDW_sd[[yr]] <- sp_sd_idw
    
  ## Calculate IDW for presence 
    presenceData <- dplyr::mutate(speciesData, presence = ifelse(wCPUE==0,0,1))
    pres_idw <- make_IDW(data=presenceData, var_name="presence", year=years[yr], input_grid=predictGrid)
    
    IDW_presence[[yr]] <- pres_idw
    
}

test_cpue <- st_join(speciesData[speciesData$YEAR==years[3],],IDW_cpue[[3]], join = st_nearest_feature ) %>%
  mutate(cpueDiff = wCPUE - var1.pred)
test_var <- st_join(speciesData[speciesData$YEAR==years[3],],IDW_var[[3]], join = st_nearest_feature ) 
#hauljoin -2922 has highest variance in 2007

qqplot(speciesData[speciesData$YEAR==years[22],]$wCPUE, IDW_cpue[[22]]$var1.pred)
abline(0,1)

# Run for each species to get local variance of 4 NN





# Function to take input raster (prediction grid) and spatially defined variables
idwSim <- function(year, sfData, var1="wCPUE", sfPredict, idp=2, plot=F) {
  yData <- filter(sfData, YEAR==year)
  
  simIDW <- idw(as.formula(paste0(var1,"~1")),yData,newdata=sfPredict, idp=idp) %>%
    st_transform(crs=st_crs(sfPredict))
  
  if (plot) {
    png(file = paste0(species,'_IDW_',var1,'_',year,'.png'), bg = "transparent")
    plot(st_geometry(simIDW), col = simIDW$var1.pred)
    plot(st_geometry(EBSbundle$EBSstrata), add=T)
    dev.off()
  }
  
  return(simIDW)
}

# TEST
tatf <- filter(atf, YEAR==2015) 
tIDW <- idwSim(2015,spp,"wCPUE",EBSpredict, plot=T)
plot(st_geometry(tIDW), col = tIDW$var1.pred)


cpueIDW <- lapply(years, FUN=idwSim, sppSim,"wCPUE",EBSpredict, plot=T)
names(cpueIDW) <- years

varIDW <- lapply(years, FUN=idwSim, sppSim,"wCPUEvar",EBSpredict, plot=T)
names(varIDW) <- years

binIDW <- lapply(years, FUN=idwSim, sppSim,"presence",EBSpredict, plot=T)
names(binIDW) <- years


  