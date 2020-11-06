library(tidyverse)
library(sf)
library(sp)
library(FNN)
library(gstat) 
library(raster)
library(fasterize)


# Load Data ---------------------------------------------------------------
load(here::here("data","EBSbundle.rdata"))
#coldPool <- read.csv(here::here("data","cpa_areas2019.csv"))



rb <- EBSbundle$EBS
haul <- EBSbundle$EBSspp
EBSpredict <- EBSbundle$EBSpredictSP 

### Kept this for transformation code
# predictGrid <- EBSbundle$EBSfullPredict %>%
#   mutate(LAT_m= LAT, LON_m=LONG) %>%
#   st_as_sf(coords=c("LONG","LAT"), crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
#   st_transform(crs = "+proj=longlat +datum=NAD83 +no_defs") %>%
#   bind_cols(LAT_degrees = st_coordinates(.)[,2], LON_degrees = st_coordinates(.)[,1]) %>%
#   st_transform(crs = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
#   st_drop_geometry() %>%
#   data.frame()
predictGrid <- EBSbundle$EBSfullPredict %>%
  mutate(LAT_m= LAT, LON_m=LONG) %>%
  st_as_sf(coords=c("LONG","LAT"), crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
  dplyr::select(LON_m,LAT_m,BOTTOM_DEPTH) %>%
  raster::rasterFromXYZ(crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")


# Make spatial objects of survey data
atf <- filter(haul, SPECIES_CODE==10110) %>%
  st_as_sf(coords = c("LON_DEGREES","LAT_DEGREES"),
           crs = st_crs("+proj=longlat +datum=NAD83 +no_defs"),
           remove = F)  %>%
  st_transform(crs = "+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") 
  
plot(predictGrid)
plot(st_geometry(atf), add=T)

IDW_result <- list()
years <- unique(atf$YEAR)

for (yr in seq_along(years)) {
  sp_yr <- dplyr::filter(atf, YEAR==years[yr])
  
  sp_idw <- gstat::idw(formula = sp_yr$wCPUE ~ 1, locations = sp_yr, idp = 8, newdata = predictGrid) %>%
    mutate(YEAR = years[yr]) %>%
    bind_cols(predictGrid)
  
  
  IDW_result[[yr]] <- sp_idw
}

test <- st_join(atf,sp_idw, join = st_nearest_feature ) %>%
  mutate(cpueDiff = wCPUE - var1.pred)

# Make a full prediction raster -------------------------------------------
  EBSraster <- rasterFromXYZ(cbind(predictGrid$LON_degrees, predictGrid$LAT_degrees))

# ggplot() +
#   geom_sf(data=predictSF, aes(color='blue')) +
#   geom_sf(data=atf, aes(fill=wCPUE))

# cod <- filter(haul, SPECIES_CODE==21720) %>%
#   st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
#            crs = 4326,
#            remove = F) %>%
#   st_transform(crs = 102006)
# 
# pol <- filter(haul, SPECIES_CODE==21740) %>%
#   st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
#            crs = 4326,
#            remove = F) %>%
#   st_transform(crs = 102006)
# 
# yfs <- filter(haul, SPECIES_CODE==10210) %>%
#   st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
#            crs = 4326,
#            remove = F) %>%
#   st_transform(crs = 102006)


years <- sort(unique(EBSbundle$EBSspp$YEAR))

# Run for each species to get local variance of 4 NN
species <- 'Arrowtooh'
spp <- atf
sppSim <- st_sf(st_sfc(),crs=102006)

for (y in years) {
  # subset by year
  ySpp <- filter(spp, YEAR==y) %>%
    st_drop_geometry()
  
  # Get local variance
  nn <- get.knn(select(ySpp,START_LONGITUDE,START_LATITUDE),k=4)
  
  for (i in 1:nrow(ySpp)) {
    # TEST FOR LOCATION of NN
    # t =rbind(ySpp[i,],ySpp[nn$nn.index[i,],])
    # plot(t$START_LONGITUDE,t$START_LATITUDE, pch=c('A','B','C','D','E'))
    
    # Calculate Var 
    nn.var <- var(c(ySpp[i,"wCPUE"],ySpp[nn$nn.index[i,],"wCPUE"]))
    
    ySpp[i,"wCPUEvar"] <- nn.var

  }
  
  ySpp <- mutate(ySpp, presence = ifelse(wCPUE==0,0,1)) %>%
    st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
             crs = 102006,
             remove = F) 
  
  sppSim <- rbind(sppSim, ySpp)
}


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


  