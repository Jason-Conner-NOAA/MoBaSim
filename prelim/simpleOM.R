library(sumfish)
library(sf)
library(sp)
library(FNN)
library(gstat) 
library(raster)

getSQL()

# Load data framework from F:/Research/Model Based Sim/MoBaSim/dataPrep generated 2/28/2020
load("F:/Research/Model Based Sim/MoBaSim/EBSbundle.rdata")

rb <- EBSbundle$EBS
haul <- EBSbundle$EBSspp
EBSpredict <- EBSbundle$EBSpredictSP 

EBSpredict <- EBSbundle$EBSpredictSF %>%
  st_drop_geometry() %>%
  st_as_sf(coords = c('LONG','LAT'),
           crs = 102006, # AK Albers EA
           remove = F)

# Make spatial objects of survey data
atf <- filter(haul, SPECIES_CODE==10110) %>%
  st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
           crs = 4326,
           remove = F) %>%
  st_transform(crs = 102006)

cod <- filter(haul, SPECIES_CODE==21720) %>%
  st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
           crs = 4326,
           remove = F) %>%
  st_transform(crs = 102006)

pol <- filter(haul, SPECIES_CODE==21740) %>%
  st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
           crs = 4326,
           remove = F) %>%
  st_transform(crs = 102006)

yfs <- filter(haul, SPECIES_CODE==10210) %>%
  st_as_sf(coords = c("START_LONGITUDE","START_LATITUDE"),
           crs = 4326,
           remove = F) %>%
  st_transform(crs = 102006)


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


  