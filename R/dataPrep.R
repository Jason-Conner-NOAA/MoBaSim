library(sumfish)
library(sf)
library(gstat)
library(raster)
library(FNN)
library(PBSmapping)
library(automap)


load("F:/R/simfish/data/Predict_data.rda")
EBSstrata <- st_read(dsn="F:/Research/Sampling Density Simulation/R/Seed401/shapefiles", layer="EBSstrata")

EBS <- getRacebase(1995:2019,'EBS_SHELF') #saveRDS(EBS,"F:/R/sumfish_EBS.RDS")
EBSdata <- sumHaul(EBS)
saveRDS(EBS,"F:/R/sumfish_EBS.RDS")


EBSsppSF <- EBSdata %>%
  filter(SPECIES_CODE %in% c(21740,21720,10210,10110)) %>%
  inner_join(EBS$species, by="SPECIES_CODE") %>%
  mutate(LAT_DEGREES= START_LATITUDE, LON_DEGREES=START_LONGITUDE) %>%
  st_as_sf(coords=c("START_LONGITUDE","START_LATITUDE"), crs = crs("+proj=longlat +datum=NAD83 +no_defs")) %>%
  st_transform(crs = crs(EBSstrata))

EBSspp <- EBSsppSF %>%
  st_drop_geometry() %>%
  data.frame() %>%
  bind_cols(data.frame(st_coordinates(EBSsppSF))) %>%
  rename(LONG=X, LAT=Y)
  
  

# unique(EBSspp$SPECIES_NAME)

# Reduce prediction grid size by a factor of 4
predictRaster <- filter(Predict_data, EBS==1)  %>%
  dplyr::select(LONG,LAT,BOTTOM_DEPTH) %>%
  rasterFromXYZ() %>%
  raster::aggregate(fact=4)

crs(predictRaster) = crs(EBSstrata)

predictArea <- area(predictRaster)

predictLatLong <- predictRaster
projection(predictLatLong) <- CRS("+proj=longlat +datum=NAD83 +no_defs")


### From Kot's code - predict gear temperature
### Now include the SST values based on surface kriging of available temperature data (BS+GOA)


dat_new <- dplyr::filter(EBSspp, SPECIES_CODE == 21740) %>%
  dplyr::select("YEAR", "GEAR_TEMPERATURE",  "HAULJOIN","LATITUDE"=LAT_DEGREES, "LONGITUDE"=LON_DEGREES)

dat_new <- dat_new[!duplicated(dat_new),]
coordinates(dat_new) <- ~ LONGITUDE + LATITUDE
crs(dat_new) = "+init=epsg:4269"
dat_new <- sp::spTransform(dat_new, crs(EBSstrata))

YEARS = sort(unique(dat_new$YEAR))
raster_temp_bottom <- list()
missingTemp <- data.frame()

for (yr in seq_along(YEARS))
{

  # For gear temperature
  dat1 <- dat2 <- subset(dat_new, subset=c(YEAR == YEARS[yr]))
  dat1 <- dat1[!is.na(dat1$GEAR_TEMPERATURE),]
  m <- autofitVariogram(GEAR_TEMPERATURE~1, dat1)
  # plot(m)
  
  v <- variogram(GEAR_TEMPERATURE~1, dat1)	#create a variogram of the sorting data
  m <- fit.variogram(v, vgm(psill=m$var_model[2,2], model=as.character(m$var_model[2,1]), range=m$var_model[2,3], nugget=m$var_model[1,2], kappa=m$var_model[2,4]))    #fit a model to the variogram
  # plot(v, model= m)
  
  bathy <- predictRaster
  g <- gstat(id = "GEAR_TEMPERATURE", formula = GEAR_TEMPERATURE~1, data=dat1, model = m, nmax=5)
  bottomInt <- raster::interpolate(bathy, g, xyOnly=TRUE, xyNames=c('LONG','LAT'), progress="text", overwrite=TRUE) #Interpolate the object to a raster
  
  # Replace missing gear_temperatures 
  checkTemp <- raster::extract(bottomInt, dat2@coords) %>%
    data.frame(krigeTemp=.) %>%
    bind_cols(dat2@data, .) %>%
    mutate(diff = GEAR_TEMPERATURE - krigeTemp,
           fillTemp = ifelse(is.na(GEAR_TEMPERATURE),krigeTemp,GEAR_TEMPERATURE)
                                     ) %>%
    dplyr::select(HAULJOIN, diff, fillTemp)
  
  missingTemp <- bind_rows(missingTemp, checkTemp)
  
  
  names(bottomInt) <- paste0('GEAR_TEMPERATURE',YEARS[yr])
  
  raster_temp_bottom[[yr]] <- bottomInt

}


bottomTemps <- append(raster_temp_bottom,predictRaster,0) %>%
  raster::stack() %>%
  append(predictArea,0) %>%
  raster::stack() %>%
  as.data.frame(xy=TRUE)


  

EBSspp <- EBSspp %>%
  inner_join(missingTemp, by = "HAULJOIN") %>%
  mutate( GEAR_TEMPERATURE = ifelse(is.na(GEAR_TEMPERATURE), fillTemp, GEAR_TEMPERATURE)) %>%
  dplyr::select(-fillTemp, -diff)

# # Loop though temperature data and get IDW for each year
# # y = 2001
# 
# EBSframe <- raster::rasterToPoints(predictRaster) %>%
#   data.frame() %>%
#   rename(LONG=1, LAT=2) %>%
#   dplyr::select(LONG,LAT,BOTTOM_DEPTH)
# 
# test <- list()
# 
# for (y in unique(EBSdata$YEAR)) {
#   EBSyear <- dplyr::filter(EBS$haul, CRUISEJOIN %in% EBS$cruise[EBS$cruise$YEAR==y,"CRUISEJOIN"] & !is.na(GEAR_TEMPERATURE)) %>%
#     sf::st_as_sf(coords = c(x = "LON_DEGREES", y = "LAT_DEGREES"), crs = sf::st_crs(4269)) %>% 
#     sf::as_Spatial() %>%
#     sp::spTransform(crs(EBSstrata))
#   
#   colName <- paste0('GEAR_TEMPERATURE',y)
#   
# ## From akgfmaps
#     # Inverse distance weighting----------------------------------------------------------------------
#     idw_fit <- gstat::gstat(formula = GEAR_TEMPERATURE~1, locations = EBSyear, nmax = 4)
#     
#     # Predict station points--------------------------------------------------------------------------
#     #stn.predict <- predict(idw_fit, EBSyear)
#     
#     extrap.grid <- predict(idw_fit, as(predictRaster, "SpatialPoints")) 
#   
#   
#   # test prediction
#   test <- get.knnx(extrap.grid@coords, EBSyear@coords, k=1)
#   
#   test2 <- bind_cols(data.frame(EBSyear), data.frame(extrap.idx =as.character(test[[1]]))) %>%
#     inner_join(rownames_to_column(data.frame(extrap.grid), var= "extrap.idx"), by = c("extrap.idx")) %>%
#     transmute(deltaTemp = GEAR_TEMPERATURE - var1.pred)
#   
#   yearData <- data.frame(extrap.grid@data$var1.pred) %>%
#     rename(!!colName := 1)
#   
#   
#   
#   EBSframe <- bind_cols(EBSframe, yearData)
#      
# }
# 
#   
# # testRow <- 3000
# # extrap.grid@coords[testRow,]
# # EBSpredict[testRow,]

EBSfullPredict <- bottomTemps %>%
  rename(LONG=x, LAT=y, Area_m = layer) %>%
  mutate(predict_id = row_number())

# eliminate depths outside of survey range
EBSpredict <- EBSfullPredict %>%
  dplyr::filter(BOTTOM_DEPTH >= 20) %>%
  dplyr::filter(BOTTOM_DEPTH <=200)

# Verify no NAs in prediction dataframe
naTest <- EBSpredict[is.na(EBSpredict)]



EBSpredictSF <- EBSpredict %>%
  mutate(LON2 = LONG, LAT2 = LAT) %>%
  st_as_sf(coords = c('LON2','LAT2'),
           crs=st_crs(EBSstrata)) %>%
  st_join(EBSstrata,.)

EBSpredictSP <- EBSpredict 
coordinates(EBSpredictSP) <- ~ LONG + LAT
projection(EBSpredictSP) <- projection(EBSstrata)


# Check Temperature Predictions -------------------------------------------





readme <- "EBSspp - catch, effort, and environmental data for the EBS survey, years 1995-2019 for 4 spp (cod, pollock, yellowfin,
            arrowtooth). Use wCPUE, units are kg/ha, and LON_DEGREES/LAT_DEGREES for station location.  EBSstrata - sf 
            package, polygons for the EBS sampling frame. EBSpredict prediction dataframe reduced by a factor of 4 from Kot's original grid
            EBSpredictSF- sf dataframe (geometry=sfc_POLYGON) with the prediction grid (including Kotaro's extrapolated depths and bottom 
            temperatures for each year. EBS - this is the raw RACE database for 1995-2019, saved for reproducibility."

EBSbundle <- list(readme = readme,
                  EBSspp = EBSspp,
                  EBSstrata = EBSstrata,
                  EBSfullPredict = EBSfullPredict,
                  EBSpredict = EBSpredict,
                  EBSpredictSF = EBSpredictSF,
                  EBSpredictSP = EBSpredictSP,
                  EBS = EBS
)
save(EBSbundle, file='EBSbundle.rdata')

df <- data.frame(st_drop_geometry(EBSbundle$EBSpredict))
write.table(df,"EBSpredict.txt", sep=',', row.names=F)




# Plots -------------------------------------------------------------------

ggplot() +
  geom_point(data=EBSpredict[which(!is.na(EBSpredict$BOTTOM_DEPTH)),], mapping = aes(x=LONG, y=LAT, color="red"))


ggplot(data=EBSpredict, aes(x=LONG,y=LAT)) +
  geom_point(aes(color = GEAR_TEMPERATURE2000))