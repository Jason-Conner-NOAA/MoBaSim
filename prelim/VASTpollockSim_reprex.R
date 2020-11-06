# Load packages
library(tidyverse)
library(TMB)               
library(VAST)
library(sf)
library(raster)
library(googledrive)


# Get data from Google ----------------------------------------------------

# this should fire up browser window
drive_auth()

dir.create(here::here("data"))

# specify url of folder you want to pull in
data <- drive_ls(as_id("https://drive.google.com/drive/folders/1ZyjInG4AaJetqsdVw4zSm9GmyNwDk4sh"))


# loop over files in folder if need to download many (or can do by file name, etc)
for(i in 1:nrow(data)) {
  drive_download(file=data[i,], 
                 path=here::here("data",data$name[i]), 
                 overwrite = TRUE)
}


# Set Directories ---------------------------------------------------------
Date = Sys.Date()
Run = 1
runDir = paste0('VAST',Date,'_V1')
dir.create(here::here(runDir))
setwd(here::here(runDir))


# Load Data ---------------------------------------------------------------
  load(here::here("data","EBSbundle.rdata"))
  #coldPool <- read.csv(here::here("data","cpa_areas2019.csv"))


# VAST settings -----------------------------------------------------------

species <- 21740
Version = get_latest_version( package="VAST" )
Region= "User"
n_x = 376   # Specify number of stations (a.k.a. "knots")
FieldConfig = matrix( c(0,0,"IID",0, 
                        1,1,"IID","Identity"), 
                      ncol=2, nrow=4, 
                      dimnames=list(c("Omega","Epsilon","Beta","Epsilon_year"),
                                    c("Component_1","Component_2"))
)
RhoConfig = c("Beta1"=4, "Beta2"=4, "Epsilon1"=0, "Epsilon2"=4)
OverdispersionConfig = c("Eta1"=0, "Eta2"=0)
ObsModel = c(10,2) #c(2,1)
Options =  c("Calculate_Range"=FALSE, "Calculate_effective_area"=FALSE, "treat_nonencounter_as_zero"=TRUE )
Aniso = FALSE
Npool = 100
BiasCorr = TRUE
Purpose = "index2"
FineScale=TRUE
TestFit = FALSE
MaxCells = 2000


# Data Formatting ---------------------------------------------------------

Data <- EBSbundle$EBSspp %>%
  dplyr::filter(SPECIES_CODE==species)

Data_Geostat <-  transmute(Data,
                           Catch_KG = wCPUE*100,
                           Year = YEAR,
                           Vessel = "missing",
                           AreaSwept_km2 = 1,
                           Lat = LAT_DEGREES,
                           Lon = LON_DEGREES,
                           Pass = 0
) %>%
  data.frame()


# VAST Spatial Setup ------------------------------------------------------

strata.limits <- data.frame(STRATA = 'All_areas')

predictSF <- EBSbundle$EBSpredict %>%
  mutate(LAT_m= LAT, LON_m=LONG) %>%
  st_as_sf(coords=c("LONG","LAT"), crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") %>%
  st_transform(crs = "+proj=longlat +datum=NAD83 +no_defs") %>%
  bind_cols(LAT_degrees = st_coordinates(.)[,2], LON_degrees = st_coordinates(.)[,1]) %>%
  st_drop_geometry() %>%
  data.frame()


input_grid=cbind(Lat=predictSF$LAT_degrees,Lon=predictSF$LON_degrees,Area_km2=predictSF$Area_m/1000000) 


# VAST Make Settings ------------------------------------------------------

settings = make_settings( Version = Version,
                          n_x= n_x, 
                          Region= Region, 
                          purpose= Purpose, 
                          Options =  Options,
                          fine_scale= FineScale, 
                          ObsModel= ObsModel,
                          FieldConfig = FieldConfig,
                          RhoConfig = RhoConfig,
                          OverdispersionConfig = OverdispersionConfig,
                          strata.limits=strata.limits, 
                          bias.correct=BiasCorr )

# Cold Pool covariate config ----------------------------------------------

# covariate_data = coldPool[ which(coldPool[,'YEAR'] %in% unique(Data_Geostat[,'Year'])), ]
# covariate_data = data.frame( "Year"=covariate_data[,"YEAR"], "Lat"=mean(Data_Geostat[,'Lat']), "Lon"=mean(Data_Geostat[,'Lat']), apply(covariate_data[-1],MARGIN=2,FUN=function(vec){(vec-mean(vec))/sd(vec)}) )

## Load covariates
# formula = ~ AREA_SUM_KM2_LTE2
# Xconfig_zcp = array(2, dim=c(2,1,1) )

# Run model
fit = fit_model( "settings"=settings, 
                 "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], 
                 "t_i"=Data_Geostat[,'Year'], 
                 "c_i"=rep(0,nrow(Data_Geostat)),
                 "b_i"=Data_Geostat[,'Catch_KG'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], 
                 "v_i"=Data_Geostat[,'Vessel'],
                 "input_grid"=input_grid,
                # knot_method = Method, 
                # Npool = Npool,
                 # "formula"=formula, 
                 # "covariate_data"=covariate_data, 
                 # "Xconfig_zcp"=Xconfig_zcp, 
                 test_fit = TestFit,
                 Aniso = Aniso, 
                 max_cells = MaxCells)

saveRDS(fit, "VASTfit.rds")

# Plot results
plot( fit )



# # Change some values to demonstrate capacity to change operating model
# Par = fit$tmb_list$Obj$env$last.par.best
# # Double decorrelation distance
# Par[c("logkappa1","logkappa2")] = Par[c("logkappa1","logkappa2")] - log(2)

# Loop through OM
for( rI in 1:1 ){
  Keep = FALSE
  while( Keep==FALSE ){
    Data_sim = fit$tmb_list$Obj$simulate( par=Par, complete=TRUE )
    Enc_t = tapply( Data_sim$b_i, INDEX=fit$data_frame$t_i, FUN=function(vec){mean(vec>0)})
    if( all(Enc_t>0 & Enc_t<1) ) Keep = TRUE
  }
  save(Data_sim, file=paste0(DateFile,"Data_sim",rI,".RData"))
}

# Loop through EM
for( rI in 1:Nrep ){
  load(file=paste0(DateFile,"Data_sim",rI,".RData"))
  RhoConfig = c("Beta1"=0, "Beta2"=0, "Epsilon1"=0, "Epsilon2"=0)
  ObsModel = c(2,0)
  
  Data_new = make_data("b_i"=Data_sim$b_i, "Version"=Version, "FieldConfig"=FieldConfig, "OverdispersionConfig"=OverdispersionConfig, "RhoConfig"=RhoConfig, "ObsModel"=ObsModel, "c_i"=rep(0,nrow(Data_Geostat)), "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year'], "spatial_list"=Spatial_List, "Options"=Options )
  TmbList_new = make_model("TmbData"=Data_new, "RunDir"=DateFile, "Version"=Version, "RhoConfig"=RhoConfig, "loc_x"=Spatial_List$loc_x, "Method"=Method)
  Obj_new = TmbList_new[["Obj"]]
  
  Opt_new = TMBhelper::fit_tmb( obj=Obj_new, lower=TmbList_new[["Lower"]], upper=TmbList_new[["Upper"]], getsd=TRUE, savedir=DateFile, bias.correct=TRUE, newtonsteps=1, bias.correct.control=list(sd=FALSE, split=NULL, nsplit=1, vars_to_correct="Index_cyl") )
  Report_new = Obj_new$report()
  
  Index = plot_biomass_index( DirName=DateFile, TmbData=Data_new, Sdreport=Opt_new[["SD"]], Year_Set=Year_Set, Years2Include=Years2Include, use_biascorr=TRUE )
  Save_new = list("Opt"=Opt_new, "Report"=Report_new, "ParHat"=Obj_new$env$parList(Opt_new$par), "Data"=Data_new, "Index"=Index)
  save(Save_new, file=paste0(DateFile,"Save_new_",rI,".RData"))
}

# Compile results
Index_array = array(NA, dim=c(Nrep,length(Year_Set),3,2), dimnames=list(paste0("Rep_",1:Nrep),Year_Set,c("Orig","True","Est"),c("Index","SE")) )
for( rI in 1:Nrep ){
  load(file=paste0(DateFile,"Save_orig.RData"))
  load(file=paste0(DateFile,"Data_sim",rI,".RData"))
  load(file=paste0(DateFile,"Save_new_",rI,".RData"))
  Index_array[rI,,'Orig','Index'] = Save_orig$Report$Index_cyl[1,,1]
  Index_array[rI,,'True','Index'] = Data_sim$Index_cyl[1,,1]
  Index_array[rI,,'Est','Index'] = Save_new$Index$Table[,'Estimate_metric_tons'] # Save_new$Report$Index_cyl[1,,1]
}

# Compile data frame of results
DF = cbind( expand.grid(dimnames(Index_array[,,'True','Index'])), "Orig"=as.vector(Index_array[,,'Orig','Index']), "True"=as.vector(Index_array[,,'True','Index']), "Est"=as.vector(Index_array[,,'Est','Index']) )

# Test hyper-stability, should be 1.00
Lm = lm( log(Est) ~ 0 + factor(Var1) + log(True), data=DF )
summary(Lm)$coef['log(True)',]
# Test bias in average, should be 1.00
Lm = lm( Est/True ~ 1, data=DF )
summary(Lm)$coef['(Intercept)',]


# Scrap -------------------------------------------------------------------

ggplot() +
  geom_point(data=EBSbundle$EBSpredict,
             aes(x=LONG,y=LAT,color=GEAR_TEMPERATURE1995))
