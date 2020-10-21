# Download release number 3.0.0; its useful for reproducibility to use a specific release number

RootDir = file.path(getwd(),"PollockSim")
Date = Sys.Date()
DateFile = paste0(RootDir,'/',Date,'_V1/')
dir.create(DateFile)

# Load packages
library(sumfish)
library(TMB)               
library(VAST)


# load data set
# see `?load_example` for list of stocks with example data 
# that are installed automatically with `FishStatsUtils`. 
species <- 21740
#racebase <- sumfish::getRacebase(1982:2019,'EBS_SHELF')

Data <- readRDS(paste0(getwd(),"/haulSum82_19.RDS")) %>%#  %>% ## Hardcoded from saved haulSum sumHaul(racebase) 
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
settings = make_settings( n_x=100, Region="User", purpose="index", fine_scale=T, ObsModel= c(2,1),
                          strata.limits=strata.limits, bias.correct=F )

EBSgrid <- read.csv(file="F:/R/VAST2019/EBSThorsonGrid.csv")
input_grid=cbind(Lat=EBSgrid$Lat,Lon=EBSgrid$Lon,Area_km2=EBSgrid$Shape_Area/1000000) 

# Run model
fit = fit_model( "settings"=settings, "Lat_i"=Data_Geostat[,'Lat'], 
                 "Lon_i"=Data_Geostat[,'Lon'], "t_i"=Data_Geostat[,'Year'], 
                 "c_i"=rep(0,nrow(Data_Geostat)), "b_i"=Data_Geostat[,'Catch_KG'], 
                 "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=Data_Geostat[,'Vessel'],
                 "input_grid"=input_grid)

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