## Script for estimating biomass survey indices in Eastern Bering Sea using the surveyIndex package.
## Also, a small simulation study is carried out conditional on the estimated fits.
## Author: Casper W. Berg, DTU Aqua.

library(surveyIndex) ## Make sure to use latest version (1.09)
library(rgdal)

## It is highly recommended to run this with an optimized BLAS library (e.g. openBLAS)
## and to use several cores (if available) in order to reduce running times.
if( suppressWarnings(require(RhpcBLASctl)) ){ blas_set_num_threads(4) } 

YEARS = 1995:2019
NYEARS = length(YEARS)

LongLat2EBSproj<-function(x,y){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")
    res <- spTransform(xy, CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
    return(as.data.frame(res))
}

EBSproj2LongLat<-function(x,y){
    xy <- data.frame(ID = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
    res <- spTransform(xy, CRS("+proj=longlat +datum=WGS84"))
    return(as.data.frame(res))
}

# load(here::here("data","EBSbundle_1_2.rdata"))
load(here::here("data","EBSbundle_truncate_5_years.rdata"))

###################
## Prepare data
###################

## transform prediction grid coordinates to longitude, latitude  
ll = EBSproj2LongLat(EBSbundle$EBSpredict$LONG,EBSbundle$EBSpredict$LAT)
EBSbundle$EBSpredict$lon = ll$X
EBSbundle$EBSpredict$lat = ll$Y
## center and scale EBSproj for numerical stability
mx = mean(EBSbundle$EBSpredict$LONG)
my = mean(EBSbundle$EBSpredict$LAT)
EBSbundle$EBSpredict$sx = (EBSbundle$EBSpredict$LONG - mx)/1000
EBSbundle$EBSpredict$sy = (EBSbundle$EBSpredict$LAT - my)/1000

d = EBSbundle$EBSspp

summary(d)

factors=c("YEAR","VESSEL","CRUISE","STRATUM","STATIONID","SPECIES_CODE","SPECIES_NAME","COMMON_NAME","SUBSAMPLE_CODE","REGION")

for(fac in factors) d[,fac]=as.factor(d[,fac])

## Add and rename some variables (e.g. names 'Year' and 'ctime' are needed for 'surveyIndex')
d$Year = d$YEAR
d$Ship = d$VESSEL
d$lat = d$START_LATITUDE
d$lon = d$START_LONGITUDE
## data lon,lat => grid projection
ll = LongLat2EBSproj( EBSbundle$EBSspp$LON_DEGREES,EBSbundle$EBSspp$LAT_DEGREES)
d$sx = (ll$X - mx)/1000
d$sy = (ll$Y - my)/1000
d$ctime <- as.numeric(as.character(d$Year))
d$hour = as.numeric(format(d$START_TIME,"%H"))
d$minute = as.numeric(format(d$START_TIME,"%M"))
d$day = as.numeric(format(d$START_TIME,"%d"))
d$month = as.numeric(format(d$START_TIME,"%m"))
d$TimeShotHour = d$hour + d$minute/60
d$timeOfYear <- (d$month - 1) * 1/12 + (d$day - 1)/365
d$Gear = "dummy"
d$Quarter = "2"

dall = d

## Data frame with one row per haul to DATRASraw alike object
df2dr<-function(x){
    x$haul.id = 1:nrow(x)
    dd = list()
    dd[[1]] = data.frame()
    dd[[2]] = x
    dd[[3]] = data.frame()
    class(dd)<-"DATRASraw"
    dd
}


## prediction a grid for a specific year 
getpd<-function(year,subsel=NULL){
    nam = c("lon","lat","sx","sy","BOTTOM_DEPTH")
    nam2 = paste0(c("GEAR_TEMPERATURE"),year)
            
    pd = data.frame(EBSbundle$EBSpredict[,c(nam,nam2)])
    colnames(pd)<-c(nam,c("GEAR_TEMPERATURE"))
    
    pd$EFFORT = 1.0 
    ## Note, Better to use median effort here if splines on effort is used.
    ## Otherwise uncertainty is inflated because it is outside normal effort range.
    
    if(!is.null(subsel))
        pd = pd[subsel,]
    pd
}

## get list of prediction grids by year
allpd = lapply(YEARS,getpd)
names(allpd) <- as.character(YEARS)

## split data by species, make into DATRASraw + Nage matrix
ds = split(dall,dall$COMMON_NAME)
ds = lapply(ds,df2dr)
## OBS, response is added here in "Nage" matrix -- use wCPUE
ds = lapply(ds,function(x) { x[[2]]$Nage <- matrix(x$wCPUE,ncol=1); colnames(x[[2]]$Nage)<-1; x } )
# dtest = lapply(ds,function(x) { x[[2]]$Nage <- matrix(x$wCPUE,ncol=1);x } )

# reduce to one species for testing JC
# ds = ds[[1]]

##############
## Fitting
##############
fm = "Year + s(sx,sy,bs=c('ts'),k=376)+s(sx,sy,bs=c('ts'),k=50,by=Year,id=1) +
      s(BOTTOM_DEPTH,bs='ts',k=10)+s(log(GEAR_TEMPERATURE+3),bs='ts',k=10)"

models = list()
fittimes = list()
specLevels = levels(dall$COMMON_NAME)

for(SPECIES in specLevels){
    cat("Fitting ",SPECIES,"\n")
        
    fittimes[[ SPECIES ]] <- system.time( models[[ SPECIES ]] <- getSurveyIdx(ds[[ SPECIES ]],1,myids=NULL,predD=allpd,cutOff=0,fam="Tweedie",modelP=fm,gamma=1,control=list(trace=TRUE,maxit=20)) )

}

## Check basis dimensions splines (spatial resolution)
sink("gamcheck.txt")
lapply(models,function(x) gam.check(x$pModels[[1]]))
sink()

## Model summaries
sink("summaries.txt")
lapply(models,function(x) summary(x$pModels[[1]]))
sink()

#################
## Plots
#################
pngscal = 2

for(SPECIES in specLevels){
    ddd = ds[[ SPECIES ]]

    png(paste0(SPECIES,"-%03d.png"),width=640*pngscal,height=480*pngscal)
    
    ## abundance maps
    surveyIdxPlots(models[[ SPECIES ]],ds[[ SPECIES ]],predD=allpd[[1]],select="absolutemap",year=YEARS,colors=rev(heat.colors(10)),
                   par=list(mfrow=n2mfrow(NYEARS),mar=c(0,0,2,0)),map.cex=1.3)
    
    ## spatial residuals
    par(mfrow=n2mfrow(NYEARS),mar=c(0,0,2,0))
    for(yy in YEARS) surveyIdxPlots(models[[ SPECIES ]],ds[[ SPECIES ]],myids=NULL,predD=allpd[[1]],select="spatialResiduals",year=yy,colors=rev(heat.colors(10)),par=list(),map.cex=1.5,legend=FALSE,axes=FALSE,main=yy)
    
    ## further residuals
    par(mfrow=c(3,3),mar=c(4,4,4,4),oma=c(2,2,4,2))
    resi = residuals(models[[ SPECIES ]],1)
    
    hist(resi,80,xlab="residuals",main="")
    title(paste0("Residuals - ",SPECIES),outer=TRUE)
    qq.gam(models[[ SPECIES ]]$pModels[[1]],rep=200,level=0.95)

    plot(factor( ds[[ SPECIES]]$Ship),resi,xlab="Vessel",ylab="residual")
    abline(h=0,col=2)

    plot(factor( ds[[ SPECIES]]$Year),resi,xlab="Year",ylab="residual")
    abline(h=0,col=2)
    
    plot(factor( paste(ds[[ SPECIES]]$Year,ds[[ SPECIES ]]$Ship,sep=":") ),resi,xlab="Year X vessel",ylab="residual")
    abline(h=0,col=2)

    ## Residuals vs. time of year
    plot(ds[[ SPECIES]]$timeOfYear,resi,xlab="Time of year",ylab="residual")
    oo = order(ds[[ SPECIES]]$timeOfYear)
    abline(h=0,col=2)
    lines(ds[[ SPECIES ]]$timeOfYear[oo],fitted(gam(resi~s(ddd$timeOfYear)))[oo],col=3,lwd=2)
    

    ## Residuals vs. time of day
    plot(ds[[ SPECIES]]$TimeShotHour,resi,xlab="Time of day",ylab="residual")
    oo = order(ds[[ SPECIES]]$TimeShotHour)
    abline(h=0,col=2)
    lines(ds[[ SPECIES ]]$TimeShotHour[oo],fitted(gam(resi~s(ddd$TimeShotHour,bs="cc")))[oo],col=3,lwd=2)
    
    ## Residuals vs. effort
    plot(ds[[ SPECIES ]]$EFFORT,resi,,xlab="EFFORT",ylab="residual")
    oo = order(ds[[ SPECIES]]$EFFORT)
    abline(h=0,col=2)
    lines(ds[[ SPECIES ]]$EFFORT[oo],fitted(gam(resi~s(ddd$EFFORT)))[oo],col=3,lwd=2)
  
    ## Depth effect
    par(mfrow=c(1,1),mar=c(4,3,3,3))
    plot.gam( models[[ SPECIES ]]$pModels[[1]] , select=NYEARS+2, residuals=TRUE,shade=TRUE,main=SPECIES)

    ## Temperature effect
    plot.gam( models[[ SPECIES ]]$pModels[[1]] , select=NYEARS+3, residuals=TRUE,shade=TRUE,main=SPECIES)    

    surveyIndex:::plot.SIlist( list(models[[ SPECIES ]] ) , main = SPECIES)
    surveyIndex:::plot.SIlist( list(models[[ SPECIES ]] ),rescale=TRUE, main = SPECIES)
    
    dev.off()
}

###############################
## Retrospective analysis 
###############################
doRETRO = FALSE  

if(doRETRO){
    retros = list()
    for(SPECIES in specLevels){
        cat("Doing retro for ",SPECIES,"\n")
        retros[[ SPECIES ]] = retro.surveyIdx(models[[SPECIES]],ds[[SPECIES]],grid=NULL,predD=allpd,npeels=3,control=list(trace=TRUE,maxit=10))
    }
    
    png("retros.png",width=640*pngscal,height=480*pngscal)
    par(mfrow=c(2,2))
    for(SPECIES in specLevels){
        plot(retros[[SPECIES]],base=models[[SPECIES]],main=SPECIES,lwd=2.5)
    }
    dev.off()
}

###############################
## Mini simulation study
###############################
REPS = 4
ests = list()

for(SPECIES in specLevels){
    ests[[ SPECIES ]] = list()
    
    ## simulate data
    csim = surveyIdx.simulate(models[[ SPECIES]],ds[[ SPECIES ]])
    sims = lapply(1:REPS,function(i) surveyIdx.simulate(models[[ SPECIES]],ds[[ SPECIES ]],FALSE,csim) )

    ## re-estimate
    tmp = ds[[ SPECIES ]]
    for(i in 1:REPS){
        tmp[[2]]$Nage = matrix(sims[[i]][[1]][,1],ncol=1)
        colnames(tmp$Nage)<-1
        
        ests[[SPECIES]][[i]]  <- getSurveyIdx(tmp,1,myids=NULL,predD=allpd,cutOff=0,fam="Tweedie",modelP=fm,gamma=1,control=list(trace=TRUE,maxit=10))
        cat(i, " done.\n")
    }
    
}

png("simest.png",width=640*pngscal,height=480*pngscal)
par(mfrow=c(2,2))
for(SPECIES in specLevels){
    surveyIndex:::plot.SIlist(ests[[SPECIES]],base=models[[SPECIES]],main=SPECIES,lwd=2)
}
dev.off()


