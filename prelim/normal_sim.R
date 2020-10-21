library(sumfish)
library(FNN)
getSQL()

rb <- getRacebase(2019,'EBS_SHELF')
haul <- sumHaul(rb)

atf <- filter(haul, SPECIES_CODE==10110)

nn <- get.knn(select(atf,START_LONGITUDE,START_LATITUDE),k=4)

atfSim <- atf
atfSim$wCPUESD <- NA

for (i in 1:nrow(atf)) {
  nn.SD <- sqrt(var(atf[nn$nn.index[i,],"wCPUE"]))
  if (nn.SD > 0) {
    ifelse(atfSim[i,"wCPUE"] != 0, atfSim[i,"wCPUESD"] <- nn.SD, NA)
    
  } else {
    message("Error - zero variance")
  }
}

simAtf <- atfSim
simAtf$simwCPUE <- NA

for (x in 1:nrow(atf)) {
  if (is.na(atfSim[x,"wCPUESD"])) {
    simAtf[x,"simwCPUE"] <- 0
  } else {
    simAtf[x,"simwCPUE"] <- rnorm(1,atfSim[x,"wCPUE"],atfSim[x,"wCPUESD"])
    
    
plot( density(rlnorm(1,(atfSim[x,"wCPUE"]),.1)  )) 