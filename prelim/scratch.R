  tGrid <- st_make_grid(st_buffer(atf, 1000),
                        cellsize = 18000,
                        what = "centers") %>%
    st_sf()
  
  plot(tGrid, axes = TRUE)
  plot(atf, add = TRUE, col = "red")

  tatf <- filter(atf, YEAR==1999)  
  
  tIDW <- idw(wCPUE~1,tatf,newdata=EBSpredict, idp=1)
              
plot(st_geometry(tIDW), col = tIDW$var1.pred)
sum(is.na(tIDW$var1.pred)) 
plot(st_geometry(tatf))

tIDW <- idwSim(tatf,"wCPUE",EBSpredict)

# sf seems to not work as expected, try sp
library(sp)

tatf <- atf
sp::coordinates(tatf) <- ~START_LONGITUDE + START_LATITUDE
sp::proj4string(tatf) <- proj4string(EBSpredict)

spIDW <- idw(wCPUE~1,tatf,newdata=EBSpredict, idp=2)
rIDW <- raster(spIDW)
