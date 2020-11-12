#' Create a systmatic grid for survey sampling
#' 
#' Takes the bounding box around the EBS survey shapefile and creates a regular square grid with dimensions sqSize X sqSize. Experimentation
#' led to adding lat and lon buffers to evenly spread edge samples from the prediction grid and to evenly divide full grid cells spatially. A
#' grid size of 32000 X 32000 m results in 16 available systematic samples per full grid cell.

require(tidyverse)
require(sf)

load(here::here("data","EBSbundle_1_2.rdata"))

# Ensure there is a shapefiles subdirectory
dir.create(here::here("data","shapefiles"))

# Import EBS polygon with no strata
EBStotal <- EBSbundle$EBSstrata
aeaProj <- CRS(proj4string(EBStotal))

sqSize <- 32000

sysGrid <- function(buffer_lat = 5000, buffer_lon = 15000) { 
  sqGrid <- st_make_grid(x = EBStotal,
                         cellsize = sqSize,
                         offset = st_bbox(EBStotal)[c("xmin", "ymin")] - c(buffer_lon, buffer_lat),
                         what = "polygons"
  ) %>%
    st_sf() %>%
    mutate(grid_id = row_number())
  
 
  
  return(sqGrid)
  
}


grid <- sysGrid() 


test = st_join(grid, sim$"2009") %>%
  dplyr::filter(!is.na(predict_id))

unique(test$grid_id)

#######################################################
ggplot() + 
 # geom_sf(data = EBStotal, fill = "blue") +
  geom_sf(data = sim$"2009", aes(color=simCPUE)) +
  geom_sf(data = grid, color="red", fill = NA)

st_write(grid,
         dsn = here::here("data","shapefiles"),
         layer = paste0("SystematicGrid_",sqSize),
         driver = "ESRI Shapefile",
         append = F,
         overwrite = T
)

st_write(sim$"2009",
         dsn = here::here("data","shapefiles"),
         layer = paste0("cod_sim_2009"),
         driver = "ESRI Shapefile"
)
