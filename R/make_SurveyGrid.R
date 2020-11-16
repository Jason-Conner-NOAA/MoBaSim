#' Create a systmatic grid for survey sampling
#' 
#' Takes the bounding box around the EBS survey shapefile and creates a regular square grid with dimensions sqSize X sqSize. Experimentation
#' led to adding lat and lon buffers to evenly spread edge samples from the prediction grid and to evenly divide full grid cells spatially. A
#' grid size of 32000 X 32000 m results in 16 available systematic samples per full grid cell. Below this code is code for assigning the index
#' for a systematic sample within a single grid, the output file is an index of grid_id, subgrid_id, and predict_id, for use in simulating
#' systematic surveys.

require(tidyverse)
require(sf)

load(here::here("data","EBSbundle_1_2.rdata"))

# Ensure there is a shapefiles subdirectory
dir.create(here::here("shapefiles"))

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
         dsn = here::here("shapefiles"),
         layer = paste0("SystematicGrid_",sqSize),
         driver = "ESRI Shapefile",
         append = F,
         overwrite = T
)

st_write(sim$"2009",
         dsn = here::here("shapefiles"),
         layer = paste0("cod_sim_2009"),
         driver = "ESRI Shapefile"
)


# Create Index for systematic points within 1 grid ------------------------
predictGrid <- EBSbundle$EBSfullPredict %>%
  mutate(LAT_m= LAT, LON_m=LONG) %>%
  st_as_sf(coords=c("LONG","LAT"), 
           crs="+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs") 

grid_mod <- st_join(predictGrid,grid ) %>%
  mutate(lon_m = st_coordinates(.)[,1],
         lat_m = st_coordinates(.)[,2]
  ) %>%
  st_drop_geometry()

grids <- unique(grid_mod$grid_id)

grid_index <- data.frame()

for (g in grids) {
         subGrid <- dplyr::filter(grid_mod, grid_id == g) 
         
         if (min(subGrid$lon_m) > -1340600 & max(subGrid$lon_m) < -260000) {
           full_grid <- data.frame(subgrid_id = seq(1:16))
           index <- subGrid %>%
             dplyr::arrange(desc(lat_m),lon_m) %>% 
             dplyr::bind_cols(full_grid) 
         } else if (min(subGrid$lon_m) < -1340600) {
           # western edge subgrid ids
           w_grid <- data.frame(subgrid_id = c(3,4,7,8,11,12,15,16))
           index <- subGrid %>%
             dplyr::arrange(desc(lat_m),lon_m) %>% 
             dplyr::bind_cols(w_grid) 
         } else if (max(subGrid$lon_m) > -260000) {
           # eastern edge subgrid ids
           e_grid <- data.frame(subgrid_id = c(1,2,3,5,6,7,9,10,11,13,14,15))
           index <- subGrid %>%
             dplyr::arrange(desc(lat_m),lon_m) %>% 
             dplyr::bind_cols(e_grid) 
         }
         
         index_all <- dplyr::select(index, grid_id, subgrid_id, predict_id) 
         grid_index <- bind_rows(grid_index, index_all)
        
       }

saveRDS(grid_index, here::here("data","grid_index.rds"))

# ggplot() + 
#   geom_sf(data = grid_mod)



# write.csv(EBSbundle$EBSfullPredict,
#          here::here("shapefiles","PredictGrid.csv")
# )