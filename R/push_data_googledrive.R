#' DO NOT USE - Uploads local files to Google Drive
#' 
#' As of 11/16/2020, the overwrite argument does not work for upload. This is an open
#' issue for *googledrive*: https://github.com/tidyverse/googledrive/issues/312
#' 
#' Instead, use the internet interface to upload and replace files.

library(googledrive)
library(tidyverse)

# this should fire up browser window
drive_auth()


local_data <- list.files(here::here("data"), full.names = TRUE) 
local_shapefiles <- list.files(here::here("shapefiles"), full.names = TRUE) 

g_data <- drive_get("MoBaSim/data")
g_shapefiles <- drive_get("MoBaSim/shapefiles")

upload_data <- purrr::map(local_data, ~ drive_upload(.x, path = as_id(g_data), verbose = TRUE, overwrite = TRUE))
upload_shapefiles <- purrr::map(local_shapefiles, ~ drive_upload(.x, path = as_id(g_shapefiles), verbose = TRUE, overwrite = TRUE))
  
  
  