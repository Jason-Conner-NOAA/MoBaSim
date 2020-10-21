library(googledrive)

dir.create(here::here("shapefiles"))
dir.create(here::here("data"))

# this should fire up browser window
drive_auth()

# specify url of folder you want to pull in
shapefiles <- drive_ls(as_id("https://drive.google.com/drive/folders/11nnb9rasPmTMRYyrUKlewLw9ek_YEvbL"))
data <- drive_ls(as_id("https://drive.google.com/drive/folders/1ZyjInG4AaJetqsdVw4zSm9GmyNwDk4sh"))


# loop over files in folder if need to download many (or can do by file name, etc)
for(i in 1:nrow(shapefiles)) {
  drive_download(file=shapefiles[i,], 
    path=here::here("shapefiles",shapefiles$name[i]), 
    overwrite = TRUE)
}

for(i in 1:nrow(data)) {
   drive_download(file=data[i,], 
                  path=here::here("data",data$name[i]), 
                  overwrite = TRUE)
}