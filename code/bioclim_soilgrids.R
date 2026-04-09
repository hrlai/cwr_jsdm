# Retrieve data from WorldClim and SoilGrids
library(terra)
library(geodata)
library(raster)
library(gdalUtils)



# bioclim -----------------------------------------------------------------

bioclim <- worldclim_global(
    var = "bio", 
    res = 10, 
    path = "data/WorldClim")





# soilgrids ---------------------------------------------------------------

# paths to the webDAV data
url <- "/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/" 
voi <- c("bdod", "cec", "clay", "soc", "phh2o", "silt", "sand")
depth <- "0-5cm"
quantile <- "mean"

variable <- paste(url, voi, sep="")
layer <- paste(voi, depth, quantile, sep="_")
vrt_layer <- paste0(variable, "/", layer, '.vrt')
filename <- paste0("data/soilgrids/", layer, ".tif")

# download, project, resample and then store locally as geoTIFF
overwrite <- FALSE
for (i in seq_along(vrt_layer)) {
  # solution from 
  # https://gis.stackexchange.com/a/361350
  # will error if TIFF exists, which is what we want
  gdalwarp(
    vrt_layer[i],  # Input VRT
    filename[i],    # Output TIFF        
    t_srs = crs(bioclim, proj = TRUE),
    multi = TRUE,
    r = "average",
    wm = 200,
    co = c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"),
    te = as.vector(ext(bioclim))[c(1,3,2,4)],
    tr = res(bioclim),  # output resolution (same as bioclim)
    verbose = T,
    overwrite = overwrite
  ) 
}

# read the rasters
soilgrids_files <- list.files("data/soilgrids", full.names = TRUE)
soilgrids_rast <- lapply(soilgrids_files, rast)
names(soilgrids_rast) <- lapply(soilgrids_rast, names)
soilgrids_rast <- rast(soilgrids_rast)

# plot(soilgrids_rast)

bioclim_soilgrids <- c(bioclim, soilgrids_rast)
