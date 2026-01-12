# Clip environmental variables for Caatinga
# Author: Oliveira, 2025
# Date: 01/12
# Version: 1.1

# 1. Load packages
library(terra)
library(dplyr)

# 2. Load data
# Load ecoregion vector (Update path if necessary)
caatinga_shp <- terra::vect("shps/caatinga_ecoregion.shp")

# Load environmental rasters (lists all .tif files in 'rasters/' folder)
predictors <- list.files("rasters/landcover_raw/", pattern = "\\.tif$", full.names = TRUE)
variables  <- terra::rast(predictors)

# 3. Spatial Processing

# Ensure Coordinate Reference Systems match
if (terra::crs(caatinga_shp) != terra::crs(variables)) {
  message("Reprojecting shapefile to match raster CRS...")
  caatinga_shp <- terra::project(caatinga_shp, terra::crs(variables))
}

# Crop and Mask variables
# Crop reduces extent; Mask sets values outside polygon to NA
message("Cropping and masking variables for Caatinga...")
vars_caatinga <- terra::crop(variables, caatinga_shp)
vars_caatinga <- terra::mask(vars_caatinga, caatinga_shp)

# 4. Export Results

# Create specific output directory
out_dir <- "rasters/bioclimatic"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  message(paste("Directory created:", out_dir))
}

# Save layers individually
message("Saving individual layers...")

for(i in 1:terra::nlyr(vars_caatinga)) {
  layer_name <- names(vars_caatinga)[i]
  file_path  <- file.path(out_dir, paste0(layer_name, ".tif"))
  
  terra::writeRaster(vars_caatinga[[i]], 
                     filename = file_path,
                     overwrite = TRUE)
}

message("Processing complete! Rasters saved in: rasters/caatinga/")