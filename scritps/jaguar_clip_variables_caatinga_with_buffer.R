# Clip selected environmental variables (post-VIF) for Caatinga with buffer
# Author: Oliveira, 2025
# Date: 01/22
# Version: 1.0 - Buffer version with VIF-selected variables only
# Description: Clips only the variables that passed VIF analysis

# 1. Load packages ----
library(terra)
library(dplyr)

# 2. Define paths and variables ----

# Input directories
landcover_dir <- "rasters/landcover_raw/"
bioclim_dir   <- "rasters/bioclimatic_raw/"

# Output directory
out_dir <- "rasters/variables_buffer"

# Shapefile with buffer (Update path to your buffer shapefile)
buffer_shp_path <- "shps/caatinga_ibge/caatinga_ibge_buffer_50km.shp"

# Variables that passed VIF analysis
bioclim_vars <- c("bio2", "bio3", "bio8", "bio12", "bio14", "bio18", "bio19")
landcover_vars <- c("class5", "class6", "class11", "class12")

# 3. Create output directory ----
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# 4. Load buffer shapefile ----
caatinga_buffer <- terra::vect(buffer_shp_path)

# 5. Process BIOCLIMATIC variables ----
message("Processing BIOCLIMATIC variables...")

for (var in bioclim_vars) {
  # Find file matching the variable name (with or without leading zero)
  var_num <- gsub("bio", "", var)
  var_padded <- sprintf("bio%02d", as.numeric(var_num))  # Add leading zero if needed
  file_pattern <- paste0("CHELSA_", var_padded, ".*\\.tif$")
  var_file <- list.files(bioclim_dir, pattern = file_pattern, full.names = TRUE)
  
  if (length(var_file) == 0) {
    warning(paste("File not found for:", var))
    next
  }
  
  if (length(var_file) > 1) {
    var_file <- var_file[1]
  }
  
  # Load raster
  raster_var <- terra::rast(var_file)
  
  # Reproject shapefile if necessary
  if (terra::crs(caatinga_buffer) != terra::crs(raster_var)) {
    buffer_proj <- terra::project(caatinga_buffer, terra::crs(raster_var))
  } else {
    buffer_proj <- caatinga_buffer
  }
  
  # Crop and Mask
  var_clipped <- terra::crop(raster_var, buffer_proj)
  var_clipped <- terra::mask(var_clipped, buffer_proj)
  
  # Save with simplified name
  out_file <- file.path(out_dir, paste0(var, ".tif"))
  terra::writeRaster(var_clipped, 
                     filename = out_file,
                     overwrite = TRUE)
  message(paste("✓", var))
}


# 6. Process LAND COVER variables ----
message("Processing LAND COVER variables...")

for (var in landcover_vars) {
  # Find file matching the variable name (with underscore before number)
  class_num <- gsub("class", "", var)
  file_pattern <- paste0("consensus_full_class_", class_num, "\\.tif$")
  var_file <- list.files(landcover_dir, pattern = file_pattern, full.names = TRUE)
  
  if (length(var_file) == 0) {
    warning(paste("File not found for:", var))
    next
  }
  
  if (length(var_file) > 1) {
    var_file <- var_file[1]
  }
  
  # Load raster
  raster_var <- terra::rast(var_file)
  
  # Reproject shapefile if necessary
  if (terra::crs(caatinga_buffer) != terra::crs(raster_var)) {
    buffer_proj <- terra::project(caatinga_buffer, terra::crs(raster_var))
  } else {
    buffer_proj <- caatinga_buffer
  }
  
  # Crop and Mask
  var_clipped <- terra::crop(raster_var, buffer_proj)
  var_clipped <- terra::mask(var_clipped, buffer_proj)
  
  # Save with simplified name
  out_file <- file.path(out_dir, paste0(var, ".tif"))
  terra::writeRaster(var_clipped, 
                     filename = out_file,
                     overwrite = TRUE)
  message(paste("✓", var))
}

# 7. Summary ----
message("\nProcessing complete!")
message(paste("Output:", out_dir))