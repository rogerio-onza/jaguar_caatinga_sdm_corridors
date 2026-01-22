# Clip and align environmental variables for Caatinga with buffer
# Author: Oliveira, 2025
# Date: 01/22
# Version: 3.0 - Integrated clip + alignment
# Description: Clips VIF-selected variables with buffer and aligns all grids

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

# ============================================================================
# PART 1: PROCESS BIOCLIMATIC VARIABLES
# ============================================================================

message("Processing BIOCLIMATIC variables...")

bio_files_created <- c()

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
  
  # Rename layer
  names(var_clipped) <- var
  
  # Save with simplified name
  out_file <- file.path(out_dir, paste0(var, ".tif"))
  terra::writeRaster(var_clipped, 
                     filename = out_file,
                     overwrite = TRUE)
  bio_files_created <- c(bio_files_created, out_file)
  message(paste("✓", var))
}

# ============================================================================
# PART 2: PROCESS LAND COVER VARIABLES (initial clip)
# ============================================================================

message("Processing LAND COVER variables...")

class_temp_files <- c()

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
  
  # Rename layer
  names(var_clipped) <- var
  
  # Save temporarily
  temp_file <- file.path(out_dir, paste0(var, "_temp.tif"))
  terra::writeRaster(var_clipped, 
                     filename = temp_file,
                     overwrite = TRUE)
  class_temp_files <- c(class_temp_files, temp_file)
  message(paste("✓", var, "(temporary)"))
}

# ============================================================================
# PART 3: ALIGN LAND COVER TO BIOCLIMATIC GRID
# ============================================================================

message("Aligning land cover to bioclimatic grid...")

# Load bioclimatic reference
bio_reference <- terra::rast(bio_files_created[1])

# Load temporary land cover files
class_stack <- terra::rast(class_temp_files)

# Resample to bioclimatic grid
class_aligned <- terra::resample(class_stack, bio_reference, method = "near")

# Save aligned files with final names
for(i in 1:terra::nlyr(class_aligned)) {
  layer_name <- names(class_aligned)[i]
  out_file <- file.path(out_dir, paste0(layer_name, ".tif"))
  
  terra::writeRaster(class_aligned[[i]], 
                     filename = out_file,
                     overwrite = TRUE)
  message(paste("✓", layer_name, "(aligned)"))
}

# Remove temporary files
file.remove(class_temp_files)

# ============================================================================
# PART 4: FINAL VERIFICATION
# ============================================================================

message("\nVerifying final stack...")

variables_files <- list.files(out_dir, pattern = "\\.tif$", full.names = TRUE)
variables <- terra::rast(variables_files)
