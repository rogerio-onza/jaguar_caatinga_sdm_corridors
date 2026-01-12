# Data and variable selection (vif < 3)
# Region: Caatinga
# Author: Oliveira, 2025
# Version: 1.0

# Load packages
library(terra)
library(dplyr)
library(flexsdm)

# Define reference grid
bio_files <- list.files("rasters/bioclimatic", pattern = "\\.tif$", full.names = TRUE)
if(length(bio_files) == 0) stop("Bioclimatic files not found in 'rasters/bioclimatic'.")

ref_raster <- terra::rast(bio_files[1])

# Process elevation (90m -> 1km, WGS84)
message("Processing elevation...")
elev_files <- list.files("rasters/elevation", pattern = "\\.tif$", full.names = TRUE)
if(length(elev_files) == 0) stop("Elevation file not found in 'rasters/elevation'.")

raw_elev <- terra::rast(elev_files[1])
elev_1km <- terra::project(raw_elev, ref_raster, method = "bilinear")
names(elev_1km) <- "elevation"

if(!dir.exists("rasters/processed_1km")) dir.create("rasters/processed_1km", recursive = TRUE)
terra::writeRaster(elev_1km, "rasters/processed_1km/elevation_1km.tif", overwrite = TRUE)

# Load and rename variables
bio_stack <- terra::rast(bio_files)
bio_names_new <- gsub("CHELSA_bio0?(\\d+)_.*", "bio\\1", names(bio_stack))
names(bio_stack) <- bio_names_new

land_files <- list.files("rasters/landcover", pattern = "\\.tif$", full.names = TRUE)
if(length(land_files) == 0) stop("Landcover files not found in 'rasters/landcover'.")

land_stack <- terra::rast(land_files)
land_names_new <- gsub("consensus_full_class_(\\d+)", "class\\1", names(land_stack))
names(land_stack) <- land_names_new

# Filter landcover: keep only ecologically relevant classes for Caatinga
# Remove: class2, class3 (forest types not present in Caatinga)
classes_to_keep <- c("class1", "class4", "class5", "class6", "class7", "class9", "class11", "class12")
available_classes <- names(land_stack)
classes_keep <- intersect(classes_to_keep, available_classes)

if(length(classes_keep) > 0) {
  land_stack <- land_stack[[classes_keep]]
  message(paste("Landcover filtered:", length(classes_keep), "classes kept"))
} else {
  stop("No relevant landcover classes found.")
}

env_stack_full <- c(bio_stack, land_stack, elev_1km)

# Define protected variables (Morato et al. 2014)
protected_vars <- c("bio14")
candidate_vars <- setdiff(names(env_stack_full), protected_vars)

# Step 1: vif on non-protected variables
message("\nRunning vif < 3 on non-protected variables...")
env_candidates <- env_stack_full[[candidate_vars]]

vif_step1 <- flexsdm::correct_colinvar(
  env_layer = env_candidates,
  method = c("vif", th = "3"),
  maxcell = 50000
)

message("\nStep 1 - Variables removed:")
print(vif_step1$removed_variables)
print(vif_step1$vif_table)

# Step 2: combine selected + protected variables
protected_stack <- env_stack_full[[protected_vars]]
final_stack <- c(protected_stack, vif_step1$env_layer)

message(paste("\nFinal stack:", nlyr(final_stack), "variables"))

# Final vif check (reporting only)
vif_final <- flexsdm::correct_colinvar(
  env_layer = final_stack,
  method = c("vif", th = "3"),
  maxcell = 50000
)

message("\nFinal vif table:")
print(vif_final$vif_table)

# Export final variables (only those with vif < 3)
message("\nExporting variables to 'rasters/variables'...")
out_dir <- "rasters/variables"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

final_selected <- vif_final$env_layer

for(i in 1:terra::nlyr(final_selected)) {
  var_name <- names(final_selected)[i]
  output_path <- file.path(out_dir, paste0(var_name, ".tif"))
  terra::writeRaster(final_selected[[i]], output_path, overwrite = TRUE)
}

message(paste("\nCompleted! Final variables:", paste(names(final_selected), collapse = ", ")))
