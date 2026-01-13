# Generate background points
# Author: Oliveira, 2025
# Date: 01/13
# Version: 1.0

# 1. Load packages
library(terra)
library(dplyr)
library(flexsdm)
library(readr)
library(ggplot2)
library(sf)

# 2. Load data
# Import raw occurrence data
data <- readr::read_csv("data/jaguar_data.csv")

# Import sampling bias raster
bias <- terra::rast("rasters/bias/onca_bias_caatinga.tif")

# Import environmental variable rasters
variables_files <- list.files("rasters/variables/", pattern = "\\.tif$", full.names = TRUE)
variables <- terra::rast(variables_files)

# Import calibration area vector (Shapefile)
caatinga_shp <- terra::vect("shps/caatinga_ecoregion.shp")

# 3. Data preparation/cleaning
# Ensure coordinate systems match before sampling
if (!terra::same.crs(variables, caatinga_shp)) {
  caatinga_shp <- terra::project(caatinga_shp, variables)
}

# 4. Analysis/modeling
# Generate background for GBIF data
# Note: Using the first raster layer as reference for resolution/extent
data_bg <- flexsdm::sample_background(
  data = data,
  x = "Longitude",
  y = "Latitude",
  n = 10000,
  method = "biased",
  rlayer = variables[[1]], 
  rbias = bias,
  calibarea = caatinga_shp
)

# 5. Visualization
# Convert SpatVector to sf for ggplot2 compatibility
caatinga_sf <- sf::st_as_sf(caatinga_shp)

# Visualizing generated background points
ggplot() +
  geom_sf(data = caatinga_sf, fill = NA, color = "black") +
  geom_point(data = data_bg, aes(x = Longitude, y = Latitude), color = "red", alpha = 0.3, size = 0.5) +
  labs(title = "Background Points Distribution - GBIF (Caatinga)",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# 6. Export results
readr::write_csv(data_bg, "data/background/jaguar_data_bg.csv")
