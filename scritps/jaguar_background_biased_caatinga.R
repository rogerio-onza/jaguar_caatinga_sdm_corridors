# Generate background points
# Author: Oliveira, 2025
# Date: 01/13
# Version: 1.2

# 1. Load packages
library(terra)
library(dplyr)
library(flexsdm)
library(readr)
library(ggplot2)
library(sf)

# 2. Load data
# Import raw occurrence data
data <- readr::read_csv("data/presence/jaguar_data.csv")

#Import sampling bias raster - 
bias <- terra::rast("rasters/bias/onca_bias_caatinga.tif")

# Import environmental variable rasters
variables_files <- list.files("rasters/variables/", pattern = "\\.tif$", full.names = TRUE)
variables <- terra::rast(variables_files)

caatinga_ecoregion   <- terra::vect("shps/caatinga_ecoregion.shp")

# Data preparation / Projection check 
# Ensure coordinate systems match the raster layers (variables)
# Check and project both vectors if necessary
if (!terra::same.crs(variables, caatinga_ecoregion)) {
  caatinga_ecoregion <- terra::project(caatinga_ecoregion, variables)
}

# 3. Analysis / Modeling (Background Sampling)

# Generate background for the IUCN BUFFER (Blue in plot)
data_bg_ <- flexsdm::sample_background(
  data = data,
  x = "Longitude",
  y = "Latitude",
  n = 10000,
  method = "biased",
  rlayer = variables[[1]],
  rbias = bias,
  calibarea = caatinga_ecoregion 
)

# Export the cleaned file
readr::write_csv(data_bg, "data/background/jaguar_data_bg.csv")

