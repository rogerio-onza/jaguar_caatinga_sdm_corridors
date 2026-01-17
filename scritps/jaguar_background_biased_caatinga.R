# Generate background points
# Author: Oliveira, 2025
# Date: 01/13
# Version: 1.1

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

# Import sampling bias raster - 
# bias <- terra::rast("rasters/bias/onca_bias_caatinga.tif")

# Import environmental variable rasters
variables_files <- list.files("rasters/variables/", pattern = "\\.tif$", full.names = TRUE)
variables <- terra::rast(variables_files)

caatinga_iucn_buffer <- terra::vect("shps/jaguar_iucn_range_buffer.shp")
caatinga_ecoregion   <- terra::vect("shps/caatinga_ecoregion.shp")

# Data preparation / Projection check 
# Ensure coordinate systems match the raster layers (variables)
# Check and project both vectors if necessary

if (!terra::same.crs(variables, caatinga_iucn_buffer)) {
  caatinga_iucn_buffer <- terra::project(caatinga_iucn_buffer, variables)
}

if (!terra::same.crs(variables, caatinga_ecoregion)) {
  caatinga_ecoregion <- terra::project(caatinga_ecoregion, variables)
}

# 3. Analysis / Modeling (Background Sampling)

# Generate background for the IUCN BUFFER (Blue in plot)
data_bg_iucn <- flexsdm::sample_background(
  data = data,
  x = "Longitude",
  y = "Latitude",
  n = 5000,
  method = "random",
  rlayer = variables[[1]], 
  calibarea = caatinga_iucn_buffer # Uses the buffer shapefile
)

# Add column to identify origin
data_bg_iucn$bg_type <- "IUCN_Buffer"

# Generate background for the entire CAATINGA ECOREGION (Red in plot)
data_bg_caatinga <- flexsdm::sample_background(
  data = data,
  x = "Longitude",
  y = "Latitude",
  n = 5000,
  method = "random",
  rlayer = variables[[1]], 
  calibarea = caatinga_ecoregion # Uses the ecoregion shapefile
)
# Add column to identify origin
data_bg_caatinga$bg_type <- "Caatinga_Ecoregion"

# 4. Visualization

# Convert SpatVector to sf for ggplot2 compatibility
caatinga_sf <- sf::st_as_sf(caatinga_ecoregion)

# Plotting
ggplot() +
  # 1. Draw the Caatinga boundary (Background)
  geom_sf(data = caatinga_sf, fill = NA, color = "black", size = 0.8) +
  
  # 2. Caatinga Points (RED)
  geom_point(data = data_bg_caatinga, aes(x = Longitude, y = Latitude), 
             color = "red", alpha = 0.3, size = 0.5) +
  
  # 3. Buffer Points (BLUE) - Plotted last to appear on top
  geom_point(data = data_bg_iucn, aes(x = Longitude, y = Latitude), 
             color = "blue", alpha = 0.3, size = 0.5) +
  
  labs(title = "Background Points Distribution - Jaguar",
       subtitle = "Red: Caatinga Ecoregion | Blue: IUCN Buffer",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

# --- 5. Merge and Export results ---

# Merge the two dataframes
data_bg_merged <- dplyr::bind_rows(data_bg_caatinga, data_bg_iucn)

# Select ONLY Longitude and Latitude columns
data_bg_final <- data_bg_merged %>% 
  dplyr::select(Longitude, Latitude)

# Export the cleaned file
readr::write_csv(data_bg_final, "data/background/jaguar_data_bg_merged_max.csv")

