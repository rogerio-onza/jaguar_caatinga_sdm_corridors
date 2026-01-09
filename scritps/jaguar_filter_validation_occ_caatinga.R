# Filter and Clean jaguar occurrences for Caatinga
# Author: Oliveira, 2025
# Date: 01/09
# Version: 1.2 - Added CoordinateCleaner validation and deduplication

# 1. Load packages
library(terra)              # Spatial data
library(readr)              # Import/Export
library(dplyr)              # Data manipulation
library(ggplot2)            # Visualization
library(tidyterra)          # Spatial plotting
library(CoordinateCleaner)  # Coordinate validation
library(countrycode)        # ISO code conversion

# 2. Load data
# Shapefile
caatinga <- terra::vect("shps/caatinga_ecoregion.shp")

# CSV Datasets
data_final <- readr::read_csv("data/dataset_final.csv")
gbif_final <- readr::read_csv("data/gbif_final.csv")
morato     <- readr::read_csv("data/jaguar_wgs84_morato_2014.csv")
salve_raw  <- readr::read_csv("data/jaguar_salve_raw.csv")

# 3. Data preparation/cleaning

# Ensure shapefile is WGS84
caatinga <- terra::project(caatinga, "EPSG:4326")

# Process: Jaguar Salve
salve_clean <- salve_raw %>%
  dplyr::filter(Ano >= 2010, Precisao_da_coordenada == "Exata")

# Fix Datum (SIRGAS2000 -> WGS84)
salve_sirgas <- salve_clean %>% dplyr::filter(Datum == "SIRGAS2000")
salve_other  <- salve_clean %>% dplyr::filter(Datum != "SIRGAS2000")

if (nrow(salve_sirgas) > 0) {
  v_sirgas <- terra::vect(salve_sirgas, geom = c("Longitude", "Latitude"), crs = "EPSG:4674")
  v_wgs84  <- terra::project(v_sirgas, "EPSG:4326")
  coords   <- terra::crds(v_wgs84)
  salve_sirgas$Longitude <- coords[, 1]
  salve_sirgas$Latitude  <- coords[, 2]
}
salve_final <- dplyr::bind_rows(salve_sirgas, salve_other)

# Merge datasets (Keep 'dataset' column for initial report)
df_list <- list(
  data_final %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "Dataset Final"),
  gbif_final %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "GBIF"),
  morato     %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "Morato 2014"),
  salve_final %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "Salve (Filtered)")
)

all_records <- dplyr::bind_rows(df_list)

# 4. Analysis (Spatial Filter)
records_vect <- terra::vect(all_records, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
caatinga_points <- terra::intersect(records_vect, caatinga)

# Convert back to DataFrame
spatial_df <- as.data.frame(caatinga_points, geom = "XY") %>%
  dplyr::rename(Longitude = x, Latitude = y)

# 5. Coordinate Cleaner & Deduplication

# Prepare data for CoordinateCleaner
# Adding temporary columns for validation
clean_input <- spatial_df %>%
  dplyr::mutate(
    countryCode = "BR",           # Initial ISO2
    Species = "Panthera onca"
  )

# Convert ISO2 to ISO3 (Required for CoordinateCleaner)
clean_input$countryCode <- countrycode::countrycode(
  clean_input$countryCode, 
  origin = "iso2c", 
  destination = "iso3c"
)

# Run CoordinateCleaner (Standard validation)
# Removed 'seas' and 'outliers' as spatial clip already handles domain logic
# retained checking for biological impossibilities (equal, zeros) and administrative centroids
flags <- CoordinateCleaner::clean_coordinates(
  x = clean_input,
  lon = "Longitude",
  lat = "Latitude",
  countries = "countryCode",
  species = "Species",
  tests = c("capitals", "centroids", "equal", "zeros"),
  verbose = FALSE
)

# Filter flagged records
data_validated <- clean_input[flags$.summary, ]

# Remove exact duplicates (Spatial Deduplication)
data_unique <- data_validated %>%
  dplyr::distinct(Longitude, Latitude, .keep_all = TRUE)

# Report processing stats
cat("Processing Summary:\n")
cat("Total inside Caatinga:", nrow(spatial_df), "\n")
cat("Post-validation (CC):", nrow(data_validated), "\n")
cat("Final Unique Records:", nrow(data_unique), "\n")

# 6. Visualization
plot_map <- ggplot() +
  geom_spatvector(data = caatinga, fill = "gray95", color = "gray50") +
  geom_point(data = data_unique, aes(x = Longitude, y = Latitude), 
             color = "darkred", alpha = 0.6, size = 1.5) +
  theme_minimal() +
  labs(title = "Panthera onca - Caatinga Final", 
       subtitle = paste("N =", nrow(data_unique)))

print(plot_map)

# 7. Export results
# Format strictly as requested
data_export <- data_unique %>%
  dplyr::mutate(countryCode = "BR") %>% # Revert to BR if you prefer ISO2 for final csv, or keep BRA
  dplyr::select(Longitude, Latitude)

readr::write_csv(data_export, "jaguar_data.csv")
cat("File saved: jaguar_data.csv\n")
