# Filter jaguar occurrences for Caatinga (No Telemetry)
# Author: Oliveira, 2025
# Date: 01/09
# Version: 1.1

# 1. Load packages
library(terra)      # Spatial data handling
library(readr)      # Data import/export
library(dplyr)      # Data manipulation
library(ggplot2)    # Visualization
library(tidyterra)  # geom_spatvector function

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
# Filter: Year >= 2010 AND Exact coordinates only
salve_clean <- salve_raw %>%
  dplyr::filter(Ano >= 2010, 
                Precisao_da_coordenada == "Exata")

# Fix Datum (SIRGAS2000 -> WGS84)
salve_sirgas <- salve_clean %>% dplyr::filter(Datum == "SIRGAS2000")
salve_other  <- salve_clean %>% dplyr::filter(Datum != "SIRGAS2000")

if (nrow(salve_sirgas) > 0) {
  # Convert to spatial vector using SIRGAS2000 (EPSG:4674)
  v_sirgas <- terra::vect(salve_sirgas, 
                          geom = c("Longitude", "Latitude"), 
                          crs = "EPSG:4674")
  
  # Project to WGS84 (EPSG:4326)
  v_wgs84  <- terra::project(v_sirgas, "EPSG:4326")
  
  # Update coordinates
  coords   <- terra::crds(v_wgs84)
  salve_sirgas$Longitude <- coords[, 1]
  salve_sirgas$Latitude  <- coords[, 2]
}
salve_final <- dplyr::bind_rows(salve_sirgas, salve_other)

# Merge datasets with source ID for tracking
df_list <- list(
  data_final %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "Dataset Final"),
  gbif_final %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "GBIF"),
  morato     %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "Morato 2014"),
  salve_final %>% dplyr::select(Longitude, Latitude) %>% dplyr::mutate(dataset = "Salve (Filtered)")
)

all_records <- dplyr::bind_rows(df_list)

# 4. Analysis (Spatial Filter)

# Convert all records to spatial vector
records_vect <- terra::vect(all_records, 
                            geom = c("Longitude", "Latitude"), 
                            crs = "EPSG:4326")

# Clip to Caatinga (intersect keeps only points inside polygon)
caatinga_points <- terra::intersect(records_vect, caatinga)

# Convert back to DataFrame
final_df_full <- as.data.frame(caatinga_points, geom = "XY") %>%
  dplyr::rename(Longitude = x, Latitude = y)

# Report counts per dataset inside Caatinga
cat("Records per dataset inside Caatinga:\n")
final_counts <- final_df_full %>% 
  dplyr::count(dataset) %>% 
  dplyr::arrange(desc(n))

print(final_counts)

# 5. Visualization
plot_map <- ggplot() +
  geom_spatvector(data = caatinga, fill = "gray95", color = "gray50") +
  geom_point(data = final_df_full, aes(x = Longitude, y = Latitude, color = dataset), 
             alpha = 0.7, size = 1.5) +
  scale_color_viridis_d() +
  theme_minimal() +
  labs(title = "Panthera onca in Caatinga", 
       subtitle = paste("Total N:", nrow(final_df_full)))

print(plot_map)

# 6. Export results
export_df <- final_df_full %>%
  dplyr::mutate(
    countryCode = "BR",
    Species = "Panthera onca"
  ) %>%
  dplyr::select(Longitude, Latitude, countryCode, Species)

readr::write_csv(export_df, "jaguar_data.csv")
cat("File saved: jaguar_data.csv\n")
