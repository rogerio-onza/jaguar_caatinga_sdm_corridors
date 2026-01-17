# Generate background points
# Author: Oliveira, 2025
# Date: 01/13
# Version: 1.0

# 1. Load packages
library(terra)
library(megaSDM)
library(ggplot2)
library(sf)
library(dplyr)
library(readr)

# Load Vectors
caatinga_iucn_buffer <- terra::vect("shps/jaguar_iucn_range_buffer.shp")
caatinga_ecoregion   <- terra::vect("shps/caatinga_ecoregion.shp")

# Import environmental variable rasters
variables_files <- list.files("rasters/variables/", pattern = "\\.tif$", full.names = TRUE)
variables <- terra::rast(variables_files)

# Prepare Environmental Data (The "Training Area")
# We crop and mask the variables to the Caatinga Ecoregion so megaSDM 
# understands this is the full available area.
env_caatinga <- terra::crop(variables, caatinga_ecoregion)
env_caatinga <- terra::mask(env_caatinga, caatinga_ecoregion)

# --- 2. Generate Background Points with megaSDM ---

# This function saves the file directly to the output folder
megaSDM::BackgroundPoints(
  spplist = c("jaguar"),      # Species name for filename
  envdata = env_caatinga,     # Raster masked to Caatinga
  output = "data/background/", # Folder to save CSV
  nbg = 105,                # Total number of points
  
  # THE KEY PARAMETER:
  # 0.6 means 60% of points come from the BUFFER.
  # The remaining 60% come from the whole envdata (Caatinga).
  spatial_weights = 0.6,       
  
  buffers = list(caatinga_iucn_buffer), # The shapefile for the 60%
  method = "random",              # "random" mimics your previous logic. Use "Varela" for environmental filtering.
  ncores = 1
)

output_folder <- "data/background/"
nome_desejado <- "jaguar_bg_nnet_gbm_rf.csv"

# 2. Renomear o arquivo gerado para o nome que você quer
arquivo_original <- file.path(output_folder, "jaguar_background.csv") # Nome padrão do pacote
arquivo_novo     <- file.path(output_folder, nome_desejado)   # Seu nome personalizado

if (file.exists(arquivo_original)) {
  file.rename(from = arquivo_original, to = arquivo_novo)
  message(paste("Arquivo renomeado com sucesso para:", nome_desejado))
} else {
  warning("O arquivo original não foi encontrado. Verifique se o nome da espécie em 'spplist' é 'Jaguar'.")
}
