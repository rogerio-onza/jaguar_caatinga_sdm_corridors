# Script: Full Benchmark - 50 Reps + Null Models + Tuning
# Region: Caatinga | Species: Panthera onca
# Author: Oliveira, 2025
# Version: 3.1 - Minor text fixes

# 1. Load packages
library(terra)
library(dplyr)
library(readr)
library(ggplot2)
library(flexsdm)
library(GeoThinneR)
library(future)

# 2. Settings
# Computational Setup
# Allocating 2 cores for parallel tuning inside flexsdm functions
future::plan("multisession", workers = 8)

# Parameters
species_name <- "Panthera onca"
input_csv    <- "jaguar_data.csv"
raster_dir   <- "rasters/variables"
output_dir   <- "results/full_benchmark_v3"

hr_value     <- 9.0
multipliers  <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
n_replicates <- 50
seed_base    <- 100

# Create Directories
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 3. PREPARE DATA ---------------------------------------------------------
message(">>> Loading Data...")

# Load Occurrences
occ_raw <- readr::read_csv(input_csv) %>%
  dplyr::select(x = Longitude, y = Latitude) %>%
  dplyr::mutate(pr_ab = 1) %>%
  dplyr::filter(!is.na(x) & !is.na(y))

# Load Rasters
env_layers <- terra::rast(list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE))
names(env_layers) <- gsub("[^[:alnum:]_]", "", names(env_layers))

# 4. DEFINE HYPERPARAMETER GRIDS ------------------------------------------
grid_rf <- expand.grid(mtry = 1:7, splitrule = "gini", min.node.size = 1)
grid_svm <- expand.grid(C = c(2, 4, 8, 16, 20), sigma = c(0.01, 0.1, 0.2, 0.3, 0.4))
grid_max <- expand.grid(regmult = c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4), classes = c("l", "lq", "lqp","lqh", "lqph", "lqhpt"))
grid_gbm <- expand.grid(n.trees = c(20, 50, 100), shrinkage = c(0.1, 0.5, 1), n.minobsinnode = c(1, 3, 5, 7, 9))

# 5. HELPER FUNCTION: FIT & EVALUATE --------------------------------------
fit_and_evaluate <- function(data_input, env_stack, algo_grids, scenario_label, rep_id, mult_val, dist_val) {
  
  metrics_list <- list()
  
  # Prepare Data (Pseudo-absences 1:1 & Block CV)
  # Using tryCatch to avoid stopping long loops on single errors
  d_sdm <- tryCatch({
    d <- flexsdm::sample_pseudoabs(data_input, x = "x", y = "y", n = nrow(data_input), method = "random", rlayer = env_stack, calibarea_method = NULL)
    d <- flexsdm::sdm_extract(d, x = "x", y = "y", env_layer = env_stack)
    flexsdm::part_sblock(d, x = "x", y = "y", pr_ab = "pr_ab", n_folds = 4)
  }, error = function(e) return(NULL))
  
  if(is.null(d_sdm)) return(NULL)
  
  # --- A. Random Forest ---
  m_rf <- try(flexsdm::fit_raf(d_sdm, "pr_ab", names(env_stack), ".part", grid_rf, thr = "max_sens_spec"), silent = TRUE)
  if (!inherits(m_rf, "try-error")) {
    metrics_list[[length(metrics_list)+1]] <- m_rf$performance %>% 
      mutate(Algorithm = "RF", Scenario = scenario_label, Replicate = rep_id, Multiplier = mult_val, Distance = dist_val)
  }
  
  # --- B. SVM ---
  m_svm <- try(flexsdm::fit_svm(d_sdm, "pr_ab", names(env_stack), ".part", grid_svm, thr = "max_sens_spec"), silent = TRUE)
  if (!inherits(m_svm, "try-error")) {
    metrics_list[[length(metrics_list)+1]] <- m_svm$performance %>% 
      mutate(Algorithm = "SVM", Scenario = scenario_label, Replicate = rep_id, Multiplier = mult_val, Distance = dist_val)
  }
  
  # --- C. MaxEnt ---
  m_max <- try(flexsdm::fit_max(d_sdm, "pr_ab", names(env_stack), ".part", grid_max, thr = "max_sens_spec", clamp = TRUE), silent = TRUE)
  if (!inherits(m_max, "try-error")) {
    metrics_list[[length(metrics_list)+1]] <- m_max$performance %>% 
      mutate(Algorithm = "MaxEnt", Scenario = scenario_label, Replicate = rep_id, Multiplier = mult_val, Distance = dist_val)
  }
  
  # --- D. GBM ---
  m_gbm <- try(flexsdm::fit_gbm(d_sdm, "pr_ab", names(env_stack), ".part", grid_gbm, thr = "max_sens_spec"), silent = TRUE)
  if (!inherits(m_gbm, "try-error")) {
    metrics_list[[length(metrics_list)+1]] <- m_gbm$performance %>% 
      mutate(Algorithm = "GBM", Scenario = scenario_label, Replicate = rep_id, Multiplier = mult_val, Distance = dist_val)
  }
  
  return(bind_rows(metrics_list))
}

# 6. MAIN LOOP ------------------------------------------------------------
all_results <- data.frame()
message(">>> Starting Benchmark Loop (50 Replicates per scenario)...")

for (m in multipliers) {
  
  dist_km <- m * hr_value
  message(paste0("\nProcessing Multiplier: ", m, "x (", dist_km, " km)"))
  
  for (rep in 1:n_replicates) {
    curr_seed <- seed_base + rep
    
    # 1. THINNING (REAL DATA)
    if (m == 0) {
      occ_thinned <- occ_raw
    } else {
      # Use GeoThinneR with seed
      # Note: GeoThinneR usually deterministic, but trials included
      thin_out <- GeoThinneR::thin_points(occ_raw, "x", "y", method = "distance", thin_dist = dist_km, trials = 1, seed = curr_seed, verbose = FALSE)
      occ_thinned <- occ_raw[thin_out$retained[[1]], ]
    }
    
    # Skip if too few points
    if (nrow(occ_thinned) < 10) next
    
    # 2. NULL MODEL DATA
    # Generate random points (same N as thinned) within raster extent
    # This simulates "what if the points were just random?"
    set.seed(curr_seed)
    occ_null <- terra::spatSample(env_layers[[1]], size = nrow(occ_thinned), method = "random", na.rm = TRUE, as.points = FALSE, xy = TRUE, values = FALSE)
    occ_null <- as.data.frame(occ_null) %>% mutate(pr_ab = 1)
    
    # 3. FIT MODELS (REAL)
    res_real <- fit_and_evaluate(occ_thinned, env_layers, NULL, "Real", rep, m, dist_km)
    
    # 4. FIT MODELS (NULL)
    res_null <- fit_and_evaluate(occ_null, env_layers, NULL, "Null", rep, m, dist_km)
    
    # Combine and Save
    batch_res <- bind_rows(res_real, res_null)
    all_results <- bind_rows(all_results, batch_res)
    
    # Partial Save (Safety)
    if (rep %% 5 == 0) {
      readr::write_csv(all_results, file.path(output_dir, "benchmark_partial_results.csv"))
      cat(".") # Progress indicator
    }
  }
}

# 7. EXPORT & VISUALIZE ---------------------------------------------------
message("\n>>> Finalizing and Plotting...")

# Save Final Table
readr::write_csv(all_results, file.path(output_dir, "benchmark_final_results.csv"))

# Comparison Plot (Null vs Real)
# We want to see if Real TSS/CBI is significantly higher than Null
plot_data <- all_results %>%
  dplyr::filter(Scenario %in% c("Real", "Null")) %>%
  dplyr::select(Algorithm, Multiplier, Scenario, TSS = tss, CBI = boyce) %>%
  tidyr::pivot_longer(cols = c(TSS, CBI), names_to = "Metric", values_to = "Value")

p_comp <- ggplot(plot_data, aes(x = factor(Multiplier), y = Value, fill = Scenario)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_grid(Metric ~ Algorithm, scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("Null" = "gray70", "Real" = "#E69F00")) +
  labs(title = "Performance: Real vs Null Models",
       subtitle = paste0("50 Replicates | HR: ", hr_value, " km"),
       x = "Multiplier (x HR)",
       y = "Metric Score") +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "Null_vs_Real_Comparison.png"), p_comp, width = 12, height = 8, dpi = 300)

message("Benchmark Complete. Check 'results/full_benchmark_v3' folder.")