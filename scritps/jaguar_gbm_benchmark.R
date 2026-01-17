# GBM SDM ANALYSIS WITH NULL MODEL COMPARISON
# Species Distribution Modeling with hyperparameter tuning and validation
# Author: Oliveira, 2025
# Date: 2026-01-16

library(tidyverse)
library(terra)
library(flexsdm)
library(GeoThinneR)

# Configuration
presence_file <- "data/presence/jaguar_data.csv"
background_file <- "data/background/jaguar_bg_nnet_gbm_rf.csv"
rasters_path <- "rasters/variables"
output_dir <- "results/gbm_benk/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed_base <- 18
base_thin_distance <- 9
thin_multipliers <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)
#thin_multipliers <- c(0,1) #test only
n_cores_geothinner <- 4
n_cores_flexsdm <- 8
n_replicates <- 50

part_method <- c(method = "rep_kfold", folds = 4, replicates = 5)

gbm_grid <- expand.grid(
  n.trees = c(20, 50, 100),
  shrinkage = c(0.1, 0.5, 1),
  n.minobsinnode = c(1, 3, 5, 7, 9),
  stringsAsFactors = FALSE
)

threshold_type <- "max_sens_spec"
selection_metric <- "BOYCE"

# Load data
occ_raw <- read_csv(presence_file, show_col_types = FALSE) %>%
  rename(x = Longitude, y = Latitude) %>%
  mutate(pr_ab = 1) %>%
  select(x, y, pr_ab)

bg_raw <- read_csv(background_file, show_col_types = FALSE) %>%
  rename(x = Longitude, y = Latitude) %>%
  mutate(pr_ab = 0) %>%
  select(x, y, pr_ab)

var_rasters <- list.files(rasters_path, pattern = "\\.tif$", full.names = TRUE) %>%
  rast()

cat("Data loaded: ", nrow(occ_raw), "presences | ", nrow(bg_raw), "background points | ",
    nlyr(var_rasters), "variables\n\n")

# Initialize results storage
all_results <- list()
result_idx <- 1

# Main analysis loop
for (mult_idx in seq_along(thin_multipliers)) {
  
  mult <- thin_multipliers[mult_idx]
  thin_distance <- base_thin_distance * mult
  current_seed <- seed_base + mult_idx
  
  cat("Multiplier", mult, "(", thin_distance, "km ) - Seed:", current_seed, "\n")
  
  # Spatial thinning
  if (thin_distance > 0) {
    real_thin <- thin_points(
      data = occ_raw,
      lon_col = "x",
      lat_col = "y",
      method = "distance",
      thin_dist = thin_distance,
      search_type = "local_kd_tree",
      trials = 10,
      n_cores = n_cores_geothinner,
      seed = current_seed,
      verbose = FALSE
    )
    # Filter original data using retained logical vector
    occ_thinned <- real_thin$original_data[real_thin$retained[[1]], ]
  } else {
    occ_thinned <- occ_raw
  }
  
  # Extract environmental values
  swd_data <- sdm_extract(
    data = occ_thinned,
    x = "x",
    y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
  bg_data <- sdm_extract(
    data = bg_raw,
    x = "x",
    y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
  n_pres <- sum(swd_data$pr_ab == 1)
  n_bg <- sum(bg_data$pr_ab == 0)
  
  cat("  Thinned:", nrow(occ_thinned), "| Extracted:", n_pres, "presences\n")
  
  # Partition presences
  set.seed(current_seed)
  presence_part <- part_random(
    data = swd_data,
    pr_ab = "pr_ab",
    method = part_method
  )
  
  part_cols <- grep("^\\.part", names(presence_part), value = TRUE)
  predictor_names <- setdiff(names(presence_part), 
                             c("x", "y", "pr_ab", part_cols))
  
  # Real model replicates
  cat("  Real models: ")
  real_metrics <- list()
  
  for (rep in 1:n_replicates) {
    rep_seed <- current_seed + rep * 100
    
    part_col <- part_cols[((rep - 1) %% length(part_cols)) + 1]
    
    presence_data <- presence_part %>%
      mutate(.part = .data[[part_col]]) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    fold_ids <- sort(unique(presence_data$.part))
    set.seed(rep_seed)
    bg_part <- bg_data %>%
      mutate(.part = sample(fold_ids, size = n(), replace = TRUE)) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    # Combine presence and background data for tune_gbm
    combined_data <- bind_rows(presence_data, bg_part)
    
    model <- tune_gbm(
      data = combined_data,
      response = "pr_ab",
      predictors = predictor_names,
      partition = ".part",
      grid = gbm_grid,
      thr = threshold_type,
      metric = selection_metric,
      n_cores = n_cores_flexsdm
    )
    
    real_metrics[[rep]] <- model$performance %>%
      mutate(
        replicate = rep,
        model_type = "Real",
        multiplier = mult
      ) %>%
      select(replicate, model_type, multiplier, n.trees, shrinkage, n.minobsinnode,
             TSS = TSS_mean, AUC = AUC_mean, BOYCE = BOYCE_mean, OR = OR_mean, MCC = MCC_mean)
    
    if (rep %% 10 == 0) cat(rep, "")
  }
  cat("done\n")
  
  # Null model replicates
  # Skip null models when multiplier = 0 (no thinning)
  if (mult == 0) {
    cat("  Null models: Skipped (multiplier = 0)\n\n")
  } else {
    cat("  Null models: ")
    null_metrics <- list()
    
    best_n.trees <- real_metrics[[1]]$n.trees
    best_shrinkage <- real_metrics[[1]]$shrinkage
    best_n.minobsinnode <- real_metrics[[1]]$n.minobsinnode
    
    for (rep in 1:n_replicates) {
      rep_seed <- current_seed + rep * 100 + 5000
      
      part_col <- part_cols[((rep - 1) %% length(part_cols)) + 1]
      
      presence_template <- presence_part %>%
        mutate(.part = .data[[part_col]]) %>%
        select(x, y, pr_ab, .part, all_of(predictor_names))
      
      fold_ids <- sort(unique(presence_template$.part))
      set.seed(rep_seed)
      
      # Null model logic depends on multiplier
      if (mult == 0) {
        # No thinning: use random background points as pseudo-presences
        null_presences <- bg_data %>%
          sample_n(size = n_pres, replace = FALSE) %>%
          mutate(
            pr_ab = 1,
            .part = sample(fold_ids, n_pres, replace = TRUE)
          ) %>%
          select(x, y, pr_ab, .part, all_of(predictor_names))
      } else {
        # With thinning: randomly sample from real presences (before thinning)
        # This maintains spatial sampling structure but randomizes which presences
        null_presences <- swd_data %>%
          sample_n(size = n_pres, replace = FALSE) %>%
          mutate(
            pr_ab = 1,
            .part = sample(fold_ids, n_pres, replace = TRUE)
          ) %>%
          select(x, y, pr_ab, .part, all_of(predictor_names))
      }
      
      set.seed(rep_seed + 1)
      bg_part <- bg_data %>%
        mutate(.part = sample(fold_ids, size = n(), replace = TRUE)) %>%
        select(x, y, pr_ab, .part, all_of(predictor_names))
      
      # Combine null presences and background data for tune_gbm
      combined_null_data <- bind_rows(null_presences, bg_part)
      
      null_model <- tryCatch({
        tune_gbm(
          data = combined_null_data,
          response = "pr_ab",
          predictors = predictor_names,
          partition = ".part",
          grid = data.frame(
            n.trees = best_n.trees,
            shrinkage = best_shrinkage,
            n.minobsinnode = best_n.minobsinnode,
            stringsAsFactors = FALSE
          ),
          thr = threshold_type,
          metric = selection_metric,
          n_cores = n_cores_flexsdm
        )
      }, error = function(e) {
        cat("X")  # Mark failed replicate
        return(NULL)
      })
      
      if (!is.null(null_model)) {
        null_metrics[[rep]] <- null_model$performance %>%
          mutate(
            replicate = rep,
            model_type = "Null",
            multiplier = mult
          ) %>%
          select(replicate, model_type, multiplier, n.trees, shrinkage, n.minobsinnode,
                 TSS = TSS_mean, AUC = AUC_mean, BOYCE = BOYCE_mean, OR = OR_mean, MCC = MCC_mean)
      } else {
        # Store NA results for failed replicates
        null_metrics[[rep]] <- tibble(
          replicate = rep,
          model_type = "Null",
          multiplier = mult,
          n.trees = best_n.trees,
          shrinkage = best_shrinkage,
          n.minobsinnode = best_n.minobsinnode,
          TSS = NA_real_,
          AUC = NA_real_,
          BOYCE = NA_real_,
          OR = NA_real_,
          MCC = NA_real_
        )
      }
      
      if (rep %% 10 == 0) cat(rep, "")
    }
    cat("done\n\n")
  } # End of else block for null models
  
  # Combine results for this multiplier
  if (mult == 0) {
    all_results[[result_idx]] <- bind_rows(real_metrics)
  } else {
    all_results[[result_idx]] <- bind_rows(real_metrics, null_metrics)
  }
  result_idx <- result_idx + 1
}

# Compile final results
final_results <- bind_rows(all_results)

# Save complete results
write_csv(final_results, paste0(output_dir, "gbm_results_complete.csv"))

# Summary statistics
summary_stats <- final_results %>%
  group_by(multiplier, model_type) %>%
  summarise(
    n_replicates = n(),
    TSS_mean = mean(TSS, na.rm = TRUE),
    TSS_sd = sd(TSS, na.rm = TRUE),
    AUC_mean = mean(AUC, na.rm = TRUE),
    AUC_sd = sd(AUC, na.rm = TRUE),
    BOYCE_mean = mean(BOYCE, na.rm = TRUE),
    BOYCE_sd = sd(BOYCE, na.rm = TRUE),
    OR_mean = mean(OR, na.rm = TRUE),
    OR_sd = sd(OR, na.rm = TRUE),
    MCC_mean = mean(MCC, na.rm = TRUE),
    MCC_sd = sd(MCC, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(summary_stats, paste0(output_dir, "gbm_results_summary.csv"))

cat("Results summary:\n")
print(summary_stats, n = Inf)

# Visualization
plot_data <- final_results %>%
  pivot_longer(
    cols = c(TSS, AUC, BOYCE, OR, MCC),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    multiplier = factor(multiplier),
    metric = factor(metric, levels = c("TSS", "AUC", "BOYCE", "OR", "MCC"))
  )

p <- ggplot(plot_data, aes(x = multiplier, y = value, fill = model_type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Real" = "#2E86AB", "Null" = "#A23B72")) +
  labs(
    title = "Model Performance Comparison: Real vs Null Models (GBM)",
    subtitle = paste0("n = ", n_replicates, " replicates per multiplier"),
    x = "Thinning Multiplier",
    y = "Metric Value",
    fill = "Model Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(paste0(output_dir, "performance_comparison_boxplots.png"), 
       p, width = 12, height = 10, dpi = 300)

cat("\nAnalysis complete. Results saved to:", output_dir, "\n")