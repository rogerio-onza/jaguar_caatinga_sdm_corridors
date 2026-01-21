# Unified SDM workflow: Tuning → Selection → Pause → Ensemble → MSDM
# Author: Oliveira, 2025
# Date: 2026-01-20
# Version: 1.1

library(tidyverse)
library(terra)
library(flexsdm)
library(GeoThinneR)

# Global configuration
seed_base <- 18
n_cores_geothinner <- 4
n_cores_flexsdm <- 8
part_method <- c(method = "rep_kfold", folds = 4, replicates = 50)
threshold_type <- "max_sens_spec"
selection_metric <- "BOYCE"

base_thin_distance <- 9
thin_multipliers <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)

presence_file <- "data/presence/jaguar_data.csv"
background_nnet_gbm_raf <- "data/background/jaguar_bg_nnet_gbm_rf.csv"
background_maxent <- "data/background/jaguar_bg_maxe.csv"
rasters_path <- "rasters/variables"

output_dirs <- list(
  net = "results/net/",
  gbm = "results/gbm/",
  raf = "results/raf/",
  maxent = "results/maxent/",
  selection = "results/selection/",
  ensemble = "results/ensemble/"
)

lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

# Hyperparameter grids
net_grid <- expand.grid(
  size = c(2, 4, 6, 8, 10),
  decay = c(0.001, 0.05, 0.1, 1, 3, 4, 5, 10),
  stringsAsFactors = FALSE
)

gbm_grid <- expand.grid(
  n.trees = c(20, 50, 100),
  shrinkage = c(0.1, 0.5, 1),
  n.minobsinnode = c(1, 3, 5, 7, 9),
  stringsAsFactors = FALSE
)

raf_grid <- expand.grid(
  mtry = c(1, 2, 3, 4, 5, 6, 7),
  ntree = c(500, 1000),
  stringsAsFactors = FALSE
)

maxent_grid <- expand.grid(
  regmult = c(0.5, 1, 1.5, 2, 3, 4),
  classes = c("l", "lq", "lqpht"),
  stringsAsFactors = FALSE
)

# Load data
occ_raw <- read_csv(presence_file, show_col_types = FALSE) %>%
  rename(x = Longitude, y = Latitude) %>%
  mutate(pr_ab = 1) %>%
  select(x, y, pr_ab)

bg_nnet_gbm_raf_raw <- read_csv(background_nnet_gbm_raf, show_col_types = FALSE) %>%
  rename(x = Longitude, y = Latitude) %>%
  mutate(pr_ab = 0) %>%
  select(x, y, pr_ab)

bg_maxent_raw <- read_csv(background_maxent, show_col_types = FALSE) %>%
  rename(x = Longitude, y = Latitude) %>%
  mutate(pr_ab = 0) %>%
  select(x, y, pr_ab)

var_rasters <- list.files(rasters_path, pattern = "\\.tif$", full.names = TRUE) %>%
  rast()

# Phase 1A: Tuning - Neural Network
n_replicates <- 50
net_all_results <- list()
result_idx <- 1

for (mult_idx in seq_along(thin_multipliers)) {
  
  mult <- thin_multipliers[mult_idx]
  thin_distance <- base_thin_distance * mult
  current_seed <- seed_base + mult_idx
  
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
    occ_thinned <- real_thin$original_data[real_thin$retained[[1]], ]
  } else {
    occ_thinned <- occ_raw
  }
  
  # Extract environmental values
  swd_data <- sdm_extract(
    data = occ_thinned,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
  bg_data <- sdm_extract(
    data = bg_nnet_gbm_raf_raw,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
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
  
  # Run 50 replicates for this multiplier
  for (rep in 1:n_replicates) {
    
    rep_seed <- current_seed + rep * 100
    
    # Cycle through partition schemes
    part_col <- part_cols[((rep - 1) %% length(part_cols)) + 1]
    
    # Prepare data with selected partition scheme
    presence_data <- presence_part %>%
      mutate(.part = .data[[part_col]]) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    fold_ids <- sort(unique(presence_data$.part))
    set.seed(rep_seed)
    bg_part <- bg_data %>%
      mutate(.part = sample(fold_ids, size = n(), replace = TRUE)) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    combined_data <- bind_rows(presence_data, bg_part)
    
    # Tune model
    net_model <- tune_net(
      data = combined_data,
      response = "pr_ab",
      predictors = predictor_names,
      partition = ".part",
      grid = net_grid,
      thr = threshold_type,
      metric = selection_metric,
      n_cores = n_cores_flexsdm
    )
    
    # Store results
    net_all_results[[result_idx]] <- list(
      replicate = rep,
      multiplier = mult,
      model_object = net_model,
      performance = net_model$performance %>% 
        mutate(replicate = rep, multiplier = mult)
    )
    
    result_idx <- result_idx + 1
  }
}

# Save Neural Network results
saveRDS(net_all_results, paste0(output_dirs$net, "net_tuning_complete.rds"))
net_performance_all <- bind_rows(lapply(net_all_results, function(x) x$performance))
write_csv(net_performance_all, paste0(output_dirs$net, "net_performance_all.csv"))

# Phase 1B: Tuning - GBM
gbm_all_results <- list()
result_idx <- 1

for (mult_idx in seq_along(thin_multipliers)) {
  
  mult <- thin_multipliers[mult_idx]
  thin_distance <- base_thin_distance * mult
  current_seed <- seed_base + mult_idx
  
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
    occ_thinned <- real_thin$original_data[real_thin$retained[[1]], ]
  } else {
    occ_thinned <- occ_raw
  }
  
  # Extract environmental values
  swd_data <- sdm_extract(
    data = occ_thinned,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
  bg_data <- sdm_extract(
    data = bg_nnet_gbm_raf_raw,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
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
  
  # Run 50 replicates for this multiplier
  for (rep in 1:n_replicates) {
    
    rep_seed <- current_seed + rep * 100
    
    # Cycle through partition schemes
    part_col <- part_cols[((rep - 1) %% length(part_cols)) + 1]
    
    # Prepare data with selected partition scheme
    presence_data <- presence_part %>%
      mutate(.part = .data[[part_col]]) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    fold_ids <- sort(unique(presence_data$.part))
    set.seed(rep_seed)
    bg_part <- bg_data %>%
      mutate(.part = sample(fold_ids, size = n(), replace = TRUE)) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    combined_data <- bind_rows(presence_data, bg_part)
    
    # Tune model
    gbm_model <- tune_gbm(
      data = combined_data,
      response = "pr_ab",
      predictors = predictor_names,
      partition = ".part",
      grid = gbm_grid,
      thr = threshold_type,
      metric = selection_metric,
      n_cores = n_cores_flexsdm
    )
    
    # Store results
    gbm_all_results[[result_idx]] <- list(
      replicate = rep,
      multiplier = mult,
      model_object = gbm_model,
      performance = gbm_model$performance %>% 
        mutate(replicate = rep, multiplier = mult)
    )
    
    result_idx <- result_idx + 1
  }
}

# Save GBM results
saveRDS(gbm_all_results, paste0(output_dirs$gbm, "gbm_tuning_complete.rds"))
gbm_performance_all <- bind_rows(lapply(gbm_all_results, function(x) x$performance))
write_csv(gbm_performance_all, paste0(output_dirs$gbm, "gbm_performance_all.csv"))

# Phase 1C: Tuning - Random Forest
raf_all_results <- list()
result_idx <- 1

for (mult_idx in seq_along(thin_multipliers)) {
  
  mult <- thin_multipliers[mult_idx]
  thin_distance <- base_thin_distance * mult
  current_seed <- seed_base + mult_idx
  
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
    occ_thinned <- real_thin$original_data[real_thin$retained[[1]], ]
  } else {
    occ_thinned <- occ_raw
  }
  
  # Extract environmental values
  swd_data <- sdm_extract(
    data = occ_thinned,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
  bg_data <- sdm_extract(
    data = bg_nnet_gbm_raf_raw,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
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
  
  # Run 50 replicates for this multiplier
  for (rep in 1:n_replicates) {
    
    rep_seed <- current_seed + rep * 100
    
    # Cycle through partition schemes
    part_col <- part_cols[((rep - 1) %% length(part_cols)) + 1]
    
    # Prepare data with selected partition scheme
    presence_data <- presence_part %>%
      mutate(.part = .data[[part_col]]) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    fold_ids <- sort(unique(presence_data$.part))
    set.seed(rep_seed)
    bg_part <- bg_data %>%
      mutate(.part = sample(fold_ids, size = n(), replace = TRUE)) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    combined_data <- bind_rows(presence_data, bg_part)
    
    # Tune model
    raf_model <- tune_raf(
      data = combined_data,
      response = "pr_ab",
      predictors = predictor_names,
      partition = ".part",
      grid = raf_grid,
      thr = threshold_type,
      metric = selection_metric,
      n_cores = n_cores_flexsdm
    )
    
    # Store results
    raf_all_results[[result_idx]] <- list(
      replicate = rep,
      multiplier = mult,
      model_object = raf_model,
      performance = raf_model$performance %>% 
        mutate(replicate = rep, multiplier = mult)
    )
    
    result_idx <- result_idx + 1
  }
}

# Save Random Forest results
saveRDS(raf_all_results, paste0(output_dirs$raf, "raf_tuning_complete.rds"))
raf_performance_all <- bind_rows(lapply(raf_all_results, function(x) x$performance))
write_csv(raf_performance_all, paste0(output_dirs$raf, "raf_performance_all.csv"))

# Phase 1D: Tuning - MaxEnt
maxent_all_results <- list()
result_idx <- 1

for (mult_idx in seq_along(thin_multipliers)) {
  
  mult <- thin_multipliers[mult_idx]
  thin_distance <- base_thin_distance * mult
  current_seed <- seed_base + mult_idx
  
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
    occ_thinned <- real_thin$original_data[real_thin$retained[[1]], ]
  } else {
    occ_thinned <- occ_raw
  }
  
  # Extract environmental values
  swd_data <- sdm_extract(
    data = occ_thinned,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
  bg_data <- sdm_extract(
    data = bg_maxent_raw,
    x = "x", y = "y",
    env_layer = var_rasters,
    filter_na = TRUE
  )
  
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
  
  # Run 50 replicates for this multiplier
  for (rep in 1:n_replicates) {
    
    rep_seed <- current_seed + rep * 100
    
    # Cycle through partition schemes
    part_col <- part_cols[((rep - 1) %% length(part_cols)) + 1]
    
    # Prepare data with selected partition scheme
    presence_data <- presence_part %>%
      mutate(.part = .data[[part_col]]) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    fold_ids <- sort(unique(presence_data$.part))
    set.seed(rep_seed)
    bg_part <- bg_data %>%
      mutate(.part = sample(fold_ids, size = n(), replace = TRUE)) %>%
      select(x, y, pr_ab, .part, all_of(predictor_names))
    
    # Tune model
    maxent_model <- tune_max(
      data = presence_data,
      response = "pr_ab",
      predictors = predictor_names,
      partition = ".part",
      background = bg_part,
      grid = maxent_grid,
      thr = threshold_type,
      metric = selection_metric,
      clamp = TRUE,
      pred_type = "cloglog",
      n_cores = n_cores_flexsdm
    )
    
    # Store results
    maxent_all_results[[result_idx]] <- list(
      replicate = rep,
      multiplier = mult,
      model_object = maxent_model,
      performance = maxent_model$performance %>% 
        mutate(replicate = rep, multiplier = mult)
    )
    
    result_idx <- result_idx + 1
  }
}

# Save MaxEnt results
saveRDS(maxent_all_results, paste0(output_dirs$maxent, "maxent_tuning_complete.rds"))
maxent_performance_all <- bind_rows(lapply(maxent_all_results, function(x) x$performance))
write_csv(maxent_performance_all, paste0(output_dirs$maxent, "maxent_performance_all.csv"))

# Phase 2: Model selection via Composite Score (3 metrics: TSS, BOYCE, OR)
net_selection_data <- net_performance_all %>%
  mutate(
    model = "Neural Network",
    hyperparams = paste0("size=", size, ", decay=", decay)
  ) %>%
  select(model, multiplier, hyperparams, size, decay, 
         TSS = TSS_mean, BOYCE = BOYCE_mean, OR = OR_mean,
         everything())

gbm_selection_data <- gbm_performance_all %>%
  mutate(
    model = "GBM",
    hyperparams = paste0("n.trees=", n.trees, ", shrinkage=", shrinkage, 
                         ", n.minobsinnode=", n.minobsinnode)
  ) %>%
  select(model, multiplier, hyperparams, n.trees, shrinkage, n.minobsinnode,
         TSS = TSS_mean, BOYCE = BOYCE_mean, OR = OR_mean,
         everything())

raf_selection_data <- raf_performance_all %>%
  mutate(
    model = "Random Forest",
    hyperparams = paste0("mtry=", mtry, ", ntree=", ntree)
  ) %>%
  select(model, multiplier, hyperparams, mtry, ntree,
         TSS = TSS_mean, BOYCE = BOYCE_mean, OR = OR_mean,
         everything())

maxent_selection_data <- maxent_performance_all %>%
  mutate(
    model = "MaxEnt",
    hyperparams = paste0("regmult=", regmult, ", classes=", classes)
  ) %>%
  select(model, multiplier, hyperparams, regmult, classes,
         TSS = TSS_mean, BOYCE = BOYCE_mean, OR = OR_mean,
         everything())

# Combine all data
all_selection_data <- bind_rows(
  net_selection_data,
  gbm_selection_data,
  raf_selection_data,
  maxent_selection_data
)

# Calculate composite scores per model
composite_scores <- all_selection_data %>%
  group_by(model) %>%
  mutate(
    TSS_norm = (TSS - min(TSS, na.rm = TRUE)) / 
      (max(TSS, na.rm = TRUE) - min(TSS, na.rm = TRUE)),
    BOYCE_norm = (BOYCE - min(BOYCE, na.rm = TRUE)) / 
      (max(BOYCE, na.rm = TRUE) - min(BOYCE, na.rm = TRUE)),
    OR_norm = 1 - ((OR - min(OR, na.rm = TRUE)) / 
                     (max(OR, na.rm = TRUE) - min(OR, na.rm = TRUE))),
    composite_score = (TSS_norm + BOYCE_norm + OR_norm) / 3
  ) %>%
  ungroup() %>%
  arrange(model, desc(composite_score))

# Save all composite scores
write_csv(composite_scores, paste0(output_dirs$selection, "all_composite_scores.csv"))

# Select best model for each algorithm
best_models <- composite_scores %>%
  group_by(model) %>%
  slice_max(order_by = composite_score, n = 1) %>%
  ungroup() %>%
  select(model, multiplier, hyperparams, TSS, BOYCE, OR, composite_score)

write_csv(best_models, paste0(output_dirs$selection, "best_models_selected.csv"))

# Mandatory pause for user validation
stop("PAUSE: Review 'results/selection/best_models_selected.csv' and re-run script to continue from Phase 3.")

# Phase 3: Extract selected models and prepare for ensemble
net_all_results <- readRDS(paste0(output_dirs$net, "net_tuning_complete.rds"))
gbm_all_results <- readRDS(paste0(output_dirs$gbm, "gbm_tuning_complete.rds"))
raf_all_results <- readRDS(paste0(output_dirs$raf, "raf_tuning_complete.rds"))
maxent_all_results <- readRDS(paste0(output_dirs$maxent, "maxent_tuning_complete.rds"))

best_models <- read_csv(paste0(output_dirs$selection, "best_models_selected.csv"), 
                        show_col_types = FALSE)

# Function to extract best model from results using exact hyperparameters
extract_selected_model <- function(all_results, best_row, model_type) {
  target_mult <- best_row$multiplier
  hyperparams_str <- best_row$hyperparams
  
  # Filter results for target multiplier
  mult_results <- all_results[sapply(all_results, function(x) x$multiplier == target_mult)]
  
  if (length(mult_results) == 0) {
    stop("Could not find results for multiplier ", target_mult)
  }
  
  # Find the model with matching hyperparameters
  for (i in seq_along(mult_results)) {
    result <- mult_results[[i]]
    perf <- result$performance
    
    # Build hyperparameter string based on model type
    if (model_type == "net") {
      current_hyper <- paste0("size=", perf$size, ", decay=", perf$decay)
    } else if (model_type == "gbm") {
      current_hyper <- paste0("n.trees=", perf$n.trees, 
                              ", shrinkage=", perf$shrinkage, 
                              ", n.minobsinnode=", perf$n.minobsinnode)
    } else if (model_type == "raf") {
      current_hyper <- paste0("mtry=", perf$mtry, ", ntree=", perf$ntree)
    } else if (model_type == "maxent") {
      current_hyper <- paste0("regmult=", perf$regmult, ", classes=", perf$classes)
    }
    
    # Check if hyperparameters match
    if (current_hyper == hyperparams_str) {
      return(result$model_object)
    }
  }
  
  stop("Could not find model with hyperparameters: ", hyperparams_str, 
       " for multiplier: ", target_mult)
}

# Extract each model
net_best <- extract_selected_model(
  net_all_results, 
  best_models %>% filter(model == "Neural Network"),
  model_type = "net"
)

gbm_best <- extract_selected_model(
  gbm_all_results,
  best_models %>% filter(model == "GBM"),
  model_type = "gbm"
)

raf_best <- extract_selected_model(
  raf_all_results,
  best_models %>% filter(model == "Random Forest"),
  model_type = "raf"
)

maxent_best <- extract_selected_model(
  maxent_all_results,
  best_models %>% filter(model == "MaxEnt"),
  model_type = "maxent"
)

# Save individual model performances
write_csv(net_best$performance, paste0(output_dirs$ensemble, "net_selected_performance.csv"))
write_csv(gbm_best$performance, paste0(output_dirs$ensemble, "gbm_selected_performance.csv"))
write_csv(raf_best$performance, paste0(output_dirs$ensemble, "raf_selected_performance.csv"))
write_csv(maxent_best$performance, paste0(output_dirs$ensemble, "maxent_selected_performance.csv"))


# Phase 4: Create ensemble model
ensemble_model <- fit_ensemble(
  models = list(net_best, gbm_best, raf_best, maxent_best),
  ens_method = "median",
  thr = threshold_type,
  thr_model = threshold_type,
  metric = selection_metric
)


# Save ensemble performance
write_csv(ensemble_model$performance, paste0(output_dirs$ensemble, "ensemble_performance.csv"))

# Create summary of all models
all_models_summary <- bind_rows(
  net_best$performance %>% mutate(model = "Neural Network"),
  gbm_best$performance %>% mutate(model = "GBM"),
  raf_best$performance %>% mutate(model = "Random Forest"),
  maxent_best$performance %>% mutate(model = "MaxEnt"),
  ensemble_model$performance %>% mutate(model = "Ensemble")
) %>%
  select(model, everything())

write_csv(all_models_summary, paste0(output_dirs$ensemble, "all_models_summary.csv"))

# Generate predictions
prediction <- sdm_predict(
  models = ensemble_model,
  pred = var_rasters,
  thr = threshold_type,
  con_thr = TRUE,
  clamp = TRUE,
  pred_type = "cloglog"
)

prediction_continuous <- prediction$median[[1]]
prediction_binary <- prediction$median[[2]]

plot(prediction_continuous)

names(prediction_continuous) <- "ensemble_continuous"
names(prediction_binary) <- "ensemble_binary"

writeRaster(prediction_continuous, 
            paste0(output_dirs$ensemble, "ensemble_continuous.tif"),
            overwrite = TRUE)
writeRaster(prediction_binary, 
            paste0(output_dirs$ensemble, "ensemble_binary.tif"),
            overwrite = TRUE)

# Phase 5

occ_thinned <- occ_raw  # no thinning

swd_data <- sdm_extract(
  data = occ_thinned,
  x = "x", y = "y",
  env_layer = var_rasters,
  filter_na = TRUE
)

bg_data <- sdm_extract(
  data = bg_nnet_gbm_raf_raw,
  x = "x", y = "y",
  env_layer = var_rasters,
  filter_na = TRUE
)

records_for_msdm <- bind_rows(
  swd_data %>% select(x, y, pr_ab),
  bg_data %>% select(x, y, pr_ab)
)

# Apply MSDM posteriori correction
msdm_result <- msdm_posteriori(
  records = records_for_msdm,
  x = "x",
  y = "y",
  pr_ab = "pr_ab",
  cont_suit = prediction_continuous,
  method = "obr",
  thr = 0.547,
  buffer = NULL
)

prediction_continuous_msdm <- msdm_result[[1]]
prediction_binary_msdm <- msdm_result[[2]]

plot(prediction_continuous_msdm)

names(prediction_continuous_msdm) <- "ensemble_continuous_msdm"
names(prediction_binary_msdm) <- "ensemble_binary_msdm"

writeRaster(prediction_continuous_msdm, 
            paste0(output_dirs$ensemble, "ensemble_continuous_msdm.tif"),
            overwrite = TRUE)
writeRaster(prediction_binary_msdm, 
            paste0(output_dirs$ensemble, "ensemble_binary_msdm.tif"),
            overwrite = TRUE)

# Save metadata
metadata <- list(
  date = Sys.time(),
  configuration = list(
    seed_base = seed_base,
    part_method = part_method,
    threshold_type = threshold_type,
    selection_metric = selection_metric,
    thin_multipliers = thin_multipliers
  ),
  selected_models = best_models,
  n_presences = sum(ensemble_model$data_ens$pr_ab == 1),
  n_background = sum(ensemble_model$data_ens$pr_ab == 0)
)

jsonlite::write_json(metadata, 
                     paste0(output_dirs$ensemble, "ensemble_metadata.json"),
                     pretty = TRUE,
                     auto_unbox = TRUE)

# End of workflow