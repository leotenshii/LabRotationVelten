# Last changes on 28.08.2024
# Author: Leoni Zimmermann

# Description ------------------------------------------------------------------
# This set of functions is designed to spped up the beginning of most of my script
# normalize_and_select_features log normalizes the data and selects features on variance
# mofa_build_model builds the mofa model and adds the metadata
# mofa_parameter_train trains the mofa model. It can also return the time needed for training
# stab_build_model builds the list used for StabMap calculations

# Functions --------------------------------------------------------------------

normalize_and_select_features <- function(data, expr_mean_threshold, p_value_threshold) {
  # The function log normalizes the input data
  data <- logNormCounts(data)
  decomp <- modelGeneVar(data)
  hvgs <- rownames(decomp)[decomp$mean > expr_mean_threshold & decomp$p.value <= p_value_threshold]
  data <- data[hvgs,]
  return(data)
}


mofa_build_model <- function(data, na_features, na_cells, metadata) {

  logcounts_allNA <- data
  logcounts_allNA[na_features, na_cells] <- NA
  
  mofa_list <- list(
    missing_feat = logcounts_allNA[na_features, ],
    all_feat = logcounts_allNA[setdiff(1:nrow(logcounts_allNA), na_features), ]
  )
  
  model <- create_mofa(mofa_list)
  
  metadata_df <- as.data.frame(metadata)
  if (has_rownames(metadata_df)) {
    metadata_df <- rownames_to_column(metadata_df, "sample")
  }
  metadata_df <- metadata_df %>%
    mutate(isNA = if_else(sample %in% colnames(data)[na_cells], "NA", "notNA"))
  

  samples_metadata(model) <- metadata_df
  
  return(model)
}


mofa_parameter_train <- function(scale_views = FALSE, num_factors = 15, spikeslab_weights = TRUE, model, save = TRUE, path = NULL, return_time = FALSE) {
  data_opts <- get_default_data_options(model)
  model_opts <- get_default_model_options(model)
  
  data_opts$scale_views <- scale_views
  model_opts$num_factors <- num_factors
  model_opts$spikeslab_weights <- spikeslab_weights
  
  MOFAobject <- prepare_mofa(model, 
                             data_options = data_opts,
                             model_options = model_opts)
  

  time_integration <- system.time(trained_model <- run_mofa(MOFAobject, save_data = save, outfile = path, use_basilisk = TRUE))
  
  if (return_time) {
    return(list(time_integration, trained_model))
  }
  else {
    return(trained_model)
  }
}


stab_build_model <- function(data, na_features, na_cells) {
  
  stab_list <- list(
    missing_feat = data[setdiff(1:nrow(data), na_features), na_cells],
    all_feat = data[, setdiff(1:ncol(data), na_cells)]
  )
  
  return(stab_list)
}

