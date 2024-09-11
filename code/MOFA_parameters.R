# Last changes on 28.08.2024
# Author: Leoni Zimmermann

#---------------------------Description-----------------------------------------
# This script runs MOFA with a range of different parameters.

#---------------------------Seed------------------------------------------------
set.seed(42)

#---------------------------Prep for running script-----------------------------
#---------------------------Libraries-------------------------------------------
suppressMessages(c(library(scater),
                   library(scran),
                   library(MOFA2),
                   library(tidyverse),
                   library(cluster),
                   library(patchwork),
                   library(fpc),
                   library(MultiAssayExperiment)))


source("~/R/Functions/adaptiveKNN.R")
source("~/R/Functions/data_prep_functions.R")
source("~/R/Functions/integration_metrics_functions.R")

#---------------------------Dataset---------------------------------------------
# Peripheral Blood Mononuclear Cells provided by 10x Genomics website
# 10x Genomics Multiome technology enables simultaneous profiling of the transcriptome 
# (using 3’ gene expression) and epigenome (using ATAC-seq) from single cells to deepen 
# our understanding of how genes are expressed and regulated across different cell types.

# Loaded from RDS cause I cant install the SingleCellMultiModal package
mae <- readRDS("/home/hd/hd_hd/hd_fb235/R/Data/data.RDS")
metadata <- mae@colData
#---------------------------Data preparation------------------------------------
# RNA
sce.rna <- normalize_and_select_features(experiments(mae)[["rna"]], 0.01, 0.05)

# ATAC
sce.atac <- normalize_and_select_features(experiments(mae)[["atac"]], 0.25, 0.05)

logcounts_all <- as.matrix(rbind(logcounts(sce.rna), logcounts(sce.atac)))
#---------------------------General Data----------------------------------------
# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all)))
colnames(all_celltypes) <- "celltype"

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

na_features <-  953:1740
na_cells <- 1:5016
#---------------------------MOFA------------------------------------------------
# Create model for MOFA
model <- mofa_build_model(data = logcounts_all, 
                 na_features = na_features,
                 na_cells = na_cells,
                 metadata = metadata)

# Get differnt model settings
data_opts <- get_default_data_options(model)
model_opts <- get_default_model_options(model)


#---------------------------Parameters------------------------------------------


param_grid <- expand.grid(
  scale_views = c(TRUE, FALSE),
  num_factors = c(10, 15, 20, 30, 40, 50, 60, 70),
  # spikeslab_factors = c(TRUE, FALSE),
  spikeslab_weights = c(TRUE, FALSE)
)

# param_grid <- expand.grid(
#   scale_views = c(FALSE),
#   num_factors = c(10, 15, 20, 30, 40, 50, 60, 70),
#   spikeslab_weights = c(FALSE)
# )




# Overview of parameters:
# scale_groups: if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE -> only have one group anyway
# scale_views: if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is FALSE
# num_factors: number of factors (In other tasks, such as imputation of missing values, even small sources of variation can be important and hence models should be trained with a large number of factors.)
# likelihoods: likelihood per view (options are ?gaussian?, ?poisson?, ?bernoulli?). By default they are learnt automatically. We advise users to use ?gaussian? whenever possible! -> kept the automatic
    
# Not sure about those:
# spikeslab_factors: use spike-slab sparsity prior in the factors? default is FALSE. -> all factors treated equally, if TRUE: some factors are exactly zero, simpler model with fewer active factors
# spikeslab_weights: use spike-slab sparsity prior in the weights? default is TRUE. -> some weights will be exactly zero, which helps in identifying which features are most relevant for each factor, if FALSE: all weights considered, more complex model
# ard_factors: use ARD prior in the factors? Default is TRUE if using multiple groups. -> automatically determining the relevance of each factor, potentially turning off irrelevant factors, if FALSE: all equal without relevance weighting
# ard_weights: use ARD prior in the weights? Default is TRUE if using multiple views. -> automatically determining the relevance of each feature of each view, potentially turning off irrelevant features, if FALSE: all equal without relevance weighting



# maxiter: number of iterations. Default is 1000. -> Kept like this
# convergence_mode: ?fast?, ?medium?, ?slow?. For exploration, the fast mode is good enough. -> Did fast
# startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence). -> See no reason to change
# freqELBO: frequency of computations of the ELBO. -> See no reason to change
# stochastic: use stochastic inference? (default is FALSE). -> TRUE would make it faster but less accurate

 




#---------------------------Result list-----------------------------------------
results <- list()

#---------------------------Loop------------------------------------------------
for (i in 1:nrow(param_grid)) {
  set.seed(42)
  params <- param_grid[i, ]
  
  # Print current parameters 

  cat("Running MOFA with scale_views =", params$scale_views,
      ", num_factors =", params$num_factors,
      # ", spikeslab_factors =", params$spikeslab_factors,
      ", spikeslab_weights =", params$spikeslab_weights)
  
  mofa_results <- mofa_parameter_train(scale_views = params$scale_views, num_factors = params$num_factors, spikeslab_weights = params$spikeslab_weights, model, return_time = TRUE, path = "/home/hd/hd_hd/hd_fb235/R/Data/model_par.hdf5")
  trained_model <- load_model("/home/hd/hd_hd/hd_fb235/R/Data/model_par.hdf5", remove_inactive_factors = FALSE)
  
  # Run MOFA with current parameters
  time_integration_all <- mofa_results[[1]]

  

  # UMAP
  trained_model <- run_umap(trained_model)
  mofa_umap_coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
  
  mofa_umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0)  %>%
    full_join(cluster) %>% 
    mutate(isNA =  c(rep(1, ncol(logcounts_all)/2), rep(0, ncol(logcounts_all) - ncol(logcounts_all)/2))) %>%  # 1 is TRUE, 0 is FALSE
    column_to_rownames("Row.names") %>%
    select(-sample)
  
  # Silhouette scores
  mofa_sil_sum <- silhouette_summary(mofa_umap$n, mofa_umap %>% select(UMAP1, UMAP2))
  
  # Clustering stats
  mofa_cluster_stats <- cluster.stats(dist(mofa_umap %>% select(UMAP1, UMAP2)), mofa_umap$n)
  dunn <- mofa_cluster_stats$dunn
  
  # Cell Type Accuracy
  rna_train <- mofa_umap$celltype[na_cells]
  names(rna_train) <- rownames(mofa_umap)[na_cells]
  atac_query <- mofa_umap$celltype[setdiff(1:ncol(logcounts_all), na_cells)]
  names(atac_query) <- rownames(mofa_umap)[setdiff(1:ncol(logcounts_all), na_cells)]
  
  mofa_knn_out = embeddingKNN(mofa_umap_coord,
                               rna_train,
                               type = "uniform_fixed",
                               k_values = 5)
  mofa_knn_acc = mean(isEqual(mofa_knn_out[names(atac_query),"predicted_labels"], atac_query), na.rm = TRUE)
  mofa_knn_acc_bal = mean(unlist(lapply(split(isEqual(mofa_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
  
  # RMSE
  imputation_time_all <- system.time(trained_model <- impute(trained_model))
  
  mofa_rmse<- rmse_imp(imp_data = trained_model@imputed_data$missing_feat[[1]][,na_cells],
                       real_data = logcounts_all, 
                       na_features = na_features, 
                       na_cells = na_cells)
  
  # Store results

  param_str <- paste(
    "scale_views", params$scale_views,
    "num_factor", params$num_factor,
    # "spikeslab_factors", params$spikeslab_factors,
    "spikeslab_weights", params$spikeslab_weights,
    sep = "_"
  )
  

  results[[param_str]] <- list(
    mofa_umap = mofa_umap,
    mofa_sil_sum = mofa_sil_sum,
    dunn = dunn,
    mofa_cluster_stats = mofa_cluster_stats,
    mofa_knn_acc = mofa_knn_acc,
    mofa_knn_acc_bal = mofa_knn_acc_bal,
    mofa_rmse = mofa_rmse,
    time_integration = time_integration_all,
    imputation_time = imputation_time_all
  )
  
  
}

#---------------------------Comparison results----------------------------------
# Create data frame

comparison <- data.frame(
  scale_views = logical(),
  num_factor = integer(),
  # spikeslab_factors = logical(),
  spikeslab_weights = logical(),
  mean_sil_score = numeric(),
  min_sil_score = numeric(),
  max_sil_score = numeric(),
  dunn = numeric(),
  mofa_knn_acc = numeric(),
  mofa_knn_acc_bal = numeric(),
  mofa_rmse = numeric(),
  time_integration = numeric(),
  imputation_time = numeric()
)

# Fill data frame
for (param_name in names(results)) {
  res <- results[[param_name]]
  mean_sil_score <- mean(res$mofa_sil_sum$score, na.rm = TRUE)
  min_sil_score <- min(res$mofa_sil_sum$score,na.rm = TRUE)
  max_sil_score <- max(res$mofa_sil_sum$score,na.rm = TRUE)
  dunn <- res$dunn
  mofa_knn_acc <- res$mofa_knn_acc
  mofa_knn_acc_bal <- res$mofa_knn_acc_bal
  mofa_rmse <- res$mofa_rmse
  time_integration <- res$time_integration[3]
  imputation_time <- res$imputation_time[3]
  
  param_values <- unlist(strsplit(param_name, "_"))

  comparison <- rbind(comparison, data.frame(
    scale_views = as.logical(param_values[3]),
    num_factor = as.numeric(param_values[6]),
    # spikeslab_factors = as.logical(param_values[9]),
    spikeslab_weights = as.logical(param_values[9]),
    mean_sil_score = mean_sil_score,
    min_sil_score = min_sil_score,
    max_sil_score = max_sil_score,
    dunn = dunn,
    mofa_knn_acc = mofa_knn_acc,
    mofa_knn_acc_bal = mofa_knn_acc_bal,
    mofa_rmse = mofa_rmse,
    time_integration = time_integration,
    imputation_time = imputation_time
  ))

}

print(comparison)

 write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/mofa_comparison1.txt")

# For this file:
# write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/mofa_comparison_numfactors2.txt")
# change to:
# param_grid <- expand.grid(
#   scale_views = c(FALSE),
#   num_factors = c(10, 15, 20, 30, 40, 50, 60, 70),
#   spikeslab_weights = c(FALSE)
# )






