# Last changes on 28.08.2024
# Author: Leoni Zimmermann

#---------------------------Description-----------------------------------------
# This script runs StabMap with a range of different parameters.

#---------------------------Libraries-------------------------------------------
suppressMessages(c(
  library(scater),
  library(scran),
  library(StabMap),
  library(tidyverse),
  library(cluster),
  library(patchwork),
  library(fpc),
  library(MultiAssayExperiment)))

source("~/R/Functions/adaptiveKNN.R")
source("~/R/Functions/data_prep_functions.R")
source("~/R/Functions/integration_metrics_functions.R")

set.seed(42)

#---------------------------Parameters------------------------------------------
na_features <-  953:1740
na_cells <- 1:5016

param_grid <- expand.grid(
  ncomponentsReference = c( 10, 15, 20, 30, 40, 50, 60, 70),
  ncomponentsSubset = c( 10, 15, 20, 30, 40, 50, 60, 70),
  maxFeatures = c(1740),
  scale.center = c(TRUE, FALSE),
  scale.scale = c(TRUE, FALSE),
  project_all = c(TRUE, FALSE)
)


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

logcounts_all_matrix <- as.matrix(rbind(logcounts(sce.rna), logcounts(sce.atac)))


# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all_matrix)))
colnames(all_celltypes) <- "celltype"

# Clusters
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

#---------------------------StabMap---------------------------------------------

stab_list <- stab_build_model(data = logcounts_all_matrix, 
                              na_features =  na_features,
                              na_cells = na_cells)


#---------------------------Result list-----------------------------------------
results <- list()

#---------------------------Loop------------------------------------------------
for (i in 1:nrow(param_grid)) {
  set.seed(42)
  params <- param_grid[i, ]
  
  # When we want both PC to be the same
  # params$ncomponentsSubset <- params$ncomponentsReference
  
  # Print current parameters 
  cat("Running stabMap with ncomponentsReference =", params$ncomponentsReference,
      ", ncomponentsSubset =", params$ncomponentsSubset,
      ", maxFeatures =", params$maxFeatures,
      ", scale.center =", params$scale.center,
      ", scale.scale =", params$scale.scale, 
      ", project_all = ", params$project_all )
  

  # Run StabMap with current parameters
time_integration <- system.time(  stab <- stabMap(
    stab_list,
    reference_list = c("all_feat"),
    ncomponentsReference = params$ncomponentsReference,
    ncomponentsSubset = params$ncomponentsSubset,
    maxFeatures = params$maxFeatures,
    scale.center = params$scale.center,
    scale.scale = params$scale.scale,
    projectAll = params$project_all,
    plot = FALSE))
  

  # UMAP
  stab_umap_coord <- as.data.frame(calculateUMAP(t(stab)))
  stab_umap <- merge( as.data.frame(metadata), stab_umap_coord, by =0 ) %>%
    mutate(isNA = if_else(Row.names %in% colnames(logcounts_all_matrix)[na_cells] , "NA", "notNA")) %>%  
    full_join(cluster) %>%
    column_to_rownames("Row.names")
  
  # Calculate silhouette scores
  stab_sil_sum <- silhouette_summary(stab_umap$n, stab_umap %>% select(UMAP1, UMAP2))
  
  
  # Calculate clustering stats
  stab_cluster_stats <- cluster.stats(dist(stab_umap %>% select(UMAP1, UMAP2)), stab_umap$n, silhouette = FALSE)
  dunn <- stab_cluster_stats$dunn
  
  # Cell Type Accuracy
  rna_train <- stab_umap$celltype[na_cells]
  names(rna_train) <- rownames(stab_umap)[na_cells]
  atac_query <- stab_umap$celltype[setdiff(1:ncol(logcounts_all_matrix), na_cells)]
  names(atac_query) <- rownames(stab_umap)[setdiff(1:ncol(logcounts_all_matrix), na_cells)]
  
  stab_knn_out = embeddingKNN(stab_umap_coord,
                              rna_train,
                              type = "uniform_fixed",
                              k_values = 5)
  stab_knn_acc = mean(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), na.rm = TRUE)
  stab_knn_acc_bal = mean(unlist(lapply(split(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
  
  # RMSE
  
  time_imputation <- system.time(imp <- imputeEmbedding(
    stab_list,
    stab,
    reference = colnames(stab_list[["all_feat"]]),
    query = colnames(stab_list[["missing_feat"]])))
  
  stab_rmse <- rmse_imp(imp_data = imp$all_feat[na_features,], 
                        real_data = logcounts_all_matrix, 
                        na_features = na_features, 
                        na_cells = na_cells)
  
  
  # Store results
  param_str <- paste(
    "ncomponentsReference", params$ncomponentsReference,
    "ncomponentsSubset", params$ncomponentsSubset,
    "maxFeatures", params$maxFeatures,
    "scale.center", params$scale.center,
    "scale.scale", params$scale.scale,
    "project_all", params$project_all,
    sep = "_")
  
  results[[param_str]] <- list(
    stab_sil_sum = stab_sil_sum,
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal,
    stab_rmse = stab_rmse,
    time_integration = time_integration[3],
    time_imputation = time_imputation[3])
  
}

#---------------------------Comparison results----------------------------------
# Create data frame
comparison <- data.frame(
  ncomponentsReference = integer(),
  ncomponentsSubset = integer(),
  maxFeatures = integer(),
  scale.center = logical(),
  scale.scale = logical(),
  mean_sil_score = numeric(),
  stab_knn_acc = numeric(),
  stab_knn_acc_bal = numeric(),
  stab_rmse = numeric(),
  time_integration = numeric(),
  time_imputation = numeric()
)

# Fill data frame
for (param_name in names(results)) {
  res <- results[[param_name]]
  mean_sil_score <- mean(res$stab_sil_sum$score, na.rm = TRUE)
  stab_knn_acc <- res$stab_knn_acc
  stab_knn_acc_bal <- res$stab_knn_acc_bal
  stab_rmse <- res$stab_rmse
  time_integration <- res$time_integration
  time_imputation <- res$time_imputation
  
  
  param_values <- unlist(strsplit(param_name, "_"))
  comparison <- rbind(comparison, data.frame(
    ncomponentsReference = as.numeric(param_values[2]),
    ncomponentsSubset = as.numeric(param_values[4]),
    maxFeatures = as.numeric(param_values[6]),
    scale.center = as.logical(param_values[8]),
    scale.scale = as.logical(param_values[10]),
    project_all = as.logical(param_values[13]),
    mean_sil_score = mean_sil_score,
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal,
    stab_rmse = stab_rmse,
    time_integration = time_integration,
    time_imputation = time_imputation
  ))
}


write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/stab_comparison1.txt")

# For this file: 
# write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/stab_comparison_numfactors1.txt")
# use
# param_grid <- expand.grid(
#   ncomponentsReference = c( 10, 15, 20, 30, 40, 50, 60, 70),
#   maxFeatures = c(1740),
#   scale.center = c(FALSE),
#   scale.scale = c(FALSE),
#   project_all = c(TRUE)
# )
# and uncomment line 84