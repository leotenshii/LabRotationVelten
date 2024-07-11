#---------------------------Seed------------------------------------------------
set.seed(42)

#---------------------------Prep for running script-----------------------------
#---------------------------Libraries-------------------------------------------
suppressMessages(c(library(scater),
library(scran),
library(StabMap),
library(tidyverse),
library(cluster),
library(patchwork),
library(fpc),
library(MultiAssayExperiment)))


source("/home/hd/hd_hd/hd_fb235/R/Scripts/adaptiveKNN.R")

## Preparing Data




#---------------------------Dataset---------------------------------------------
# Peripheral Blood Mononuclear Cells provided by 10x Genomics website
# 10x Genomics Multiome technology enables simultaneous profiling of the transcriptome 
# (using 3’ gene expression) and epigenome (using ATAC-seq) from single cells to deepen 
# our understanding of how genes are expressed and regulated across different cell types.

# mae <- scMultiome("pbmc_10x", mode = "*", dry.run = FALSE, format = "MTX")

# Loaded from RDS cause I cant install the SingleCellMultiModal package
mae <- readRDS("/home/hd/hd_hd/hd_fb235/R/Data/data.RDS")
metadata <- mae@colData




#---------------------------Data preparation------------------------------------
# RNA
sce.rna <- normalize_and_select_features(experiments(mae)[["rna"]], 0.01, 0.05)


# ATAC
sce.atac <- normalize_and_select_features(experiments(mae)[["atac"]], 0.01, 0.05)

logcounts_all <- rbind(logcounts(sce.rna), logcounts(sce.atac))

## Some general dataframes etc.


# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename(celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters
cluster <- as.data.frame(metadata) %>% group_by(celltype) %>% summarise(n = n()) %>% mutate(k = 1:14) %>% select(-n)


#---------------------------StabMap---------------------------------------------

# Seperation ATAC Multiome
names <- c(rep("ATAC", ncol(logcounts_all)/2), rep("Multiome", ncol(logcounts_all)/2))

# Assay Types
assayType = ifelse(rownames(logcounts_all) %in% rownames(sce.rna),
                   "rna", "atac")

# List for StabMap
assay_list = list(
  ATAC = logcounts_all[assayType %in% c("rna"), names %in% c("ATAC")],
  Multiome = logcounts_all[assayType %in% c("rna", "atac"), names %in% c("Multiome")]
)



#---------------------------Parameters------------------------------------------
param_grid <- expand.grid(
  ncomponentsReference = c( 10, 15, 20, 30, 40, 50, 60, 70), # More ncomponents can capture more variance but may include noise. -> Should be similar?
  # ncomponentsSubset = c( 15, 30, 40, 50, 70),
  maxFeatures = c(1000), # Too many maybe overfitting
  scale.center = c( FALSE), # False when mean or variance carry some important information?
  scale.scale = c( FALSE),
  project_all = c(FALSE) # might help in refining the alignment
)


# restrictFeatures should be false for single hop (-> documentation)


#---------------------------Result list-----------------------------------------
results <- list()

#---------------------------Loop------------------------------------------------
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  params$ncomponentsSubset <- params$ncomponentsReference
  
  # Print current parameters 
  cat("Running stabMap with ncomponentsReference =", params$ncomponentsReference,
      ", ncomponentsSubset =", params$ncomponentsSubset,
      ", maxFeatures =", params$maxFeatures,
      ", scale.center =", params$scale.center,
      ", scale.scale =", params$scale.scale, 
      ", project_all = ", params$project_all )
  

  # Run StabMap with current parameters
time_integration <- system.time(  stab <- stabMap(
    assay_list,
    reference_list = c("Multiome"),
    ncomponentsReference = params$ncomponentsReference,
    ncomponentsSubset = params$ncomponentsSubset,
    maxFeatures = params$maxFeatures,
    scale.center = params$scale.center,
    scale.scale = params$scale.scale,
    projectAll = params$project_all,
    plot = FALSE
  ))
  

  # UMAP
  stab_umap_coord <- as.data.frame(calculateUMAP(t(stab)))
  stab_umap <- merge(as.data.frame(metadata), stab_umap_coord, by = 0) %>%
    mutate(isNA = ifelse(Row.names %in% colnames(assay_list$ATAC), 1, 0)) %>%
    full_join(cluster) %>%
    column_to_rownames("Row.names")
  
  # Calculate silhouette scores
  stab_sil_sum <- silhoutte_summary(stab_umap$k, stab_umap %>% select(V1, V2))
  
  
  # Calculate clustering stats
  stab_cluster_stats <- cluster.stats(dist(stab_umap %>% select(V1, V2)), stab_umap$k, silhouette = FALSE)
  dunn <- stab_cluster_stats$dunn
  
  # Cell Type Accuracy
  rna_train <- stab_umap$celltype[1:5016]
  names(rna_train) <- rownames(stab_umap)[1:5016]
  atac_query <- stab_umap$celltype[5016:10032]
  names(atac_query) <- rownames(stab_umap)[5016:10032]
  
  stab_knn_out = embeddingKNN(stab_umap_coord,
                              rna_train,
                              type = "uniform_fixed",
                              k_values = 5)
  stab_knn_acc = mean(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), na.rm = TRUE)
  stab_knn_acc_bal = mean(unlist(lapply(split(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
  
  # RSME
  
  time_imputation <- system.time(imp <- imputeEmbedding(
    assay_list,
    stab,
    reference = colnames(assay_list[["Multiome"]]),
    query = colnames(assay_list[["ATAC"]])))
  
  stab_imp_comp <- as.data.frame(as.matrix(imp$Multiome)[953:1740,1:(ncol(logcounts_all)/2)]) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
    full_join(as.data.frame(as.matrix(logcounts_all[953:1740, 1:(ncol(logcounts_all)/2)])) %>%
                rownames_to_column("feature") %>%
                pivot_longer(-feature, names_to = "sample", values_to = "actual"))
  
  stab_rsme <- sqrt(mean((stab_imp_comp$actual - stab_imp_comp$predicted)^2))
  
  # Comparison Baseline
  stab_baseline_comp <- as.data.frame(as.matrix(imp$Multiome)[953:1740,]) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
    full_join(as.data.frame(rowMeans(as.data.frame(as.matrix(logcounts_all[assayType %in% c("atac"), 5017:10032])))) %>%
                rownames_to_column("feature")) %>%
    rename(baseline = `rowMeans(as.data.frame(as.matrix(logcounts_all[assayType %in% c("atac"), 5017:10032])))`)
  
  stab_mae <- mean(abs(stab_baseline_comp$predicted - stab_baseline_comp$baseline))
  
  # Store results
  param_str <- paste(
    "ncomponentsReference", params$ncomponentsReference,
    "ncomponentsSubset", params$ncomponentsSubset,
    "maxFeatures", params$maxFeatures,
    "scale.center", params$scale.center,
    "scale.scale", params$scale.scale,
    "project_all", params$project_all,
    sep = "_"
  )
  
  results[[param_str]] <- list(
    #stab_umap = stab_umap,
    stab_sil_sum = stab_sil_sum,
    dunn = dunn,
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal,
    stab_rsme = stab_rsme,
    stab_mae = stab_mae,
    time_integration = time_integration[3],
    time_imputation = time_imputation[3]
  )
  
  
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
  min_sil_score = numeric(),
  max_sil_score = numeric(),
  dunn = numeric(),
  stab_knn_acc = numeric(),
  stab_knn_acc_bal = numeric(),
  stab_rsme = numeric(),
  stab_mae = numeric(),
  time_integration = numeric(),
  time_imputation = numeric()
)

# Fill data frame
for (param_name in names(results)) {
  res <- results[[param_name]]
  mean_sil_score <- mean(res$stab_sil_sum$score, na.rm = TRUE)
  min_sil_score <- min(res$stab_sil_sum$score,na.rm = TRUE)
  max_sil_score <- max(res$stab_sil_sum$score,na.rm = TRUE)
  stab_knn_acc <- res$stab_knn_acc
  stab_knn_acc_bal <- res$stab_knn_acc_bal
  stab_rsme <- res$stab_rsme
  stab_mae <- res$stab_mae
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
    min_sil_score = min_sil_score,
    max_sil_score = max_sil_score,
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal,
    stab_rsme = stab_rsme,
    stab_mae = stab_mae,
    time_integration = time_integration,
    time_imputation = time_imputation
  ))
}

print(comparison)

write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/stab_comparison_numfactors.txt")
