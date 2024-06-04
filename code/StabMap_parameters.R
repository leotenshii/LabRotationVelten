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




#---------------------------Data preperation------------------------------------
# Normalization RNA
sce.rna <- experiments(mae)[["rna"]]
sce.rna <- logNormCounts(sce.rna)

# Feature selection
decomp <- modelGeneVar(sce.rna)
hvgs <- rownames(decomp)[decomp$mean>0.01 & decomp$p.value <= 0.05]

sce.rna <- sce.rna[hvgs,]


# Normalization ATAC
sce.atac <- experiments(mae)[["atac"]]
sce.atac <- logNormCounts(sce.atac)

# Feature selection using highly variable peaks
decomp <- modelGeneVar(sce.atac)
hvgs <- rownames(decomp)[decomp$mean>0.25
                         & decomp$p.value <= 0.05]

sce.atac <- sce.atac[hvgs,]

logcounts_all <- rbind(logcounts(sce.rna), logcounts(sce.atac))


## Some general dataframes etc.


# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename(celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters
cluster <- as.data.frame(metadata) %>% group_by(celltype) %>% summarise(n = n()) %>% mutate(k = 1:14) %>% select(-n)



## StabMap


#---------------------------StabMap---------------------------------------------

# Seperation ATAC Multiome
names <- c(rep("ATAC", ncol(logcounts_all)/2), rep("Multiome", ncol(logcounts_all)/2))

# Assay Types
assayType = ifelse(rownames(logcounts_all) %in% rownames(sce.rna),
                   "rna", "atac")

# List for StabMap
assay_list = list(
  ATAC = logcounts_all[assayType %in% c("atac"), names %in% c("ATAC")],
  Multiome = logcounts_all[assayType %in% c("rna", "atac"), names %in% c("Multiome")]
)



#---------------------------Parameters------------------------------------------
param_grid <- expand.grid(
  ncomponentsReference = c( 30, 50, 70),
  ncomponentsSubset = c( 30,  50, 70),
  maxFeatures = c(900, 1000, 1100),
  scale.center = c(TRUE, FALSE),
  scale.scale = c(TRUE, FALSE)
)

#---------------------------Result list-----------------------------------------
results <- list()

#---------------------------Loop------------------------------------------------
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Print current parameters 
  cat("Running stabMap with ncomponentsReference =", params$ncomponentsReference,
      ", ncomponentsSubset =", params$ncomponentsSubset,
      ", maxFeatures =", params$maxFeatures,
      ", scale.center =", params$scale.center,
      ", scale.scale =", params$scale.scale, "\n")
  
  # Run StabMap with current parameters
  stab <- stabMap(
    assay_list,
    reference_list = c("Multiome"),
    ncomponentsReference = params$ncomponentsReference,
    ncomponentsSubset = params$ncomponentsSubset,
    maxFeatures = params$maxFeatures,
    scale.center = params$scale.center,
    scale.scale = params$scale.scale,
    plot = FALSE
  )
  
  # UMAP
  stab_umap_coord <- as.data.frame(calculateUMAP(t(stab)))
  stab_umap <- merge(as.data.frame(metadata), stab_umap_coord, by = 0) %>%
    mutate(isNA = ifelse(Row.names %in% colnames(assay_list$ATAC), 1, 0)) %>%
    full_join(cluster) %>%
    column_to_rownames("Row.names")
  
  # Calculate silhouette scores
  stab_sil <- silhouette(stab_umap$k, dist(stab_umap %>% select(V1, V2)))
  stab_sil_sum <- stab_sil %>%
    as.data.frame() %>%
    group_by(cluster) %>%
    summarise(score = mean(sil_width), frac_pos = sum(sil_width > 0) / n(), pos_score = sum((sil_width > 0) * sil_width) / sum(sil_width > 0))
  
  # Calculate clustering stats
  stab_cluster_stats <- cluster.stats(dist(stab_umap %>% select(V1, V2)), stab_umap$k)
  
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
  
  imp = imputeEmbedding(
    assay_list,
    stab,
    reference = colnames(assay_list[["Multiome"]]),
    query = colnames(assay_list[["ATAC"]]))
  
  stab_imp_comp <- as.data.frame(as.matrix(imp$Multiome)[953:1740,]) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
    full_join(as.data.frame(as.matrix(logcounts_all[953:1740, 1:(ncol(logcounts_all)/2)])) %>%
                rownames_to_column("feature") %>%
                pivot_longer(-feature, names_to = "sample", values_to = "actual"))
  
  stab_rsme <- sqrt(mean((stab_imp_comp$actual - stab_imp_comp$predicted)^2))
  
  # Store results
  param_str <- paste(
    "ncomponentsReference", params$ncomponentsReference,
    "ncomponentsSubset", params$ncomponentsSubset,
    "maxFeatures", params$maxFeatures,
    "scale.center", params$scale.center,
    "scale.scale", params$scale.scale,
    sep = "_"
  )
  
  results[[param_str]] <- list(
    stab_umap = stab_umap,
    stab_sil_sum = stab_sil_sum,
    stab_cluster_stats = stab_cluster_stats,
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal,
    stab_rsme = stab_rsme
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
  min_sil_score_name = character(),
  max_sil_score = numeric(),
  stab_knn_acc = numeric(),
  stab_knn_acc_bal = numeric(),
  stab_rsme = numeric()
)

# Fill data frame
for (param_name in names(results)) {
  res <- results[[param_name]]
  mean_sil_score <- mean(res$stab_sil_sum$score, na.rm = TRUE)
  min_sil_score <- min(res$stab_sil_sum$score,na.rm = TRUE)
  min_sil_score_name <- cluster$celltype[which.min(res$stab_sil_sum$score)]
  max_sil_score <- max(res$stab_sil_sum$score,na.rm = TRUE)
  stab_knn_acc <- res$stab_knn_acc
  stab_knn_acc_bal <- res$stab_knn_acc_bal
  stab_rsme <- res$stab_rsme
  
  param_values <- unlist(strsplit(param_name, "_"))
  comparison <- rbind(comparison, data.frame(
    ncomponentsReference = as.numeric(param_values[2]),
    ncomponentsSubset = as.numeric(param_values[4]),
    maxFeatures = as.numeric(param_values[6]),
    scale.center = as.logical(param_values[8]),
    scale.scale = as.logical(param_values[10]),
    mean_sil_score = mean_sil_score,
    min_sil_score = min_sil_score,
    min_sil_score_name = min_sil_score_name,
    max_sil_score = max_sil_score,
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal,
    stab_rsme = stab_rsme
  ))
}

print(comparison)

write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/comparison.txt")
