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


source("/home/hd/hd_hd/hd_fb235/R/Scripts/adaptiveKNN.R")


#---------------------------Dataset---------------------------------------------
# Peripheral Blood Mononuclear Cells provided by 10x Genomics website
# 10x Genomics Multiome technology enables simultaneous profiling of the transcriptome 
# (using 3’ gene expression) and epigenome (using ATAC-seq) from single cells to deepen 
# our understanding of how genes are expressed and regulated across different cell types.

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


#---------------------------General Data----------------------------------------

# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename(celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% group_by(celltype) %>% summarise(n = n()) %>% mutate(k = 1:14) %>% select(-n)

#---------------------------MOFA------------------------------------------------

# Put NAs in data, takes a long time
# logcounts_allNA <- logcounts_all
# logcounts_allNA[953:1740, 1:ncol(logcounts_all)/2] <- NA

logcounts_allNA <- readRDS("~/R/Data/logcountsNA.RDS")

# Create list for MOFA
mofa_list <- list(
  RNA =  logcounts_allNA[1:952,],
  ATAC = logcounts_allNA[953: 1740,]
)

# MOFA model 
model <- create_mofa(mofa_list)
samples_metadata(model) <- as.data.frame(metadata) %>% 
  rownames_to_column("sample") %>% 
  mutate(isNA =  c(rep(1, ncol(logcounts_all)/2), rep(0, ncol(logcounts_all) - ncol(logcounts_all)/2))) # 1 is TRUE (so it is NA), 0 is FALSE

MOFAobject <- prepare_mofa(model)

# Get differnt model settings
data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
# train_opts <- get_default_training_options(MOFAobject)

#---------------------------Parameters------------------------------------------
param_grid <- expand.grid(
  scale_views = c(TRUE,FALSE),
  num_factors = c(10, 15, 20, 30, 40, 50)
)

#---------------------------Result list-----------------------------------------
results <- list()

#---------------------------Loop------------------------------------------------
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  # Print current parameters 
  cat("Running MOFA with scale_views =", params$scale_views,
      ", num_factors =", params$num_factors, "\n")
  

  # Run MOFA with current parameters
  data_opts$scale_views <- params$scale_views
  model_opts$num_factors <- params$num_factors
  
  MOFAobject <- prepare_mofa(model, 
                             data_options = data_opts,
                             model_options = model_opts)
  
  trained_model <- run_mofa(MOFAobject, use_basilisk = TRUE)
  
  # UMAP
  trained_model <- run_umap(trained_model)
  mofa_umap_coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
  
  mofa_umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0)  %>%
    full_join(cluster) %>% 
    mutate(isNA =  c(rep(1, ncol(logcounts_all)/2), rep(0, ncol(logcounts_all) - ncol(logcounts_all)/2))) %>%  # 1 is TRUE, 0 is FALSE
    column_to_rownames("Row.names") %>%
    select(-sample)
  
  # Silhouette scores
  mofa_sil <- silhouette(mofa_umap$k, dist(mofa_umap %>% select(UMAP1, UMAP2)))
  mofa_sil_sum <- mofa_sil %>% 
    as.data.frame() %>% group_by(cluster) %>% summarise(score = mean(sil_width), 
                                                        frac_pos = sum(sil_width > 0)/n(),
                                                        pos_score = sum((sil_width>0)*sil_width)/sum(sil_width > 0))
  
  # Clustering stats
  mofa_cluster_stats <- cluster.stats(dist(mofa_umap %>% select(UMAP1, UMAP2)), mofa_umap$k)
  
  # Cell Type Accuracy
  rna_train <- mofa_umap$celltype[1:5016]
  names(rna_train) <- rownames(mofa_umap)[1:5016]
  atac_query <- mofa_umap$celltype[5016:10032]
  names(atac_query) <- rownames(mofa_umap)[5016:10032]
  
  mofa_knn_out = embeddingKNN(mofa_umap_coord,
                               rna_train,
                               type = "uniform_fixed",
                               k_values = 5)
  mofa_knn_acc = mean(isEqual(mofa_knn_out[names(atac_query),"predicted_labels"], atac_query), na.rm = TRUE)
  mofa_knn_acc_bal = mean(unlist(lapply(split(isEqual(mofa_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
  
  # RSME
  trained_model <- impute(trained_model)
  
  mofa_imp_comp <- as.data.frame(trained_model@imputed_data$ATAC[[1]][,1:(ncol(logcounts_all)/2)]) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
    full_join(as.data.frame(as.matrix(logcounts_all[953:1740, 1:(ncol(logcounts_all)/2)])) %>%
                rownames_to_column("feature") %>%
                pivot_longer(-feature, names_to = "sample", values_to = "actual"))
  
  mofa_rsme <- sqrt(mean((mofa_imp_comp$actual - mofa_imp_comp$predicted)^2))
  
  # Store results
  param_str <- paste(
    "scale_views", params$scale_views,
    "num_factor", params$num_factor,
    sep = "_"
  )
  
  results[[param_str]] <- list(
    mofa_umap = mofa_umap,
    mofa_sil_sum = mofa_sil_sum,
    mofa_cluster_stats = mofa_cluster_stats,
    mofa_knn_acc = mofa_knn_acc,
    mofa_knn_acc_bal = mofa_knn_acc_bal,
    mofa_rsme = mofa_rsme
  )
  
  
}

#---------------------------Comparison results----------------------------------
# Create data frame
comparison <- data.frame(
  scale_views = logical(),
  num_factor = integer(),
  mean_sil_score = numeric(),
  min_sil_score = numeric(),
  min_sil_score_name = character(),
  max_sil_score = numeric(),
  mofa_knn_acc = numeric(),
  mofa_knn_acc_bal = numeric(),
  mofa_rsme = numeric()
)

# Fill data frame
for (param_name in names(results)) {
  res <- results[[param_name]]
  mean_sil_score <- mean(res$mofa_sil_sum$score, na.rm = TRUE)
  min_sil_score <- min(res$mofa_sil_sum$score,na.rm = TRUE)
  min_sil_score_name <- cluster$celltype[which.min(res$mofa_sil_sum$score)]
  max_sil_score <- max(res$mofa_sil_sum$score,na.rm = TRUE)
  mofa_knn_acc <- res$mofa_knn_acc
  mofa_knn_acc_bal <- res$mofa_knn_acc_bal
  mofa_rsme <- res$mofa_rsme
  
  param_values <- unlist(strsplit(param_name, "_"))
  comparison <- rbind(comparison, data.frame(
    scale_views = as.logical(param_values[3]),
    num_factor = as.numeric(param_values[6]),
    mean_sil_score = mean_sil_score,
    min_sil_score = min_sil_score,
    min_sil_score_name = min_sil_score_name,
    max_sil_score = max_sil_score,
    mofa_knn_acc = mofa_knn_acc,
    mofa_knn_acc_bal = mofa_knn_acc_bal,
    mofa_rsme = mofa_rsme
  ))
}

print(comparison)

write.table(comparison, file = "/home/hd/hd_hd/hd_fb235/R/Data/mofa_comparison.txt")




