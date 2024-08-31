# Last changes on 28.08.2024
# Author: Leoni Zimmermann

#---------------------------Description-----------------------------------------
# This script runs MOFA and StabMap with different bridge sizes (number of overlapping ATAC features).

#---------------------------Seed------------------------------------------------
set.seed(42)
#---------------------------Libraries-------------------------------------------
suppressMessages(c(library(scater),
                   library(scran),
                   library(MOFA2),
                   library(tidyverse),
                   library(cluster),
                   library(patchwork),
                   library(fpc),
                   library(MultiAssayExperiment),
                   library(StabMap)))


source("~/R/Functions/data_prep_functions.R")
source("~/R/Functions/integration_metrics_functions.R")
source("~/R/Functions/adaptiveKNN.R")


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
na_cells <- 1:5016

# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all)))
colnames(all_celltypes) <- "celltype"

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

# Random NA features
na_features_list <- list()
bridge_size <- c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730)
for (i in 1:length(bridge_size)) {
  current_bridge_size <- bridge_size[i]
  for (j in 1:5) {  
    na_features_rand <- sample(1740, (1740-current_bridge_size))
    na_features_list <- append(na_features_list, list(na_features_rand))
  }
}



run_bridge_analysis <- function(method, bridge_size, outfile, data, metadata, na_cells, all_celltypes, cluster, random_na_features = FALSE, random_na_features_list = NULL, seperate_rmse = FALSE, num_factor = NULL) {
  results <- data.frame()
  
  for (i in 1:length(bridge_size)) {
    current_bridge_size <- bridge_size[i]
    print(paste("Current bridge size:",  current_bridge_size))
    
    if (random_na_features) {
      na_features <- random_na_features_list[[i]]
    } else {
      na_features <- current_bridge_size:1740
    }
    
    
    if (method == "mofa") {
      model <- mofa_build_model(data = data,
                                na_features = na_features,
                                na_cells = na_cells,
                                metadata = metadata)
      
      trained_model <- mofa_parameter_train(num_factors = num_factor, spikeslab_weights = FALSE,
                                            model = model)
      
      trained_model <- run_umap(trained_model)
      
      coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
      umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0) %>%
        full_join(cluster) %>% 
        column_to_rownames("Row.names") %>%
        select(-sample)
      
      trained_model <- impute(trained_model)
      imp_data <- trained_model@imputed_data$missing_feat[[1]][,na_cells]
      
    } else if (method == "stab") {
      stab_list <- stab_build_model(data = data, 
                                    na_features = na_features,
                                    na_cells = na_cells)
      
      stab <- stabMap(stab_list,
                      ncomponentsReference = 70,
                      ncomponentsSubset = 70,
                      reference_list = c("all_feat"),
                      plot = FALSE,
                      scale.center = FALSE,
                      scale.scale = FALSE,
                      maxFeatures = 1740)
      
      coord <- as.data.frame(calculateUMAP(t(stab)))
      umap <- merge( coord, as.data.frame(metadata), by = 0) %>%
        full_join(cluster) %>%
        column_to_rownames("Row.names")
      
      imp <- imputeEmbedding(stab_list, stab, 
                             reference = colnames(stab_list[["all_feat"]]), 
                             query = colnames(stab_list[["missing_feat"]]))
      imp_data <- imp$all_feat[na_features,]
      
    } else {
      stop("Invalid method. Choose 'mofa' or 'stab'.")
    }
    
    # Common analysis steps
    sil_sum <- silhouette_summary(umap$n, umap %>% select(1:2))
    mean_sil_score <- mean(sil_sum$score)
    
    rna_train <- umap$celltype[na_cells]
    names(rna_train) <- rownames(umap)[na_cells]
    atac_query <- umap$celltype[setdiff(1:ncol(data), na_cells)]
    names(atac_query) <- rownames(umap)[setdiff(1:ncol(data), na_cells)]
    
    knn_out <- embeddingKNN(coord, rna_train, type = "uniform_fixed", k_values = 5)
    knn_acc_bal <- mean(unlist(lapply(split(isEqual(knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
    
    rmse <- rmse_imp(imp_data = imp_data,
                     real_data = data,
                     na_features = na_features,
                     na_cells = na_cells)
    
    rmse_atac <- NA
    rmse_rna <- NA
    
    if (seperate_rmse) {
      rmse_tbl <- rmse_imp_table(imp_data = imp_data,
                                 real_data = data,
                                 na_features = na_features,
                                 na_cells = na_cells)
      
      rmse_tbl_atac <- rmse_tbl %>%
        filter(feature %in% rownames(sce.atac))
      
      rmse_tbl_rna <- rmse_tbl %>%
        filter(feature %in% rownames(sce.rna))
      
      rmse_atac <- sqrt(mean((rmse_tbl_atac$actual - rmse_tbl_atac$predicted)^2, na.rm = TRUE))
      rmse_rna <- sqrt(mean((rmse_tbl_rna$actual - rmse_tbl_rna$predicted)^2, na.rm = TRUE))
    }
    
    results <- rbind(results, c(current_bridge_size, mean_sil_score, knn_acc_bal, rmse, rmse_atac, rmse_rna))
    
  }
  
  colnames(results) <- c("bridge_size", "mean_sil_score", "knn_acc_bal", "rmse", "rmse_atac", "rmse_rna")
  write.table(results, file = outfile, row.names = FALSE)
}

# run_bridge_analysis("stab", c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730), "/home/hd/hd_hd/hd_fb235/R/Data/stab_bridge.txt",
#                     metadata = metadata, cluster = cluster, na_cells = na_cells, all_celltypes = all_celltypes, data = logcounts_all, seperate_rmse = TRUE)
# 
# run_bridge_analysis("mofa", num_factor = 70, c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_70.txt",
#                     metadata = metadata, cluster = cluster, na_cells = na_cells, all_celltypes = all_celltypes, data = logcounts_all, seperate_rmse = TRUE)
# 
# run_bridge_analysis("mofa", num_factor = 15, c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_15.txt",
#                     metadata = metadata, cluster = cluster, na_cells = na_cells, all_celltypes = all_celltypes, data = logcounts_all)
# 
# 
# run_bridge_analysis("stab", rep(c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730), each = 5), "/home/hd/hd_hd/hd_fb235/R/Data/stab_bridge_rand.txt",
#                     metadata = metadata, cluster = cluster, na_cells = na_cells, all_celltypes = all_celltypes, data = logcounts_all, random_na_features = TRUE, random_na_features_list = na_features_list)
# 
# run_bridge_analysis("mofa", num_factor = 70, rep(c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730), each = 5), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_70_rand.txt",
#                     metadata = metadata, cluster = cluster, na_cells = na_cells, all_celltypes = all_celltypes, data = logcounts_all, random_na_features = TRUE, random_na_features_list = na_features_list)

run_bridge_analysis("mofa", num_factor = 15, rep(c(10, 52, 104, 156, 218, 438, 653, 870 ,953 ,1088, 1305, 1523, 1730), each = 5), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_15_rand.txt",
                    metadata = metadata, cluster = cluster, na_cells = na_cells, all_celltypes = all_celltypes, data = logcounts_all, random_na_features = TRUE, random_na_features_list = na_features_list)




