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
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename( celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

na_features_1 <- 953:1740
na_cells_1 <- 1:5016
na_cells_2 <- 5017:10032

#---------------------------MOFA------------------------------------------------
mofa_build_model_2 <- function(data, na_features_1, na_cells_1, na_features_2, na_cells_2, groupname_with_na, groupname_without_na, metadata) {
  
  logcounts_allNA <- data
  logcounts_allNA[na_features_1, na_cells_1] <- NA
  logcounts_allNA[na_features_2, na_cells_2] <- NA
  
  mofa_list <- list(
    query = logcounts_allNA[na_features_1, ],
    reference = logcounts_allNA[setdiff(1:nrow(logcounts_allNA), na_features_1), ]
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

stab_build_model_2 <- function(data, na_features_1, na_cells_1, na_features_2, na_cells_2) {
  
  
  stab_list <- list(
    query = logcounts_all[setdiff(1:nrow(data), na_features_1), na_cells_1],
    reference = logcounts_all[setdiff(1:nrow(data), na_features_2),na_cells_2]
  )
  
  return(stab_list)
}


run_bridge_analysis <- function(method, bridge_size, outfile, data, metadata, na_cells_1, na_features_1, na_cells_2, all_celltypes, cluster) {
  results <- data.frame()
  
  for (i in 1:length(bridge_size)) {
    current_bridge_size <- bridge_size[i]
    print(paste("Current bridge size:", current_bridge_size))
    
    na_features_2 <-  1:current_bridge_size
    
    if (method == "mofa") {
      model <- mofa_build_model_2(data = data,
                                na_features_1 = na_features_1,
                                na_cells_1 = na_cells_1,
                                na_features_2 = na_features_2,
                                na_cells_2 = na_cells_2,
                                metadata = metadata)
      
      trained_model <- mofa_parameter_train(num_factors = 70, spikeslab_weights = FALSE,
                                            model = model)
      
      trained_model <- run_umap(trained_model)
      
      coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
      umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0) %>%
        full_join(cluster) %>% 
        column_to_rownames("Row.names") %>%
        select(-sample)
      
      trained_model <- impute(trained_model)
      imp_data_1 <- trained_model@imputed_data$query[[1]][,na_cells_1]
      imp_data_2 <- trained_model@imputed_data$reference[[1]][,na_cells_2]

      
    } else if (method == "stab") {
      stab_list <- stab_build_model_2(data = data,
                                    na_features_1 = na_features_1,
                                    na_cells_1 = na_cells_1,
                                    na_features_2 = na_features_2,
                                    na_cells_2 = na_cells_2)
      
      stab <- stabMap(stab_list,
                      ncomponentsReference = 70,
                      ncomponentsSubset = 70,
                      reference_list = c("reference"),
                      plot = FALSE,
                      scale.center = FALSE,
                      scale.scale = FALSE,
                      maxFeatures = 1740)
      
      coord <- as.data.frame(calculateUMAP(t(stab)))
      umap <- merge( coord, as.data.frame(metadata), by = 0) %>%
        full_join(cluster) %>%
        column_to_rownames("Row.names")
      
      imp <- imputeEmbedding(stab_list, stab)
      imp_data_1 <- imp$query[-na_features_1,]
      imp_data_2 <- imp$reference[-na_features_2,]
      
    } else {
      stop("Invalid method. Choose 'mofa' or 'stab'.")
    }
    
    # Common analysis steps
    sil_sum <- silhouette_summary(umap$n, umap %>% select(1:2))
    mean_sil_score <- mean(sil_sum$score)
    
    rna_train <- umap$celltype[na_cells_2]
    names(rna_train) <- rownames(umap)[na_cells_2]
    atac_query <- umap$celltype[setdiff(1:ncol(data), na_cells_2)]
    names(atac_query) <- rownames(umap)[setdiff(1:ncol(data), na_cells_2)]
    
    knn_out <- embeddingKNN(coord, rna_train, type = "uniform_fixed", k_values = 5)
    knn_acc_bal <- mean(unlist(lapply(split(isEqual(knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
    
    rmse_1 <- rmse_imp(imp_data = imp_data_1,
                     real_data = data,
                     na_features = -na_features_1,
                     na_cells = na_cells_1)
    
    rmse_2 <- rmse_imp(imp_data = imp_data_2,
                     real_data = data,
                     na_features = -na_features_2,
                     na_cells = na_cells_2)
    
    results <- rbind(results, c(current_bridge_size, mean_sil_score, knn_acc_bal, rmse_1, rmse_2))
  }
  
  colnames(results) <- c("bridge_size", "mean_sil_score", "knn_acc_bal", "rmse_1", "rmse_2")
  write.table(results, file = outfile, row.names = FALSE)
}
#---------------------------MOFA functions--------------------------------------
run_bridge_analysis("mofa", c(900,17 ), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_sep_3.txt",
                    metadata = metadata, cluster = cluster, na_cells_1 = na_cells_1, na_cells_2 = na_cells_2, na_features_1 = na_features_1, all_celltypes = all_celltypes, data = logcounts_all)

run_bridge_analysis("stab", c(900,17 ), "/home/hd/hd_hd/hd_fb235/R/Data/stab_bridge_sep_3.txt",
                    metadata = metadata, cluster = cluster, na_cells_1 = na_cells_1, na_cells_2 = na_cells_2, na_features_1 = na_features_1, all_celltypes = all_celltypes, data = logcounts_all)

#----
# Function to run the MOFA analysis
# run_mofa_analysis <- function(bridge_size, outfile) {
#   results <- data.frame()
#   
#   for (i in 1:length(bridge_size)) {
#     current_bridge_size <- bridge_size[i] 
#     print(current_bridge_size)
#     na_features_2 <-  1:current_bridge_size
# 
#     
#     model <- mofa_build_model(data = logcounts_all,
#                               na_features_1 = na_features_1,
#                               na_cells_1 = na_cells_1,
#                               na_features_2 = na_features_2,
#                               na_cells_2 = na_cells_2,
#                               metadata = metadata)
#     
#     trained_model <- mofa_parameter_train(num_factors = 70, spikeslab_weights = FALSE, 
#                                           model = model)
#     
#     # UMAP
#     trained_model <- run_umap(trained_model)
#     mofa_umap_coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
#     
#     # Add metadata to the UMAP results
#     mofa_umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0)  %>%
#       full_join(cluster) %>% 
#       column_to_rownames("Row.names") %>%
#       select(-sample)
#     
#     # MOFA celltype
#     mofa_sil_sum <- silhouette_summary(mofa_umap$n, mofa_umap %>% select(UMAP1, UMAP2))
#     mean_sil_score <- mean(mofa_sil_sum$score)
#     print(mean_sil_score)
#     
#     # Predict "ATAC only" cell's cell types
#     rna_train <- mofa_umap$celltype[na_cells_2]
#     names(rna_train) <- rownames(mofa_umap)[na_cells_2]
#     atac_query <- mofa_umap$celltype[setdiff(1:ncol(logcounts_all), na_cells_2)]
#     names(atac_query) <- rownames(mofa_umap)[setdiff(1:ncol(logcounts_all), na_cells_2)]
#     
#     mofa_knn_out <- embeddingKNN(mofa_umap_coord, rna_train, type = "uniform_fixed", k_values = 5)
#     mofa_knn_acc_bal = mean(unlist(lapply(split(isEqual(mofa_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
#     
#     print(mofa_knn_acc_bal)
#     
#     # Imputation
#     trained_model <- impute(trained_model)
#     
#     mofa_rmse_1 <- rmse_imp(imp_data = trained_model@imputed_data$query[[1]][,na_cells_1],
#                          real_data = logcounts_all,
#                          na_features = na_features_1,
#                          na_cells = na_cells_1)
#     
#     mofa_rmse_2 <- rmse_imp(imp_data = trained_model@imputed_data$reference[[1]][,na_cells_2],
#                             real_data = logcounts_all,
#                             na_features = na_features_2,
#                             na_cells = na_cells_2)
#     
#     
#     
#     results <- rbind(results, c(current_bridge_size, mean_sil_score, mofa_knn_acc_bal, mofa_rmse_1, mofa_rmse_2))
#   }
#   
#   colnames(results) <- c("bridge_size", "mean_sil_score", "mofa_knn_acc_bal", "mofa_rmse_1", "mofa_rmse_2")
#   write.table(results, file = outfile, row.names = FALSE)
# }

# Run the analysis for num_factors = 70
run_mofa_analysis(c(10 , 240 ), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_70_RNAout.txt")


#---------------------------StabMap---------------------------------------------



# run_stab_analysis <- function(bridge_size, outfile) {
#   results <- data.frame()
#   for (i in 1:length(bridge_size)) {
# 
#     current_bridge_size <- bridge_size[i]
#     print(current_bridge_size)
#     na_features_2 <-  1:current_bridge_size
# 
#     stab_list <- stab_build_model(data = logcounts_all,
#                                   na_features_1 = na_features_1,
#                                   na_cells_1 = na_cells_1,
#                                   na_features_2 = na_features_2,
#                                   na_cells_2 = na_cells_2)
# 
# 
#     stab <- stabMap(stab_list,
#                     ncomponentsReference = 70,
#                     ncomponentsSubset = 70,
#                     reference_list = c("reference"),
#                     plot = FALSE,
#                     scale.center = FALSE,
#                     scale.scale = FALSE,
#                     maxFeatures = 1740)
# 
#     # UMAP
#     stab_umap_coord <- as.data.frame(calculateUMAP(t(stab)))
# 
# 
#     # Add metadata to the UMAP results (celltype, if it was an NA cell, cluster)
#     stab_umap <- merge( as.data.frame(metadata), stab_umap_coord, by =0 ) %>%
#       full_join(cluster) %>%
#       column_to_rownames("Row.names")
# 
#     stab_sil_sum <- silhouette_summary(stab_umap$n, stab_umap %>% select(V1, V2))
#     mean_sil_score <- mean(stab_sil_sum$score)
# 
#     rna_train <- stab_umap$celltype[na_cells]
#     names(rna_train) <- rownames(stab_umap)[na_cells]
#     atac_query <- stab_umap$celltype[setdiff(1:ncol(logcounts_all), na_cells)]
#     names(atac_query) <- rownames(stab_umap)[setdiff(1:ncol(logcounts_all), na_cells)]
# 
#     # StabMap
#     stab_knn_out = embeddingKNN(stab_umap_coord,
#                                 rna_train,
#                                 type = "uniform_fixed",
#                                 k_values = 5)
#     stab_knn_acc_bal = mean(unlist(lapply(split(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
# 
#     imp = imputeEmbedding(
#       stab_list,
#       stab)
# 
#     stab_rmse_1 <- rmse_imp(imp_data = imp$query[-na_features_1,],
#                           real_data = logcounts_all,
#                           na_features = -na_features_1,
#                           na_cells = na_cells_1)
# 
#     stab_rmse_2 <- rmse_imp(imp_data = imp$reference[-na_features_2,],
#                           real_data = logcounts_all,
#                           na_features = -na_features_2,
#                           na_cells = na_cells_2)
# 
# 
# 
# 
# 
# 
# 
#     results <- rbind(results, c(current_bridge_size, mean_sil_score, stab_knn_acc_bal, stab_rmse_1, stab_rmse_2))
#   }
# 
#   colnames(results) <- c("bridge_size", "mean_sil_score", "stab_knn_acc_bal", "stab_rmse_1", "stab_rmse_2")
#   write.table(results, file = outfile, row.names = FALSE)
# }


# run_stab_analysis(c(10, 240 ), "/home/hd/hd_hd/hd_fb235/R/Data/stab_bridge_RNAout.txt")


