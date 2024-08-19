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
bridge_size <- c(10, 240, 490, 740, 990, 1240, 1490, 1640, 1730)

# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename(celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14)

# Generate na_features_list
na_features_list <- list()
for (i in 1:length(bridge_size)) {
  current_bridge_size <- bridge_size[i]
  for (j in 1:5) {  # Generate 5 different sets for each bridge size
    na_features_rand <- sample(1740, current_bridge_size)
    na_features_list <- append(na_features_list, list(na_features_rand))
  }
}

# #---------------------------MOFA------------------------------------------------
# results_mofa <- data.frame(row.names = c("bridge_size", "mean_sil_score", "mofa_knn_acc", "mofa_rsme"))
# 
# # Function to run the MOFA analysis
# run_mofa_analysis <- function(bridge_size, outfile) {
#   results <- data.frame()
#   
#   for (i in 1:length(bridge_size)) {
#     current_bridge_size <- bridge_size[i]
#     
#     for (j in 1:5) {  # Loop to repeat the analysis 5 times
#       # Get the corresponding na_features from the pre-generated list
#       na_features <- na_features_list[[((i-1)*5) + j]]
#       
#       model <- mofa_build_model(data = logcounts_all,
#                                 na_features = na_features,
#                                 na_cells = na_cells,
#                                 metadata = metadata)
#       
#       trained_model <- mofa_parameter_train(num_factors = 70, spikeslab_weights = FALSE, 
#                                             model = model)
#       
#       # UMAP
#       trained_model <- run_umap(trained_model)
#       mofa_umap_coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
#       
#       # Add metadata to the UMAP results
#       mofa_umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0)  %>%
#         full_join(cluster) %>%
#         column_to_rownames("Row.names") %>%
#         select(-sample)
#       
#       # MOFA celltype
#       mofa_sil_sum <- silhouette_summary(mofa_umap$n, mofa_umap %>% select(UMAP1, UMAP2))
#       mean_sil_score <- mean(mofa_sil_sum$score)
#       print(mean_sil_score)
#       
#       # Predict "ATAC only" cell's cell types
#       rna_train <- mofa_umap$celltype[na_cells]
#       names(rna_train) <- rownames(mofa_umap)[na_cells]
#       atac_query <- mofa_umap$celltype[setdiff(1:ncol(logcounts_all), na_cells)]
#       names(atac_query) <- rownames(mofa_umap)[setdiff(1:ncol(logcounts_all), na_cells)]
#       
#       mofa_knn_out <- embeddingKNN(mofa_umap_coord, rna_train, type = "uniform_fixed", k_values = 5)
#       mofa_knn_acc_bal <- mean(unlist(lapply(split(isEqual(mofa_knn_out[names(atac_query), "predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
#       
#       print(mofa_knn_acc_bal)
#       
#       # Imputation
#       trained_model <- impute(trained_model)
#       
#       mofa_rsme <- rsme_imp(imp_data = trained_model@imputed_data$groupname_with_na[[1]][, na_cells],
#                             real_data = logcounts_all, 
#                             na_features = na_features, 
#                             na_cells = na_cells)
#       
#       print(mofa_rsme)
#       
#       # Store the results of the current iteration
#       results <- rbind(results, c(current_bridge_size, mean_sil_score, mofa_knn_acc_bal, mofa_rsme))
#     }
#   }
#   
#   colnames(results) <- c("bridge_size", "mean_sil_score", "mofa_knn_acc_bal", "mofa_rsme")
#   write.table(results, file = outfile, row.names = FALSE)
# }
# 
# # Run the analysis for num_factors = 70
# run_mofa_analysis(bridge_size, "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge_70_rand.txt")

#---------------------------StabMap---------------------------------------------
results_stab <- data.frame(row.names = c("bridge_size", "mean_sil_score", "stab_knn_acc_bal", "stab_rsme"))

run_stab_analysis <- function(bridge_size, outfile) {
  results <- data.frame()
  
  for (i in 1:length(bridge_size)) {
    current_bridge_size <- bridge_size[i]
    
    for (j in 1:5) {  # Loop to repeat the analysis 5 times
      # Get the corresponding na_features from the pre-generated list
      na_features <- na_features_list[[((i-1)*5) + j]]
      
      stab_list <- stab_build_model(data = logcounts_all, 
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
      
      # UMAP
      stab_umap_coord <- as.data.frame(calculateUMAP(t(stab)))
      
      # Add metadata to the UMAP results (celltype, if it was an NA cell, cluster)
      stab_umap <- merge(as.data.frame(metadata), stab_umap_coord, by = 0) %>%
        full_join(cluster) %>%
        column_to_rownames("Row.names")
      
      stab_sil_sum <- silhouette_summary(stab_umap$n, stab_umap %>% select(V1, V2))
      mean_sil_score <- mean(stab_sil_sum$score)
      
      rna_train <- stab_umap$celltype[na_cells]
      names(rna_train) <- rownames(stab_umap)[na_cells]
      atac_query <- stab_umap$celltype[setdiff(1:ncol(logcounts_all), na_cells)]
      names(atac_query) <- rownames(stab_umap)[setdiff(1:ncol(logcounts_all), na_cells)]
      
      # StabMap
      stab_knn_out <- embeddingKNN(stab_umap_coord, rna_train, type = "uniform_fixed", k_values = 5)
      stab_knn_acc_bal <- mean(unlist(lapply(split(isEqual(stab_knn_out[names(atac_query), "predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
      
      imp <- imputeEmbedding(stab_list, stab, reference = colnames(stab_list[["all_feat"]]), query = colnames(stab_list[["missing_feat"]]))
      
      stab_rsme <- rsme_imp(imp_data = imp$all_feat[na_features,], 
                            real_data = logcounts_all, 
                            na_features = na_features, 
                            na_cells = na_cells)
      
      # Store the results of the current iteration
      results <- rbind(results, c(current_bridge_size, mean_sil_score, stab_knn_acc_bal, stab_rsme))
    }
  }
  
  colnames(results) <- c("bridge_size", "mean_sil_score", "stab_knn_acc_bal", "stab_rsme")
  write.table(results, file = outfile, row.names = FALSE)
}

# Run the analysis for StabMap
run_stab_analysis(bridge_size, "/home/hd/hd_hd/hd_fb235/R/Data/stab_bridge_rand.txt")
