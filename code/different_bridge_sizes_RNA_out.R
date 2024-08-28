# Last changes on 28.08.2024
# Author: Leoni Zimmermann

#---------------------------Description-----------------------------------------
# In this script, not only have the ATAC features of the query data been removed, but also some RNA features have been excluded from the reference dataset.
# The aim of this script is to investigate how the integration and imputation of data is affected by different bridge sizes (i.e. the number of overlapping RNA features between the reference and the query data).

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
# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all)))
colnames(all_celltypes) <- "celltype"

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

na_cells <- 1:5016
na_features_ATAC <- 953:1740
na_cells_ATAC <- 1:5016
na_cells_RNA <- 5017:10032

#---------------------------Model functions-------------------------------------
# I need to rewrite thme to be able to handle to regions of values missing. They work otherwise the same as before
mofa_build_model_2 <- function(data, na_features_ATAC, na_cells_ATAC, na_features_RNA, na_cells_RNA, metadata) {
  
  logcounts_allNA <- data
  logcounts_allNA[na_features_ATAC, na_cells_ATAC] <- NA
  logcounts_allNA[na_features_RNA, na_cells_RNA] <- NA
  
  mofa_list <- list(
    query = logcounts_allNA[na_features_ATAC, ],
    reference = logcounts_allNA[setdiff(1:nrow(logcounts_allNA), na_features_ATAC), ])
  
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

stab_build_model_2 <- function(data, na_features_ATAC, na_cells_ATAC, na_features_RNA, na_cells_RNA) {
  stab_list <- list(
    query = logcounts_all[setdiff(1:nrow(data), na_features_ATAC), na_cells_ATAC],
    reference = logcounts_all[setdiff(1:nrow(data), na_features_RNA),na_cells_RNA])
  
  return(stab_list)
}
#---------------------------Bridge function-------------------------------------
bridge_analysis <- function(method, bridge_size, outfile, data, metadata, na_cells_ATAC, na_features_ATAC, na_cells_RNA, all_celltypes, cluster,num_factor) {
  results <- data.frame()
  
  for (i in 1:length(bridge_size)) {
    current_bridge_size <- bridge_size[i]
    print(paste("Current bridge size:", current_bridge_size))
    
    na_features_RNA <-  1:(952 - current_bridge_size)
    
    if (method == "mofa") {
      model <- mofa_build_model_2(data = data,
                                na_features_ATAC = na_features_ATAC,
                                na_cells_ATAC = na_cells_ATAC,
                                na_features_RNA = na_features_RNA,
                                na_cells_RNA = na_cells_RNA,
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
      imp_data_ATAC <- trained_model@imputed_data$query[[1]][,na_cells_ATAC]
      imp_data_RNA <- trained_model@imputed_data$reference[[1]][,na_cells_RNA]
      
    } else if(method == "stab") {
      stab_list <- stab_build_model_2(data = data,
                                    na_features_ATAC = na_features_ATAC,
                                    na_cells_ATAC = na_cells_ATAC,
                                    na_features_RNA = na_features_RNA,
                                    na_cells_RNA = na_cells_RNA)
      
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
      imp_data_ATAC <- imp$reference[(current_bridge_size + 1):(current_bridge_size + 788),]
      imp_data_RNA <- imp$query[na_features_RNA,]
      
    } 
    
    # Common analysis steps
    sil_sum <- silhouette_summary(umap$n, umap %>% select(1:2))
    mean_sil_score <- mean(sil_sum$score)
    
    rna_train <- umap$celltype[na_cells_RNA]
    names(rna_train) <- rownames(umap)[na_cells_RNA]
    atac_query <- umap$celltype[setdiff(1:ncol(data), na_cells_RNA)]
    names(atac_query) <- rownames(umap)[setdiff(1:ncol(data), na_cells_RNA)]
    
    knn_out <- embeddingKNN(coord, rna_train, type = "uniform_fixed", k_values = 5)
    knn_acc_bal <- mean(unlist(lapply(split(isEqual(knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
    
    rmse_ATAC <- rmse_imp(imp_data = imp_data_ATAC,
                     real_data = data,
                     na_features = na_features_ATAC,
                     na_cells = na_cells_ATAC)
    
    rmse_RNA <- rmse_imp(imp_data = imp_data_RNA,
                     real_data = data,
                     na_features = na_features_RNA,
                     na_cells = na_cells_RNA)
    
    results <- rbind(results, c(current_bridge_size, mean_sil_score, knn_acc_bal, rmse_ATAC, rmse_RNA))
  }
  
  colnames(results) <- c("bridge_size", "mean_sil_score", "knn_acc_bal", "rmse_ATAC", "rmse_RNA")
  write.table(results, file = outfile, row.names = FALSE)
}


bridge_analysis("mofa", c(10, 56, 112, 167, 233, 476, 709, 942), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_RNA_out_15.txt",
                    metadata = metadata, cluster = cluster, na_cells_ATAC = na_cells_ATAC, na_cells_RNA = na_cells_RNA, na_features_ATAC = na_features_ATAC, all_celltypes = all_celltypes, data = logcounts_all, num_factor = 15)

bridge_analysis("mofa", c(10, 56, 112, 167, 233, 476, 709, 942), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_RNA_out_70.txt",
                    metadata = metadata, cluster = cluster, na_cells_ATAC = na_cells_ATAC, na_cells_RNA = na_cells_RNA, na_features_ATAC = na_features_ATAC, all_celltypes = all_celltypes, data = logcounts_all, num_factor = 70)

bridge_analysis("stab", c(10, 56, 112, 167, 233, 476, 709, 942 ), "/home/hd/hd_hd/hd_fb235/R/Data/stab_RNA_out.txt",
                    metadata = metadata, cluster = cluster, na_cells_ATAC = na_cells_ATAC, na_cells_RNA = na_cells_RNA, na_features_ATAC = na_features_ATAC, all_celltypes = all_celltypes, data = logcounts_all)
