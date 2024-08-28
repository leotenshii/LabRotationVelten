# Last changes on 28.08.2024
# Author: Leoni Zimmermann

#---------------------------Description-----------------------------------------
# In this script, cells with missing features are randomly sampled resulting in five different sets of each 5016 cells missing ATAC features. 
# The data is then processed with MOFA and StabMap to investigate whether the presence of missing features affects the integration and imputation metrics.
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

# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all)))
colnames(all_celltypes) <- "celltype"

# Clusters
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

na_features <- 953:1740

na_cells_list <- list()

for (j in 1:5) {  
  na_cells_rand <- sample(10032, 5016)
  na_cells_list <- append(na_cells_list, list(na_cells_rand))
}

#---------------------------Function--------------------------------------------

# Function to run the MOFA analysis
run_NAcells_analysis <- function(method, num_factors, outfile) {
  results <- data.frame()
  
  for (i in 1:length(na_cells_list)) {
    
    na_cells <- na_cells_list[[i]]

    if (method == "mofa") {
      # Create and train model
      model <- mofa_build_model(data = logcounts_all, 
                                metadata = metadata, 
                                na_cells = na_cells, 
                                na_features = na_features)
      
      trained_model <- mofa_parameter_train(num_factors = num_factors, 
                                            model = model, 
                                            spikeslab_weights = FALSE)
      
      # UMAP
      trained_model <- run_umap(trained_model)
      umap_coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
      
      # Adding metadata 
      umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0)  %>%
        full_join(cluster) %>% 
        column_to_rownames("Row.names") %>%
        select(-sample)
      
      # Imputation
      trained_model <- impute(trained_model)
      imp_data <- trained_model@imputed_data$missing_feat[[1]][, na_cells]
    
    } else if(method == "stab") {
      
      # Create and train model
      assay_list <- stab_build_model(data = logcounts_all, 
                                     na_features = na_features, 
                                     na_cells = na_cells)
      
      stab = stabMap(assay_list,
                     ncomponentsReference = num_factors,
                     ncomponentsSubset = num_factors,
                     reference_list = c("all_feat"),
                     plot = FALSE,
                     scale.center = FALSE,
                     scale.scale = FALSE)
      
      # UMAP
      umap_coord <- as.data.frame(calculateUMAP(t(stab)))
      
      # Adding metadata
      umap <- merge( as.data.frame(metadata), umap_coord, by =0 ) %>%
        full_join(cluster) %>%
        column_to_rownames("Row.names")
      
      # Imputation
      imp = imputeEmbedding(
        assay_list,
        stab,
        reference = colnames(assay_list[["all_feat"]]),
        query = colnames(assay_list[["missing_feat"]]))
      
      imp_data <- imp$all_feat[na_features,]
      
    }
    
    # Silhouette score
    sil_sum <- silhouette_summary(umap$n, umap %>% select(UMAP1, UMAP2))
    mean_sil_score <- mean(sil_sum$score)
    
    # Cell type accuracy
    rna_train <- umap$celltype[-na_cells]
    names(rna_train) <- rownames(umap)[-na_cells]
    atac_query <- umap$celltype[na_cells]
    names(atac_query) <- rownames(umap)[na_cells]
    
    knn_out <- embeddingKNN(umap_coord, rna_train, type = "uniform_fixed", k_values = 5)
    knn_acc_bal <- mean(unlist(lapply(split(isEqual(knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
    
    # RMSE calculation
    rmse<- rmse_imp(imp_data = imp_data,
                         real_data = logcounts_all, 
                         na_features = na_features, 
                         na_cells = na_cells)
    
    # Put results together
    results <- rbind(results, c(mean_sil_score, knn_acc_bal, rmse))
  }
  
  colnames(results) <- c("mean_sil_score", "knn_acc_bal", "rmse")
  write.table(results, file = outfile, row.names = FALSE)
}

run_NAcells_analysis("mofa", 15, "/home/hd/hd_hd/hd_fb235/R/Data/mofa_random_cells_NA_15.txt")
run_NAcells_analysis("mofa", 70, "/home/hd/hd_hd/hd_fb235/R/Data/mofa_random_cells_NA_70.txt")

run_NAcells_analysis("stab", 70, "/home/hd/hd_hd/hd_fb235/R/Data/r.txt")





  
