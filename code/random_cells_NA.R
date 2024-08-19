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


# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename(celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters
cluster <- as.data.frame(metadata) %>% 
  group_by(celltype) %>% 
  summarise(n = n()) %>% 
  mutate(n = 1:14) 

#---------------------------MOFA------------------------------------------------
results_mofa <- data.frame(row.names = c("loop", "mean_sil_score", "mofa_knn_acc_bal", "mofa_rmse"))


# Function to run the MOFA analysis
run_mofa_analysis <- function(num_factors, outfile) {
  results <- data.frame()
  
  for (i in 1:5) {
    NA_cells <- sample(10032, 5016)
    # Put NAs in data
    logcounts_all[953:1740, NA_cells] <- NA
    
    # Create list for MOFA
    mofa_list <- list(
      RNA = logcounts_all[1:952, ],
      ATAC = logcounts_all[953:1740, ]
    )
    
    # MOFA model and training
    model <- create_mofa(mofa_list)
    samples_metadata(model) <- as.data.frame(metadata) %>% 
      rownames_to_column("sample")
    plot_data_overview(model)
    
    model_opts <- get_default_model_options(model)
    model_opts$num_factors <- num_factors
    MOFAobject <- prepare_mofa(model, model_options = model_opts)
    
    # Integration
    run_mofa(MOFAobject, outfile = "R/Data/model.hdf5", use_basilisk = TRUE)
    trained_model <- load_model("R/Data/model.hdf5")
    
    # UMAP
    trained_model <- run_umap(trained_model)
    mofa_umap_coord <- trained_model@dim_red$UMAP %>% select(UMAP1, UMAP2)
    
    # Add metadata to the UMAP results
    mofa_umap <- merge(trained_model@dim_red$UMAP, all_celltypes, by = 0)  %>%
      full_join(cluster) %>% 
      column_to_rownames("Row.names") %>%
      select(-sample)
    
    #MOFA celltype
    mofa_sil_sum <- silhouette_summary(mofa_umap$n, mofa_umap %>% select(UMAP1, UMAP2))
    
    mean_sil_score <- mean(mofa_sil_sum$score)
    
    # Predict "ATAC only" cell's cell types
    rna_train <- mofa_umap$celltype[-NA_cells]
    names(rna_train) <- rownames(mofa_umap)[-NA_cells]
    atac_query <- mofa_umap$celltype[NA_cells]
    names(atac_query) <- rownames(mofa_umap)[NA_cells]
    
    mofa_knn_out <- embeddingKNN(mofa_umap_coord, rna_train, type = "uniform_fixed", k_values = 5)
    mofa_knn_acc_bal <- mean(unlist(lapply(split(isEqual(mofa_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
    
    # Imputation
    trained_model <- impute(trained_model)
    
    mofa_imp_comp <- as.data.frame(trained_model@imputed_data$ATAC[[1]][, NA_cells]) %>%
      rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
      full_join(as.data.frame(as.matrix(logcounts_all[953:1740, NA_cells])) %>%
                  rownames_to_column("feature") %>%
                  pivot_longer(-feature, names_to = "sample", values_to = "actual"))
    
    mofa_rmse <- sqrt(mean((mofa_imp_comp$actual - mofa_imp_comp$predicted)^2))
    
    results <- rbind(results, c(mean_sil_score, mofa_knn_acc_bal, mofa_rmse))
  }
  
  colnames(results) <- c("mean_sil_score", "mofa_knn_acc_bal", "mofa_rmse")
  write.table(results, file = outfile, row.names = FALSE)
}

# Run the analysis for num_factors = 70
# run_mofa_analysis(70, "/home/hd/hd_hd/hd_fb235/R/Data/mofa_random_cells_NA_70.txt")

# Run the analysis for num_factors = 15
# run_mofa_analysis(15, "/home/hd/hd_hd/hd_fb235/R/Data/mofa_random_cells_NA_15.txt")




#---------------------------StabMap---------------------------------------------
results_stab <- data.frame(row.names = c("loop", "mean_sil_score", "mofa_knn_acc_bal", "stab_rmse"))

# Seperation ATAC Multiome
names <- c(rep("ATAC", ncol(logcounts_all)/2), rep("Multiome", ncol(logcounts_all)/2))

# Assay Types
  assayType = ifelse(rownames(logcounts_all) %in% rownames(sce.rna),
                     "rna", "atac")
for (i in 1:5) {
  names <- sample(names)
  
  # List for StabMap
  assay_list = list(
    ATAC = logcounts_all[assayType %in% c("rna"), names %in% c("ATAC")], # 952x5016 -> only RNA in first half of cells
    Multiome = logcounts_all[assayType %in% c("rna", "atac"), names %in% c("Multiome")] #1740x5016 -> ATAC and RNA, second half of cells
  )
  
  
  # StabMap
  mosaicDataUpSet(assay_list)
  stab = stabMap(assay_list,
                 ncomponentsReference = 70,
                 ncomponentsSubset = 70,
                 reference_list = c("Multiome"),
                 plot = FALSE,
                 scale.center = FALSE,
                 scale.scale = FALSE)
  
  # UMAP
  stab_umap_coord <- as.data.frame(calculateUMAP(t(stab)))
  
  
  # Add metadata to the UMAP results (celltype, if it was an NA cell, cluster)
  stab_umap <- merge( as.data.frame(metadata), stab_umap_coord, by =0 ) %>%
    full_join(cluster) %>%
    column_to_rownames("Row.names")
  
  #StabMap celltype
  stab_sil_sum <- silhouette_summary(stab_umap$n, stab_umap %>% select(V1, V2))

  mean_sil_score <- mean(stab_sil_sum$score)
  
  rna_train <- stab_umap$celltype[names == "Multiome"]
  names(rna_train) <- rownames(stab_umap)[names == "Multiome"]
  atac_query <- stab_umap$celltype[names != "Multiome"]
  names(atac_query) <- rownames(stab_umap)[names != "Multiome"]
  
  # StabMap
  stab_knn_out = embeddingKNN(stab_umap_coord,
                              rna_train,
                              type = "uniform_fixed",
                              k_values = 5)
  stab_knn_acc_bal <- mean(unlist(lapply(split(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))

  #---------------------------Imputation------------------------------------------
  imp = imputeEmbedding(
    assay_list,
    stab,
    reference = colnames(assay_list[["Multiome"]]),
    query = colnames(assay_list[["ATAC"]]))
  
  stab_imp_comp <- as.data.frame(as.matrix(imp$Multiome)[953:1740,]) %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
    full_join(as.data.frame(as.matrix(logcounts_all[953:1740, names != "Multiome"])) %>%
                rownames_to_column("feature") %>%
                pivot_longer(-feature, names_to = "sample", values_to = "actual"))
  
  stab_rmse <- sqrt(mean((stab_imp_comp$actual - stab_imp_comp$predicted)^2))
  
  results_stab <- rbind(results_stab,c(mean_sil_score, stab_knn_acc_bal, stab_rmse))
}
  
  colnames(results_stab) <- c("mean_sil_score", "stab_knn_acc_bal", "stab_rmse")
  
  write.table(results_stab, file = "/home/hd/hd_hd/hd_fb235/R/Data/stab_random_cells_NA.txt")
  

