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


source("/home/hd/hd_hd/hd_fb235/R/Scripts/adaptiveKNN.R")


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
sce.atac <- normalize_and_select_features(experiments(mae)[["atac"]], 0.01, 0.05)

logcounts_all <- rbind(logcounts(sce.rna), logcounts(sce.atac))


#---------------------------General Data----------------------------------------

# Celltypes of all samples
all_celltypes <- as.data.frame(setNames(metadata$celltype, colnames(logcounts_all))) %>%
  rename( celltype = "setNames(metadata$celltype, colnames(logcounts_all))")

# Clusters as numbers
cluster <- as.data.frame(metadata) %>% group_by(celltype) %>% summarise(n = n()) %>% mutate(k = 1:14) %>% select(-n)

#---------------------------MOFA------------------------------------------------
results_mofa <- data.frame(row.names = c("bridge_size", "mean_sil_score", "mofa_knn_acc", "mofa_rsme"))


# Function to run the MOFA analysis
run_mofa_analysis <- function(bridge_size, outfile) {
  results <- data.frame()
  
  for (i in 1:length(bridge_size)) {
    print(length(bridge_size))
    current_bridge_size <- bridge_size[i]  # Use a different variable name
    print(current_bridge_size)
    na_features <- sample(1740, current_bridge_size)
    
    logcounts_allNA <- as.matrix(logcounts_all)
    logcounts_allNA[na_features, 1:5016] <- NA
    
    # Create list for MOFA
    mofa_list <- list(
      with_na_features = logcounts_allNA[na_features, ],
      no_na_features = logcounts_allNA[-na_features, ]
    )
    
    # MOFA model and training
    model <- create_mofa(mofa_list)
    samples_metadata(model) <- as.data.frame(metadata) %>% 
      rownames_to_column("sample")
    plot_data_overview(model)
    
    model_opts <- get_default_model_options(model)
    model_opts$num_factors <- 15
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
    
    # MOFA celltype
    mofa_sil_sum <- silhoutte_summary(mofa_umap$k, mofa_umap %>% select(UMAP1, UMAP2))
    mean_sil_score <- mean(mofa_sil_sum$score)
    
    # Predict "ATAC only" cell's cell types
    rna_train <- mofa_umap$celltype[5016:10032]
    names(rna_train) <- rownames(mofa_umap)[5016:10032]
    atac_query <- mofa_umap$celltype[1:5016]
    names(atac_query) <- rownames(mofa_umap)[1:5016]
    
    mofa_knn_out <- embeddingKNN(mofa_umap_coord, rna_train, type = "uniform_fixed", k_values = 5)
    mofa_knn_acc <- mean(isEqual(mofa_knn_out[names(atac_query), "predicted_labels"], atac_query), na.rm = TRUE)
    
    # Imputation
    trained_model <- impute(trained_model)
    
    mofa_imp_comp <- as.data.frame(trained_model@imputed_data$with_na_features[[1]][, 1:5016]) %>%
      rownames_to_column("feature") %>%
      pivot_longer(-feature, names_to = "sample", values_to = "predicted") %>%
      full_join(as.data.frame(as.matrix(logcounts_all[na_features, 1:5016])) %>%
                  rownames_to_column("feature") %>%
                  pivot_longer(-feature, names_to = "sample", values_to = "actual"))
    
    mofa_rsme <- sqrt(mean((mofa_imp_comp$actual - mofa_imp_comp$predicted)^2))
    
    results <- rbind(results, c(1740 - current_bridge_size, mean_sil_score, mofa_knn_acc, mofa_rsme))
  }
  
  colnames(results) <- c("bridge_size", "mean_sil_score", "mofa_knn_acc", "mofa_rsme")
  write.table(results, file = outfile, row.names = FALSE)
}

# Run the analysis for num_factors = 70
run_mofa_analysis(c(1305, 870), "/home/hd/hd_hd/hd_fb235/R/Data/mofa_bridge.txt")


# Run the analysis for num_factors = 15
# run_mofa_analysis(15, "/home/hd/hd_hd/hd_fb235/R/Data/mofa_random_cells_NA_15.txt")




#---------------------------StabMap---------------------------------------------
results_stab <- data.frame(row.names = c("loop", "mean_sil_score", "stab_knn_acc", "stab_rsme"))

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
  
  stab_sil <- silhouette(stab_umap$k, dist(stab_umap %>% select(V1, V2)))
  stab_sil_sum <- stab_sil %>% 
    as.data.frame() %>% group_by(cluster) %>% summarise(score = mean(sil_width), 
                                                        frac_pos = sum(sil_width > 0)/n(),
                                                        pos_score = sum((sil_width>0)*sil_width)/sum(sil_width > 0)) #average quality of clustering for the well-clustered data points within each cluster
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
  stab_knn_acc = mean(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), na.rm = TRUE)
  
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
  
  stab_rsme <- sqrt(mean((stab_imp_comp$actual - stab_imp_comp$predicted)^2))
  
  results_stab <- rbind(results_stab,c(mean_sil_score, stab_knn_acc, stab_rsme))
}

colnames(results_stab) <- c("mean_sil_score", "stab_knn_acc", "stab_rsme")

write.table(results_stab, file = "/home/hd/hd_hd/hd_fb235/R/Data/stab_random_cells_NA.txt")


