#---------------------------Seed------------------------------------------------
set.seed(42)

#---------------------------Parameters------------------------------------------
param_grid <- expand.grid(
  ncomponentsReference = c( 50),
  ncomponentsSubset = c(40, 50),
  maxFeatures = c(1000),
  scale.center = c(TRUE),
  scale.scale = c(TRUE)
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
  
  US1 <- ggplot(stab_umap) +
    geom_point(aes(x = V1, y = V2, color = celltype), size = .1) +
    theme_light() +
    ggtitle("StabMap celltype")
  US2 <- ggplot(stab_umap) +
    geom_point(aes(x = V1, y = V2, color = isNA), size = .1) +
    theme_light() +
    ggtitle("StabMap NA")
  US3 <- ggplot(stab_umap) +
    geom_point(aes(x = V1, y = V2, color = broad_celltype), size = .1) +
    theme_light() +
    ggtitle("StabMap broad celltype")
  
  # Calculate silhouette scores
  stab_sil <- silhouette(stab_umap$k, dist(stab_umap %>% select(V1, V2)))
  stab_sil_sum <- stab_sil %>%
    as.data.frame() %>%
    group_by(cluster) %>%
    summarise(score = mean(sil_width), frac_pos = sum(sil_width > 0) / n(), pos_score = sum((sil_width > 0) * sil_width) / sum(sil_width > 0))
  
  # Calculate clustering stats
  stab_cluster_stats <- cluster.stats(dist(stab_umap %>% select(V1, V2)), stab_umap$k)
  
  # Cell Type Accuracy
  stab_knn_out = embeddingKNN(stab_umap_coord,
                              rna_train,
                              type = "uniform_fixed",
                              k_values = 5)
  stab_knn_acc = mean(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), na.rm = TRUE)
  stab_knn_acc_bal = mean(unlist(lapply(split(isEqual(stab_knn_out[names(atac_query),"predicted_labels"], atac_query), atac_query), mean, na.rm = TRUE)))
  
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
    US1 = US1,
    US2 = US2,
    US3 = US3, 
    stab_knn_acc = stab_knn_acc,
    stab_knn_acc_bal = stab_knn_acc_bal
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
  stab_knn_acc_bal = numeric()
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
    stab_knn_acc_bal = stab_knn_acc_bal
  ))
}

print(comparison)
