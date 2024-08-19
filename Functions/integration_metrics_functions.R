# Functions ---------------------------------------------------------------


silhouette_summary <- function(cluster_numbers, coordinates) {
  sil <- silhouette(cluster_numbers, dist(coordinates))
  sil_sum <- sil %>% 
    as.data.frame() %>% group_by(cluster) %>% summarise(score = mean(sil_width), 
                                                        frac_pos = sum(sil_width > 0)/n(),
                                                        pos_score = sum((sil_width>0)*sil_width)/sum(sil_width > 0)) 
  return(sil_sum)
}

rmse_imp <- function(imp_data, real_data, na_features, na_cells) {

  imp_data_df <- as.data.frame(imp_data)
  
  real_data_df <- as.data.frame(real_data[na_features, na_cells])
  
  imp_data_long <- imp_data_df %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted")
  
  real_data_long <- real_data_df %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "actual")
  
  comp <- full_join(imp_data_long, real_data_long, by = c("feature", "sample"))
  
  rmse <- sqrt(mean((comp$actual - comp$predicted)^2, na.rm = TRUE))
  
  return(rmse)
}

rmse_imp_table <- function(imp_data, real_data, na_features, na_cells) {
  
  imp_data_df <- as.data.frame(imp_data)
  
  real_data_df <- as.data.frame(real_data[na_features, na_cells])
  
  imp_data_long <- imp_data_df %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "predicted")
  
  real_data_long <- real_data_df %>%
    rownames_to_column("feature") %>%
    pivot_longer(-feature, names_to = "sample", values_to = "actual")
  
  comp <- full_join(imp_data_long, real_data_long, by = c("feature", "sample"))
  
  return(comp)
}

