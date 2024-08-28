# Last changes on 28.08.2024
# Author: Leoni Zimmermann

# Description ------------------------------------------------------------------
# The functions are used for calculating integration and imputation metrics.
# silhouette_summary returns per cluster the mean silhouette width, proportion of observations that are well-clustered (having a positive silhouette width) 
#       and the average of the positive silhouette widths, highlighting the quality of clustering among well-clustered points only 
# rmse_imp returns the RMSE calculated from the imputated and real data
# rmse_imp_table returns the table with the predicted (imputed) and actual values. 
#       I use this to go from the table to two seperate RMSE calculation for different features (e.g. ATAC RMSE and RNA RMSE). Could also have been integrated in rmse_imp but I just didn't :)

# Functions --------------------------------------------------------------------
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

