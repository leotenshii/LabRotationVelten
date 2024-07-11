# Libraries ---------------------------------------------------------------
library(cluster)


# Functions ---------------------------------------------------------------

silhoutte_summary <- function(cluster_numbers, coordinates) {
  sil <- silhouette(cluster_numbers, dist(coordinates))
  sil_sum <- sil %>% 
    as.data.frame() %>% group_by(cluster) %>% summarise(score = mean(sil_width), 
                                                        frac_pos = sum(sil_width > 0)/n(),
                                                        pos_score = sum((sil_width>0)*sil_width)/sum(sil_width > 0)) 
  return(sil_sum)
}
