# Libraries ---------------------------------------------------------------
library(scuttle)
library(scran)


# Functions ---------------------------------------------------------------


normalize_and_select_features <- function(data, expr_mean_threshold, p_value_threshold) {
  # The function log normalizes the input data and selects highly variable genes on given thresholds

  data <- logNormCounts(data)
  decomp <- modelGeneVar(data)
  hvgs <- rownames(decomp)[decomp$mean > expr_mean_threshold & decomp$p.value <= p_value_threshold]
  data <- data[hvgs,]
  return(data)
}

