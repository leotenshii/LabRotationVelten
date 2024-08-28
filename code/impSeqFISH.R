# Some parts of this script are adapted from
# https://github.com/MarioniLab/StabMap2021/blob/main/example/MGD_chimera_StabMap_SeqFISH_example.Rmd
# written by Shila Ghazanfar

#---------------------------Seed------------------------------------------------
set.seed(42)
#---------------------------Libraries-------------------------------------------
library(scran)
library(tidyverse)
library(MOFA2)

# library(StabMap)

source("~/R/Functions/data_prep_functions.R")
source("~/R/Functions/integration_metrics_functions.R")

#---------------------------Load Data-------------------------------------------
scrna <- readRDS("~/R/Data/MGD_chimera_scrna.Rds")
seqFISH <- readRDS("~/R/Data/seqFISH.Rds")


# Feature selection for scRNA-seq data
decomp <- modelGeneVar(scrna, block = scrna$sample)
hvgs <- rownames(decomp)[decomp$mean>0.01 & decomp$p.value <= 0.05]
length(hvgs) #3399
scrna <- scrna[hvgs,]

rownames(scrna) <- make.unique(rowData(scrna)$SYMBOL)

#------------- StabMap----------------------------------------------------------
# Not needed for vignette
# assay_list = list(
#   seqFISH = seqFISH,
#   scrna = as.matrix(logcounts(scrna))
# )
# 
# 
# stab = stabMap(assay_list,
#                  reference_list = c("scrna", "seqFISH"),
#                  projectAll = TRUE,
#                  plot = FALSE,
#                  scale.center = TRUE,
#                  scale.scale = TRUE)
# 
# 
# imp = imputeEmbedding(
#   assay_list,
#   stab,
#   reference = colnames(assay_list[["scrna"]]),
#   query = colnames(assay_list[["seqFISH"]]))
# 
# 
# saveRDS(imp, "~/R/Data/impSeqFISH.RDS")

# impSeqFISH <- readRDS("~/R/Data/impSeqFISH.RDS")
# impSeqFISH_embryo1 <- impSeqFISH$scrna[, grepl("^embryo1", colnames(impSeqFISH$scrna))]
# 
# impSeqFISH_embryo1_long <- pivot_longer(as.data.frame(impSeqFISH_embryo1) %>% rownames_to_column("gene"), cols = -gene, names_to = "uniqueID", values_to = "expr")
# 
# impSeqFISH_embryo1_joined <- inner_join(impSeqFISH_embryo1_long, as.data.frame(meta))
# 
# ggplot(impSeqFISH_embryo1_joined %>% filter(gene == "T") %>% mutate(expr = ifelse(celltype_mapped_refined == "Splanchnic mesoderm", expr, NA)), aes(x=x_global, y=y_global, col = expr)) + geom_point( size = 0.5) + theme_void()         

#---------------------------MOFA------------------------------------------------
# Put samples in two different groups
sample_group <- c(rep("seqFISH", ncol(seqFISH)), rep("scrna", ncol(scrna)))

# Combine scRNA and seqFISH tables
df1 <- as.data.frame(seqFISH)
df2 <- as.data.frame(as.matrix(logcounts(scrna)))

df1$Feature <- rownames(df1)
df2$Feature <- rownames(df2)

combined_df <- merge(df1, df2, by = "Feature", all = TRUE)

rownames(combined_df) <- combined_df$Feature
combined_df$Feature <- NULL

combined_matrix <- as.matrix(combined_df)

# Create and train MOFA model
mofa_model <- create_mofa(list(combined_matrix), 
                          group_names = sample_group)

trained_model <- mofa_parameter_train(num_factors = 70, model = mofa_model, spikeslab_weights = FALSE) # takes a few hours and ~35GB RAM
trained_model <- impute(trained_model)

saveRDS(trained_model, "~/R/Data/trained_modelSeqFISH.RDS")
#-------------------------------------------------------------------------------
meta <- colData(seqFISH)

trained_modelSeqFISH <- readRDS("~/R/Data/trained_modelSeqFISH.RDS")

trained_modelSeqFISH_long <- trained_modelSeqFISH@imputed_data$view_1$group1 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "uniqueID", values_to = "expr") %>%
  inner_join(as.data.frame(meta), by = "uniqueID")

# Filter data and plot
filtered_data <- trained_modelSeqFISH_long %>%
  filter(gene == "T") %>%
  mutate(expr = ifelse(celltype_mapped_refined == "Splanchnic mesoderm", expr, NA))

ggplot(filtered_data, aes(x = x_global, y = -y_global, color = expr)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(na.value = "lightgray") +
  theme_void()

