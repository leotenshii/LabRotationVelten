# This script is partly adapted from 
# https://github.com/MarioniLab/StabMap2021/blob/main/example/MGD_Multihop_example.Rmd
# written by Shila Ghazanfar
#---------------------------Seed------------------------------------------------
set.seed(42)
#---------------------------Libraries-------------------------------------------
library(scater)
library(scran)
library(StabMap)
library(MOFA2)
library(shades)
library(tidyverse)
library(patchwork)
library(ggpubr)

source("~/R/Functions/adaptiveKNN.R")
source("~/R/Functions/data_prep_functions.R")

#---------------------------Data------------------------------------------------
# Data can be loaded from the MouseGastrulationData Bioconductor package. I wasn't able to install it from the server, therefore loading it from RDS
# Loading it from MouseGastrulationData:

# library(MouseGastrulationData)
# stages <- c("E8.5")
# mt <- MouseGastrulationData::AtlasSampleMetadata
# samples <- mt[mt[, "stage"] %in% stages, "sample"]
# atlas_raw <- EmbryoAtlasData(type = "processed", samples = samples)

# Load from RDS
atlas_raw = readRDS("~/R/Data/atlas.RDS")
atlas <- logNormCounts(atlas_raw)

# Feature selection
decomp <- modelGeneVar(atlas)
hvgs <- rownames(decomp)[decomp$mean > 0.05 & decomp$total > 0]
length(hvgs)
atlas <- atlas[hvgs,]

#---------------------------Simulation------------------------------------------
# Build simulation set up. Take the dataset and split into k groups.
# Then select up to G/k features from the HVGs and determine the p percent overlap
# between contiguous datasets. Then jointly embed the data using StabMap or MOFA, and 
# record the cell type classification accuracy.

k = 8
nGenes_all = c(100, 200, 500, 1000)
nCells_all = c(500, 1000, 2000)
labels = "celltype"

stab_acc_df_all = NULL
mofa_acc_df_all = NULL

for (i in 1:5) {
  
  for (nCells in nCells_all) {
    
    for (nGenes in nGenes_all) {
      
      print(i)
      print(nCells)
      print(nGenes)
      
      # check there are enough features, this should be FALSE
      ceiling((k+1)*(nGenes/2)) > length(hvgs)
      
      # Splitting features
      features_to_split = sample(hvgs, ceiling((k+1)*(nGenes/2)))
      feature_split_all = split(features_to_split, rep(seq_len(k+1), length.out = length(features_to_split)))
      feature_split = list()
      for (kk in 1:k) {
        feature_split[[kk]] <- unlist(feature_split_all[kk:(kk+1)])
      }
      names(feature_split) <- paste0("Dataset_",seq_len(k))
      
      # Splitting cells
      atlas_cells = sample(colnames(atlas))[seq_len(nCells*8)]
      cell_split = split(atlas_cells, rep(seq_len(k), length.out = length(atlas_cells)))
      
      # StabMap
      assay_list = mapply(function(feat,cel) {
        as.matrix(logcounts(atlas)[feat,cel])
      }, feature_split, cell_split, SIMPLIFY = FALSE)

      type = setNames(rep(paste0("Dataset_",names(cell_split)), times = unlist(lapply(cell_split, length))), unlist(cell_split))
      
      nPCs = ifelse(nGenes <=50, 10, 50)
      
      StabMap_embedding = stabMap(assay_list,
                                  reference_list = "Dataset_1",
                                  ncomponentsReference = nPCs,
                                  ncomponentsSubset = nPCs,
                                  projectAll = TRUE,
                                  scale.center = FALSE,
                                  scale.scale = FALSE)
      dim(StabMap_embedding)
      
      # now perform cell type classification using the embedding
      # predict cell type labels using knn with k = 5
      referenceLabels = colData(atlas)[colnames(assay_list[["Dataset_1"]]),labels]
      names(referenceLabels) = colnames(assay_list[["Dataset_1"]])
      
      queryLabels = colData(atlas)[,labels]
      names(queryLabels) = colnames(atlas)
      
      # Celltype accuracy
      # only use cells with labels
      
      data_all = StabMap_embedding
      labels_train = referenceLabels[!is.na(referenceLabels) & names(referenceLabels) %in% rownames(data_all)]
      
      stab_knn_out = embeddingKNN(data_all,
                             labels_train,
                             type = "uniform_fixed",
                             k_values = 5)
      
      stab_acc = isEqual(stab_knn_out[names(queryLabels),"predicted_labels"], queryLabels)
      
      stab_prop.acc = prop.table(table(stab_acc, type[names(queryLabels)]), 2)["1",]
      stab_prop.acc[gtools::mixedorder(names(stab_prop.acc))]
      barplot(stab_prop.acc[gtools::mixedorder(names(stab_prop.acc))][-1])
      
      # Put results in dataframe
      stab_acc_df = data.frame(nCells = nCells,
                          nGenes = nGenes,
                          acc = stab_prop.acc,
                          Dataset = names(stab_prop.acc),
                          rep = paste0(c(nCells,nGenes,i), collapse = "_"))
      
      stab_acc_df_all <- rbind(stab_acc_df_all, stab_acc_df)
      
      
      # MOFA
      combined <- bind_rows(
        lapply(names(assay_list), function(name) {
          df <- as.data.frame(assay_list[[name]])
          return(df)
        }))
      
      a <- list()
      
      # Loop to slice the combined dataframe into matrices of 500 rows each
      for (i in 1:k) {
        # Extract rows from (i*500 - 499) to (i*500) and convert to a matrix
        matrix_slice <- as.matrix(combined[((i - 1) * nGenes + 1):(i * nGenes), ])
        
        # Append the matrix slice to the list
        a[[i]] <- matrix_slice
      }
      
      model <- create_mofa(a)
      
      mofa_parameter_train(num_factors = ifelse(nGenes <=50, 10, 50), model = model, spikeslab_weights = FALSE, path = "/home/hd/hd_hd/hd_fb235/R/Data/multihop_model2.hdf5")
      trained_model <- load_model("/home/hd/hd_hd/hd_fb235/R/Data/multihop_model2.hdf5", remove_inactive_factors = FALSE)
      
      # Celltype accuracy
      factors <- get_factors(trained_model)
      
      mofa_knn_out = embeddingKNN(factors$group1,
                                  labels_train,
                                  type = "uniform_fixed",
                                  k_values = 5)
      
      mofa_acc = isEqual(mofa_knn_out[names(queryLabels),"predicted_labels"], queryLabels)
      
      mofa_prop.acc = prop.table(table(mofa_acc, type[names(queryLabels)]), 2)["1",]
      mofa_prop.acc[gtools::mixedorder(names(mofa_prop.acc))]
      barplot(mofa_prop.acc[gtools::mixedorder(names(mofa_prop.acc))][-1])
      
      # Put results in dataframe
      mofa_acc_df = data.frame(nCells = nCells,
                               nGenes = nGenes,
                               acc = mofa_prop.acc,
                               Dataset = names(mofa_prop.acc),
                               rep = paste0(c(nCells,nGenes,i), collapse = "_"))
      
      mofa_acc_df_all <- rbind(mofa_acc_df_all, mofa_acc_df)
    }
    
  }
  
}

#---------------------------Save results----------------------------------------
write.table(mofa_acc_df_all, file = "~/R/Data/mofa_multihop3.txt")
#write.table(stab_acc_df_all, file = "~/R/Data/stab_multihop2.txt")

stab_acc_df_all <- read.table("~/R/Data/stab_multihop.txt")
mofa_acc_df_all <- read.table("~/R/Data/mofa_multihop2.txt")

# stab_acc_df_all <- read.table("~/R/Data/stab_multihop2.txt")
# mofa_acc_df_all <- read.table("~/R/Data/mofa_multihop2.txt")

#---------------------------Plots-----------------------------------------------
g1 = ggplot(subset(mofa_acc_df_all, Dataset != "Dataset_1"), aes(x = Dataset, y = acc)) + 
  geom_line(aes(colour = factor(nGenes), group = rep), alpha = 0.1) + 
  geom_point(aes(colour = factor(nGenes), group = rep), alpha = 0.1) + 
  geom_smooth(aes(fill = factor(nGenes), group = factor(nGenes), colour = factor(nGenes)), alpha = 0.3) +
  facet_wrap(~nCells, nrow = 1, labeller = labeller(nCells = function(x) paste0(x, " cells per dataset"))) + 
  ylab("Accuracy") +
  xlab("Dataset") +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  lightness(scale_colour_brewer(palette = 1, aesthetics = c("colour", "fill")), scalefac(0.7)) +
  NULL

g2 = ggplot(subset(stab_acc_df_all, Dataset != "Dataset_1"), aes(x = Dataset, y = acc)) + 
  geom_line(aes(colour = factor(nGenes), group = rep), alpha = 0.1) + 
  geom_point(aes(colour = factor(nGenes), group = rep), alpha = 0.1) + 
  geom_smooth(aes(fill = factor(nGenes), group = factor(nGenes), colour = factor(nGenes)), alpha = 0.3) +
  facet_wrap(~nCells, nrow = 1, labeller = labeller(nCells = function(x) paste0(x, " cells per dataset"))) + 
  ylab("Accuracy") +
  xlab("Dataset") +
  theme_classic() + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  lightness(scale_colour_brewer(palette = 1, aesthetics = c("colour", "fill")), scalefac(0.7)) +
  NULL

g_leg = as_ggplot(get_legend(g1 + theme(legend.position = "right") + guides(colour = guide_legend(title = "Number of genes \n per dataset"),
                                                                           fill = guide_legend(title = "Number of genes \n per dataset"))))
print(g1+g2+ g_leg) + plot_layout(widths = c(2,2,0.5))


#---------------------------Runs------------------------------------------------
# Run 1:
# StabMap_embedding = stabMap(assay_list,
#                             reference_list = "Dataset_1",
#                             ncomponentsReference = nPCs,
#                             ncomponentsSubset = nPCs,
#                             projectAll = TRUE)
# 
# mofa_parameter_train(model = model, spikeslab_weights = FALSE, path = "/home/hd/hd_hd/hd_fb235/R/Data/multihop_model.hdf5")

# Run 2:
# StabMap_embedding = stabMap(assay_list,
#                             reference_list = "Dataset_1",
#                             ncomponentsReference = nPCs,
#                             ncomponentsSubset = nPCs,
#                             projectAll = TRUE,
#                             scale.center = FALSE,
#                             scale.scale = FALSE)
# 
# mofa_parameter_train(num_factors = ifelse(nGenes <=50, 10, 50), model = model, spikeslab_weights = FALSE, path = "/home/hd/hd_hd/hd_fb235/R/Data/multihop_model.hdf5")

# Run 3 with list of matices in mofa


