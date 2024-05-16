#---------------------------Seed------------------------------------------------
set.seed(42)
#---------------------------Libraries-------------------------------------------
library(MOFA2)
library(scater)
library(scran)
library(StabMap)
library(tidyverse)
library(cluster)

#library(SingleCellMultiModal)

#---------------------------Dataset---------------------------------------------
# Peripheral Blood Mononuclear Cells provided by 10x Genomics website
# 10x Genomics Multiome technology enables simultaneous profiling of the transcriptome 
# (using 3’ gene expression) and epigenome (using ATAC-seq) from single cells to deepen 
# our understanding of how genes are expressed and regulated across different cell types.

#mae <- scMultiome("pbmc_10x", mode = "*", dry.run = FALSE, format = "MTX")

# Loaded from RDS cause I cant install the SingleCellMultiModal package
mae <- readRDS("R/Data/data.RDS")
metadata <- mae@colData
metadata_ct <- metadata$celltype


#---------------------------Data preperation------------------------------------
# Normalization RNA
sce.rna <- experiments(mae)[["rna"]]

# Normalisation
sce.rna <- logNormCounts(sce.rna)

# Feature selection
decomp <- modelGeneVar(sce.rna)
hvgs <- rownames(decomp)[decomp$mean>0.01 & decomp$p.value <= 0.05]

sce.rna <- sce.rna[hvgs,]


# Normalization ATAC
sce.atac <- experiments(mae)[["atac"]]

# Normalise
sce.atac <- logNormCounts(sce.atac)

# Feature selection using highly variable peaks
# And adding matching peaks to genes
decomp <- modelGeneVar(sce.atac)
hvgs <- rownames(decomp)[decomp$mean>0.25
                         & decomp$p.value <= 0.05]

sce.atac <- sce.atac[hvgs,]

logcounts_all = rbind(logcounts(sce.rna), logcounts(sce.atac))

#---------------------------MOFA------------------------------------------------

# Put NAs in data 
logcounts_allNA <- logcounts_all

#logcounts_allNA[953:1740, 1:ncol(logcounts_all)/2] <- NA

logcounts_allNA <- readRDS("R/Data/logcountsNA.RDS")

# Create list for MOFa
mofa_list <- list(
  RNA =  logcounts_allNA[1:952,],
  ATAC = logcounts_allNA[953: 1740,]
)

# MOFA model and training
model <- create_mofa(mofa_list)
samples_metadata(model) <- as.data.frame(metadata) %>% rownames_to_column("sample")
plot_data_overview(model)
MOFAobject <- prepare_mofa(model)

run_mofa(MOFAobject, outfile = "~/R/Data/model.hdf5", use_basilisk = TRUE)

trained_model <- load_model("~/R/Data/model.hdf5")

# UMAP
trained_model <- run_umap(trained_model)
plot_dimred(trained_model, method = "UMAP", color_by = "celltype")


#---------------------------StabMap---------------------------------------------

# Seperation ATAC Multiome
atac <- rep("ATAC", 5016)
multiomes <- rep("Multiome", 5016)

names <- c(atac, multiomes)

# Assay Types
assayType = ifelse(rownames(logcounts_all) %in% rownames(sce.rna),
                   "rna", "atac")

# List for StabMap
assay_list = list(
  ATAC = logcounts_all[assayType %in% c("atac"), names %in% c("ATAC")],
  Multiome = logcounts_all[assayType %in% c("rna", "atac"), names %in% c("Multiome")]
)

# StabMap
mosaicDataUpSet(assay_list)
plot(mosaicDataTopology(assay_list))

stab = stabMap(assay_list,
               reference_list = c("Multiome"),
               plot = FALSE)

# UMAP
stab_umap = calculateUMAP(t(stab))

CellType <- setNames(metadata_ct, colnames(logcounts_all))

full_umap <- merge( as.data.frame(CellType), stab_umap, by =0 )

ggplot(full_umap) +
  geom_point(aes(x = V1, y = V2, color = CellType)) +
  theme_light()

#---------------------------Comparison------------------------------------------
plot_dimred(trained_model, method = "UMAP", color_by = "celltype", dot_size =1)
ggplot(full_umap) +
  geom_point(aes(x = V1, y = V2, color = CellType), size = .1) +
  theme_light()

#---------------------------Silhoutte-------------------------------------------

#StabMap

full_umap <- full_join(full_umap %>% group_by(CellType) %>% summarise(n = n()) %>% mutate(k = 1:14) %>% select(-n), full_umap)
full_umap_coords <- full_umap %>% select("V1", "V2")

sil_Stab <- silhouette(full_umap$k, dist(full_umap_coords))
sil_Stab_sum <- sil_Stab %>% 
  as.data.frame() %>% group_by(cluster) %>% summarise(score = mean(sil_width), 
                                                      frac_pos = sum(sil_width > 0)/n(),
                                                      pos_score = sum((sil_width>0)*sil_width)/sum(sil_width > 0))


#MOFA
MOFA_cluster <- full_join(as.data.frame(trained_model@dim_red), as.data.frame(metadata) %>% rownames_to_column("sample"), by = join_by(UMAP.sample == sample)) %>%
  full_join(.,y = full_umap %>% group_by(CellType) %>% summarise(n = n()) %>% mutate(k = 1:14) %>% select(-n), by = join_by(celltype  == CellType))

sil_MOFA <- silhouette(MOFA_cluster$k, dist(as.data.frame(trained_model@dim_red)))
sil_MOFA_sum <- sil_MOFA %>% 
  as.data.frame() %>% group_by(cluster) %>% summarise(score = mean(sil_width), 
                                                      frac_pos = sum(sil_width > 0)/n(),
                                                      pos_score = sum((sil_width>0)*sil_width)/sum(sil_width > 0))


#---------------------------Imputation------------------------------------------
imp = imputeEmbedding(
  assay_list,
  stab,
  reference = colnames(assay_list[["Multiome"]]),
  query = colnames(assay_list[["ATAC"]]))

trained_model <- impute(trained_model)
# imp$Multiome[953:963, 1:5]
# trained_model@imputed_data$ATAC[[1]][1:10, 1:5]

