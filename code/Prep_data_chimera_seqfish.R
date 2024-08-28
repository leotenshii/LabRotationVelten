# This script is adapted from 
# https://github.com/MarioniLab/StabMap2021/blob/main/example/MGD_chimera_generate.Rmd and
# https://github.com/MarioniLab/StabMap2021/blob/main/example/MGD_seqFISH_generate.Rmd
# written by Shila Ghazanfar
#---------------------------Seed------------------------------------------------
set.seed(42)
#---------------------------Libraries-------------------------------------------
library(scater)
#-----------------------------scRNA-seq Data------------------------------------
# Three dataset will be combined:
# E8.5 mouse gastrulation atlas, WT/WT chimera1 and T−/−/WT chimera scRNA-seq data

# Data can be loaded from the MouseGastrulationData Bioconductor package. I wasn't able to install it from the server, therefore loading it from RDS
# Loading it from MouseGastrulationData:

# library(MouseGastrulationData)
# stages <- c("E8.5")
# mt <- MouseGastrulationData::AtlasSampleMetadata
# samples <- mt[mt[, "stage"] %in% stages, "sample"]
# atlas_raw <- EmbryoAtlasData(type = "processed", samples = samples)
# wt_samples <- c(5,6,7,8,9,10)
# t_samples <- c(1,2,5,6,7,8,9,10)
# 
# wt_chim_raw <- MouseGastrulationData::WTChimeraData(
#   type = "processed",
#   samples = wt_samples
# )
# t_chim_raw <- MouseGastrulationData::TChimeraData(
#   type = "processed",
#   samples = t_samples
# )

# Load from RDS
atlas_raw = readRDS("~/R/Data/atlas.RDS")
atlas_raw = logNormCounts(atlas_raw)

wt_chim_raw <- readRDS("~/R/Data/wt_chim.RDS")
wt_chim_raw <- logNormCounts(wt_chim_raw)

t_chim_raw <- readRDS("~/R/Data/t_chim.RDS")
t_chim_raw <- logNormCounts(t_chim_raw)


# Generate a joint SCE with all three datasets
scrna_cells <- c(paste0("atlas_", colnames(atlas_raw)),
                paste0("wt_chim_", colnames(wt_chim_raw)),
                paste0("t_chim_", colnames(t_chim_raw)))

scrna_sample <- c(paste0("atlas_", atlas_raw$sample),
                 paste0("wt_chim_", wt_chim_raw$sample),
                 paste0("t_chim_", t_chim_raw$sample))

scrna_celltype <- c(atlas_raw$celltype,
                   wt_chim_raw$celltype.mapped,
                   t_chim_raw$celltype.mapped)

scrna_experiment <- c(rep("atlas", ncol(atlas_raw)),
                     rep("wt_chim", ncol(wt_chim_raw)),
                     rep("t_chim", ncol(t_chim_raw)))

scrna_tomato <- c(rep(NA, ncol(atlas_raw)),
                 wt_chim_raw$tomato,
                 t_chim_raw$tomato)

scrna_cData <- data.frame(
  cell = scrna_cells,
  sample = scrna_sample,
  celltype = scrna_celltype,
  experiment = scrna_experiment,
  tomato = scrna_tomato)

scrna_exprs <- cbind(
  rbind(logcounts(atlas_raw), "tomato-td" = 0),
  logcounts(wt_chim_raw),
  logcounts(t_chim_raw))

colnames(scrna_exprs) <- scrna_cData$cell

scrna <- SingleCellExperiment(assays = list(logcounts = scrna_exprs),
                             colData = scrna_cData,
                             rowData = rowData(wt_chim_raw))

# save the object as RDS
saveRDS(scrna, file = "~/R/Data/MGD_chimera_scrna.Rds")

#-----------------------------seqFISH Data--------------------------------------
# Data can be downloaded from https://content.cruk.cam.ac.uk/jmlab/SpatialMouseAtlas2020/
meta <- readRDS("~/R/Data/metadata.Rds")
exprs <- readRDS("~/R/Data/exprs.Rds")

sce <- SingleCellExperiment(assays = list(logcounts = exprs),
                           colData = meta)
sce <- sce[, sce$embryo == "embryo1"]

saveRDS(sce, file = "~/R/Data/seqFISH.Rds")