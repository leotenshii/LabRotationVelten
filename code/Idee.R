library(MOFA2)
library(SingleCellMultiModal)
library(scater)
library(scran)
library(StabMap)

# Peripheral Blood Mononuclear Cells provided by 10x Genomics website
# 10x Genomics Multiome technology enables simultaneous profiling of the transcriptome 
# (using 3’ gene expression) and epigenome (using ATAC-seq) from single cells to deepen 
# our understanding of how genes are expressed and regulated across different cell types.
mae <- scMultiome("pbmc_10x", mode = "*", dry.run = FALSE, format = "MTX")


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

#----------------------------MOFA-----------------------------------------------


log_counts100 <-logcounts_all[, 1:100]

log_counts100NA <- log_counts100

log_counts100NA[953:1740, 1:50] <- NA

mofa_list <- list(
  RNA = as.matrix( log_counts100NA[1:952,]),
  ATAC = as.matrix(log_counts100NA[953: 1740,])
)


model <- create_mofa(mofa_list)
plot_data_overview(model)

MOFAobject <- prepare_mofa(model)


run_mofa(MOFAobject, outfile = "Z:/Leoni/model.hdf5", use_basilisk = TRUE)

trained_model <- load_model("Z:/Leoni/model.hdf5")


trained_model <- run_umap(trained_model)

#color <- c(rep("RNA", 50), rep("Multiome", 50))

plot_dimred(trained_model, method = "UMAP")


#----------------------------StabMap--------------------------------------------
dataType = setNames(sample(c("RNA", "Multiome"), ncol(log_counts100),
                           prob = c(0.5,0.5), replace = TRUE),
                    colnames(log_counts100))



assay_list = list(
  RNA = log_counts100[assayType %in% c("rna"), random_names %in% c("RNA")],
  Multiome = log_counts100[assayType %in% c("rna", "atac"), random_names %in% c("Multiome")]
)

mosaicDataUpSet(assay_list)
plot(mosaicDataTopology(assay_list))

stab = stabMap(assay_list,
               reference_list = c("Multiome"),
               plot = FALSE)

stab_umap = calculateUMAP(t(stab))
plot(stab_umap, pch = 16, cex = 0.3, col = factor(dataType[rownames(stab)]))



#---------------------------Comparison------------------------------------------
plot_dimred(trained_model, method = "UMAP")
plot(stab_umap, pch = 16, cex = 0.3, col = factor(dataType[rownames(stab)]))

#---------------------------Imputation------------------------------------------
imp = imputeEmbedding(
  assay_list,
  stab,
  reference = colnames(assay_list[["Multiome"]]),
  query = colnames(assay_list[["RNA"]]))

trained_model <- impute(trained_model)
# trained_model@data$ATAC[[1]][1:10, 1:5]
# trained_model@imputed_data$ATAC[[1]][1:10, 1:5]

#---------------
rnas <- rep("RNA", 50)
multiomes <- rep("Multiome", 50)

# Combine the names
random_names <- c(rnas, multiomes)

# Shuffle the names

random_names <- sample(random_names)
