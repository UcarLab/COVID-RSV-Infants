library(future)
library(Signac)
library(Seurat)

options(future.globals.maxSize = 2000 * 1024^2)

args = commandArgs(trailingOnly = TRUE)
pbmc_rna_f = args[1]
pbmc_atac_f = args[2]
prefix = args[3]

pbmc_rna <- readRDS(pbmc_rna_f)
pbmc_atac <- readRDS(pbmc_atac_f)

features_var = rownames(pbmc_rna)
if (length(args) > 3){
  features_var = read.csv(args[4], header=F)$V1
}

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc_atac,
  reduction = 'cca',
  features = features_var
)

saveRDS(transfer.anchors, paste0(prefix,"_transferanchors.rds"))


predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$Annotation,
  weight.reduction = pbmc_atac[['PCA_Harmony_']],
  dims = 1:101
)

saveRDS(predicted.labels, paste0(prefix,"_predictedlabels.rds"))
