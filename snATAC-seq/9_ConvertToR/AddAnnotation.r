library(future)
library(Signac)
library(Seurat)

options(future.globals.maxSize = 2000 * 1024^2)

args = commandArgs(trailingOnly = TRUE)
pbmc_rna_f = args[1]
annotationfile = args[2]

pbmc_rna <- readRDS(pbmc_rna_f)

data_annotation <- read.csv(annotationfile, header = F, row.names = 1)
pbmc_rna <- AddMetaData(object = pbmc_rna, metadata = data_annotation, col.name = 'Annotation')

saveRDS(pbmc_rna, pbmc_rna_f)