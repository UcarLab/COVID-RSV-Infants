library(Seurat)

args <- commandArgs(trailingOnly = TRUE)

infile='/Users/athib/Desktop/CovidRSV/scRNA2/COVIDRSV_Pass2_Cleaned_Azimuth.rds'

azimuthannotated <- readRDS(infile)

pdf("/Users/athib/Desktop/CovidRSV/scRNA2/AzimuthUmap_l1.pdf", width = 20, height = 10)
Seurat::DimPlot(azimuthannotated, group.by = 'predicted.celltype.l1')
dev.off()

pdf("/Users/athib/Desktop/CovidRSV/scRNA2/AzimuthUmap_l2.pdf", width = 20, height = 10)
Seurat::DimPlot(azimuthannotated, group.by = 'predicted.celltype.l2')
dev.off()

pdf("/Users/athib/Desktop/CovidRSV/scRNA2/AzimuthUmap_l3.pdf", width = 20, height = 10)
Seurat::DimPlot(azimuthannotated, group.by = 'predicted.celltype.l3')
dev.off()

write.csv(azimuthannotated$predicted.celltype.l2, "/Users/athib/Desktop/CovidRSV/scRNA2/AzimuthPredictions.csv")
