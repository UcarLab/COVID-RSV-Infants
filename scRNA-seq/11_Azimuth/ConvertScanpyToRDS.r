library(Seurat)
library(patchwork)

##
infile = "/Users/athib/Desktop/CovidRSV/scRNA2/COVIDRSV_Pass2_Cleaned.h5ad"
h5seuratfile = "/Users/athib/Desktop/CovidRSV/scRNA2/COVIDRSV_Pass2_Cleaned.h5seurat"

outfile = "/Users/athib/Desktop/CovidRSV/scRNA2/COVIDRSV_Pass2_Cleaned.rds"

remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
Convert(infile, dest = "h5seurat", overwrite = TRUE)

mydata <- LoadH5Seurat(h5seuratfile, meta.data = FALSE, misc = FALSE)

DimPlot(mydata, label = TRUE, label.size = 3) + NoLegend()

saveRDS(mydata, outfile)


