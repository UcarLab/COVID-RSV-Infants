library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

samplename=args[1]
inputmatrixfiledirectory=args[2]
inputdoubletsummary=args[3]
outputdirectory=args[4]

###

matrixdata <- Read10X(data.dir = inputmatrixfiledirectory)
scrublet <- read.table(inputdoubletsummary, sep="\t", row.names=1, header=T)
colnames(scrublet) <- c("Doublet_score", "Is_doublet")
  
scrublet$Doublet_score = as.numeric(scrublet$Doublet_score)
  
obj <- CreateSeuratObject(counts = matrixdata, project = samplename, min.cells = 3, min.features = 200)
  
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]") 
obj <- AddMetaData(obj, scrublet)  


pdf(paste0(outputdirectory, "/", samplename, "_violin.pdf"), width=8, height=7)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_nfeature_zoomed.pdf"), width=8, height=7)
print(VlnPlot(obj, features = c("nFeature_RNA"), y.max = 1000))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_ncount_zoomed.pdf"), width=8, height=7)
print(VlnPlot(obj, features = c("nCount_RNA"), y.max = 5000))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_mt_zoomed.pdf"), width=8, height=7)
print(VlnPlot(obj, features = c("percent.mt"), y.max = 10))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_rp_zoomed.pdf"), width=8, height=7)
print(VlnPlot(obj, features = c("percent.rb"), y.max = 10))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_ncount_mt.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_nfeature_mt.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.mt")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_ncount_nfeature.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_ncount_rb.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.rb")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_nfeature_rb.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "percent.rb")+ggtitle(samplename))
dev.off()


pdf(paste0(outputdirectory, "/", samplename, "_rb_mt.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "percent.rb", feature2 = "percent.mt")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_ncount_doublet.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "Doublet_score")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_nfeature_doublet.pdf"), width=8, height=7)
print(FeatureScatter(obj, feature1 = "nFeature_RNA", feature2 = "Doublet_score")+ggtitle(samplename))
dev.off()

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000) 
all.genes <- rownames(obj)
obj <- ScaleData(obj, features=all.genes)
obj <- RunPCA(obj, features = VariableFeatures(object = obj))
obj <- FindNeighbors(obj, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.8)
obj <- RunUMAP(obj, dims=1:20, verbose=F)

pdf(paste0(outputdirectory, "/", samplename, "_umap_beforeQC.pdf"), width=8, height=7)
print(DimPlot(obj, label.size=4, repel=T, label = T)+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_umap_doublet.pdf"), width=8, height=7)
print(FeaturePlot(obj, features="Doublet_score")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_umap_mt.pdf"), width=8, height=7)
print(FeaturePlot(obj, features="percent.mt")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_umap_nCount.pdf"), width=8, height=7)
print(FeaturePlot(obj, features="nCount_RNA")+ggtitle(samplename))
dev.off()

pdf(paste0(outputdirectory, "/", samplename, "_umap_nFeature.pdf"), width=8, height=7)
print(FeaturePlot(obj, features="nFeature_RNA")+ggtitle(samplename))
dev.off()



