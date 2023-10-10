library(Seurat)
library(Signac)
library(Matrix)



args = commandArgs(trailingOnly = TRUE)
anndata_x_f = args[1]
anndata_vars_f = args[2]
anndata_obs_f = args[3]
anndata_pca_f = args[4]
outfile = args[5]

print("Read Matrix")
annadata_x <- Matrix::readMM(anndata_x_f)

print("Read Vars")
anndata_vars <- read.table(anndata_vars_f)

print("Read Obs")
anndata_obs <- read.table(anndata_obs_f)

print("Read PCA")
anndata_pca <- read.csv(anndata_pca_f,header = F)

print("Add row/column data")
rownames(annadata_x) = anndata_vars$V1
colnames(annadata_x) = anndata_obs$V1

print("Convert Matrix")
annadata_x_reformatted <- as(object = annadata_x, Class = 'CsparseMatrix')

print("Create assay object")
assayobj <- CreateAssayObject(
  data = annadata_x_reformatted,
  min.cells = 0,
  min.features = 0,
  row.names = anndata_vars
)



seurat_obj <- CreateSeuratObject(
  assayobj,
  project = "SeuratObject",
  assay = "GENEACTIVITY",
)


cell.embeddings <- as.matrix(anndata_pca)
rownames(cell.embeddings) <- colnames(seurat_obj)
colnames(cell.embeddings) <- paste0("PCA_", seq(1,ncol(anndata_pca)))

reduction.data <- CreateDimReducObject(
  embeddings = cell.embeddings,
  assay = 'GENEACTIVITY',
  key = "PCA_",
  global = TRUE,
)

seurat_obj[['PCA_']] <- reduction.data
saveRDS(seurat_obj, outfile)

