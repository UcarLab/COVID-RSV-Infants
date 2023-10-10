import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import matplotlib.pyplot as plt
import anndata as ad
import argparse
import peakoverlap
from sklearn.preprocessing import binarize

parser = argparse.ArgumentParser(description="Create new anndata object with harmonized matrix and reduced gene activity scores.")
parser.add_argument("peakanndata")
parser.add_argument("geneactivityanndata")
parser.add_argument("npcs")
parser.add_argument("exclusionlist")
parser.add_argument("outfile")

#Parse arguments
args = parser.parse_args()
peak_adata_file = args.peakanndata
ga_adata_file = args.geneactivityanndata
npcs = int(args.npcs)
exclusionlist = args.exclusionlist
outfile = args.outfile

#Read anndata objects
peak_adata = sc.read(peak_adata_file)
ga_adata = sc.read(ga_adata_file)

#filter peaks that overlap exclusion list regions
unionpeaks = []
for curpeak in peak_adata.var.index:
    splitpeak = curpeak.split("_")
    unionpeaks.append([splitpeak[0], int(splitpeak[1]), int(splitpeak[2])])
unionpeaks = np.array(unionpeaks, dtype=np.object)
exclusionlist_sites = pd.read_csv(exclusionlist, sep='\t', header=None)

counts,_ =  peakoverlap.getOverlapCount(unionpeaks, [exclusionlist_sites.values])
exclusionselection = (counts == 0)
peak_adata = peak_adata[:, exclusionselection]
​

#Concatenate TFIDF normalized peaks with normalized gene activities
combined_adata = peak_adata.copy()
combined_adata.X = ​combined_adata.X.astype('float32') 
combined_adata.obs['passed_filters'] = peak_adata.obs['passed_filters']
combined_adata.obs['mitochondrial'] = peak_adata.obs['mitochondrial']
combined_adata.obs['batch'] = peak_adata.obs['batch']

#Compute LSA with the TFIDF matrix (TFIDF+SVD)
sc.tl.pca(combined_adata, zero_center=False, n_comps=npcs) #Zero_center:False turns this into truncated SVD

#Writing the combined data file to look at PCA loadings
combined_adata.write(outfile+"_all.h5ad")

#Add the harmony matrix to the reduced gene activity matrix
ga_adata.obsm['X_pca'] = combined_adata.obsm['X_pca']
ga_adata.obs['Read Depth'] = np.log1p(combined_adata.obs['passed_filters'])
ga_adata.obs['Mitochondrial Reads'] = np.log1p(combined_adata.obs['mitochondrial'])
​ga_adata.X = ​ga_adata.X.astype('float32') 

#Write the output. X_pca_harmony should be sufficient for neighbors, umap, pass 1 subclustering
ga_adata.write(outfile)

