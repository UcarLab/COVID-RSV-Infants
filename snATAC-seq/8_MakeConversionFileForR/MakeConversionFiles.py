import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io
import argparse

parser = argparse.ArgumentParser(description="Create adata object.")
parser.add_argument("inputfile")
parser.add_argument("subsetvar") #Var used for subsetting
parser.add_argument("clusters")
parser.add_argument("annotationvar")
parser.add_argument("pcavar")
parser.add_argument("outputprefix")

args = parser.parse_args()

inputfile = args.inputfile
subsetvar = args.subsetvar
annotationvar = args.annotationvar
pcavar = args.pcavar
clusters = args.clusters.split(",")
outputprefix = args.outputprefix

adata = sc.read(inputfile)

try:
    subset = adata[adata.obs[subsetvar].isin(clusters)].raw.to_adata().copy()
except:
    subset = adata[adata.obs[subsetvar].isin(clusters)].copy()

#sc.pp.highly_variable_genes(subset, batch_key='Batch')
#subset.raw = subset
#subset = subset[:, subset.var.highly_variable]

io.mmwrite(outputprefix+".mtx", subset.X.T)

pd.DataFrame(subset.obsm[pcavar]).to_csv(outputprefix+"_pca.csv", index=None, header=None)
pd.DataFrame(subset.obs.index).to_csv(outputprefix+"_ObsIndex.csv", index=None, header=None)
pd.DataFrame(subset.var.index).to_csv(outputprefix+"_VarIndex.csv", index=None, header=None)
pd.DataFrame(subset.obs[annotationvar]).to_csv(outputprefix+"_Annotation.csv", header=None)

#pd.DataFrame(subset.var.index[subset.var.highly_variable]).to_csv(outputprefix+"_HighlyVariable.csv", index=None, header=None)
