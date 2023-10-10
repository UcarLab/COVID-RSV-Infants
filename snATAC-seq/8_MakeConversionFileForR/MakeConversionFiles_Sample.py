import numpy as np
import pandas as pd
import scanpy as sc
from scipy import io
import argparse

parser = argparse.ArgumentParser(description="Create adata object.")
parser.add_argument("inputfile")
parser.add_argument("outdir")

args = parser.parse_args()

inputfile = args.inputfile
outdir = args.outdir

adata = sc.read(inputfile)

for cursample in adata.obs['sample'].cat.categories:
    subset = adata[adata.obs['sample'].isin([cursample])].raw.to_adata().copy()

    io.mmwrite(outdir+"/"+cursample+".mtx", subset.X.T)
    pd.DataFrame(subset.obs.index).to_csv(outdir+"/"+cursample+"_ObsIndex.csv", index=None, header=None)
    pd.DataFrame(subset.var.index).to_csv(outdir+"/"+cursample+"_VarIndex.csv", index=None, header=None)

