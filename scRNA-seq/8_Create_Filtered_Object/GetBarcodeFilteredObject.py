import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

parser = argparse.ArgumentParser(description='Subclustering!')
parser.add_argument("inputfile")
parser.add_argument("cellfile")
parser.add_argument("outputdir")

args = parser.parse_args()
inputfile = args.inputfile
cellfile = args.cellfile
outputdir = args.outputdir

##

sc.set_figure_params(dpi=200)
sc.settings.verbosity = 3
sc.settings.figdir = outputdir+"/"

##

adata = sc.read(inputfile)
adata = adata.raw.to_adata().copy()

barcodes = pd.read_csv(cellfile, header=None).values[:,0]

subset = adata[barcodes]

sc.pp.highly_variable_genes(subset, batch_key="Batch")
subset.raw = subset
subset = subset[:, subset.var.highly_variable]

sc.pp.scale(subset, max_value=10)
sc.tl.pca(subset, svd_solver='arpack',n_comps=100)
sce.pp.harmony_integrate(subset, 'Batch')
sc.pp.neighbors(subset, n_pcs=100, use_rep="X_pca_harmony")
sc.tl.umap(subset)
sc.tl.leiden(subset)

subset.write(outputdir+'/Subcluster.h5ad')

