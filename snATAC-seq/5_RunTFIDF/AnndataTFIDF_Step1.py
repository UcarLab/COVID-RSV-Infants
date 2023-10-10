import pandas as pd
import numpy as np
import scanpy as sc
import scanpy.external as sce
import argparse
from scipy.sparse import diags
from sklearn.preprocessing import binarize
import scipy

parser = argparse.ArgumentParser(description="Create adata object.")
parser.add_argument("anndatafile")
parser.add_argument("outputanndata")

args = parser.parse_args()

adata = args.anndatafile
outadata = args.outputanndata


def getTFIDF(binarizedcounts, scale = 1e4, transpose=False):
    #Reimplementation of Stuart & Butler (2018) TFIDF method (method=1)

    #transpose=True if Rows=barcodes, cols=peaks (Scanpy adata)
    #transpose=False if Rows=Peaks, cols=barcodes (Signac matrix)
    if transpose:
        binarizedcounts = binarizedcounts.T
        
    npeaks = np.sum(binarizedcounts, axis=0)
    tfdiag = diags(np.array(1.0 / np.array(npeaks)).flatten())
    tf = np.dot(binarizedcounts, tfdiag)
    idf = np.array(np.shape(binarizedcounts)[1]/np.sum(binarizedcounts, axis=1))[:,0]
    idfdiag = diags(idf)

    rv = np.dot(idfdiag, tf)
    rv = np.log1p(rv * scale)
    if transpose:
        return rv.T
    return rv


adata = sc.read(adata)

#Binarize to determine percentage of cells within each feature
binarizedcounts = binarize(adata.X, copy=True)

#Binarize, Run TDFIDF, Run PCA, Write output
tfidfmat = getTFIDF(binarizedcounts, transpose=True)
scipy.sparse.save_npz(outadata, tfidfmat)
