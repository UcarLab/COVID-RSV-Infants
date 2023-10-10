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
parser.add_argument("tfidffile")
parser.add_argument("outfile")

args = parser.parse_args()

adatafile = args.anndatafile
tfidffile  = args.tfidffile
outfile = args.outfile

adata = sc.read(adatafile)
tfidfmat = scipy.sparse.load_npz(tfidffile)

print("Assigning matrix")
adata.X = tfidfmat

print("Writing file")
adata.write(outfile)