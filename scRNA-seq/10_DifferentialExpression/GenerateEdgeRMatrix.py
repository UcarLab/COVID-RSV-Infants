import numpy as np
import pandas as pd
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Gets a raw count matrix from scRNA-seq based on barcodes')
parser.add_argument("sampleinfo")
parser.add_argument("barcodes")
parser.add_argument("outputfile")

args = parser.parse_args()
sampleinfofile = args.sampleinfo
barcodefile = args.barcodes 
outputfile = args.outputfile

samples = pd.read_csv(sampleinfofile, sep="\t", header=None).values
samplebarcodes = pd.read_csv(barcodefile, index_col=0)

def getOriginalBarcodes(barcodes, sample, samplevar='sample'):
    barcodes_with_sample = barcodes[barcodes[samplevar] == sample].index
    rv = []
    for curbarcode in barcodes_with_sample:
        rv.append(curbarcode.split("-"+sample)[0])
    return rv


alldata = []
samplenames = []
for cursample in samples:
    samplename = cursample[0]
    samplenames.append(samplename)

    curmat = cursample[1]
    curdata = sc.read_10x_mtx(curmat)
    curbarcodes = getOriginalBarcodes(samplebarcodes, samplename)
    curdata = curdata[curbarcodes].copy()

    alldata.append(curdata)

    
mergeddata = alldata[0].concatenate(alldata[1:], batch_categories=samplenames)

rmat = np.empty((len(mergeddata.var), len(samplenames)))
rmat[:] = np.nan

index=0
for cursample in samplenames:
    cursel = mergeddata.obs['batch']==cursample
    
    if sum(cursel) > 0:
        rmat[:,index] = mergeddata.X[cursel].sum(0)
    index += 1
    
rmat_df = pd.DataFrame(rmat, index=mergeddata.var.index, columns=samplenames)
rmat_df.to_csv(outputfile)

