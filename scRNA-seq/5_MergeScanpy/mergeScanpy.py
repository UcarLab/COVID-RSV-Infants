import numpy as np
import pandas as pd
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Script for merging samples in scanpy.')
parser.add_argument("inputfile")
parser.add_argument("outputfile")

args = parser.parse_args()
inputfile = args.inputfile
outputfile = args.outputfile

#TODO read tab delimited file with 1) samplename 2) matrix folder, 3) barcode file
samples = pd.read_csv(inputfile, sep="\t", header=None).values

alldata = []
samplenames = []

for cursample in samples:
    samplename = cursample[0]
    samplenames.append(samplename)

    curmat = cursample[1]
    curbarcodes = pd.read_csv(cursample[2], header=None).values[:,0]

    curdata = sc.read_10x_mtx(curmat)

    curdata = curdata[curbarcodes].copy()

    #Log normalize the data sample by sample
    sc.pp.normalize_total(curdata, target_sum=1e4)
    sc.pp.log1p(curdata)
    
    curdata.raw = curdata #Save the raw data

    alldata.append(curdata)

mergeddata = alldata[0].concatenate(alldata[1:], batch_categories=samplenames)
mergeddata.write(outputfile)
