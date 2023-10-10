import numpy as np
import pandas as pd
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Script for merging adata objects.')
parser.add_argument("inputfile")
parser.add_argument("outputfile")

args = parser.parse_args()
inputfile = args.inputfile
outputfile = args.outputfile

#List of 1) sample ids, 2) adata objects
samples = pd.read_csv(inputfile, sep="\t", header=None).values

alldata = []
samplenames = []

for cursample in samples:
    samplename = cursample[0]
    adatapath = cursample[1]
    curdata = sc.read(adatapath)
    alldata.append(curdata)
    samplenames.append(samplename)

mergeddata = alldata[0].concatenate(alldata[1:], batch_categories=samplenames)
mergeddata.write(outputfile)