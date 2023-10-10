import numpy as np
import pandas as pd
import tabix
from scipy.sparse import csr_matrix
from scipy.sparse import diags
import anndata as ad
import argparse
import scanpy as sc

#Arguments:
#1. fragmentfile
#2. Single cell CSV File
#3. Peaks
#4. output file

parser = argparse.ArgumentParser(description="Create adata object.")
parser.add_argument("fragmentfile")
parser.add_argument("singlecellcsv")
parser.add_argument("genefile")
parser.add_argument("outputfile")

#Optional arguments
#1. Barcode column (default="barcode")
#2. Is cell (default="is__cell_barcode")

parser.add_argument('--bc', dest='barcodecol', type=str, default="barcode",
                    help='Name of the barcode column in the single cell csv file. (Default: barcode)')
parser.add_argument('--iscell', dest='iscellcol', type=str, default="is__cell_barcode",
                    help='Column specifying which barcodes are cells (value 1 = is a cell). (Default: is__cell_barcode)')

args = parser.parse_args()

fragmentfile = args.fragmentfile
singlecellcsv = args.singlecellcsv
genefile = args.genefile
outputfile = args.outputfile

barcodecol = args.barcodecol
iscellcol = args.iscellcol

#Functions
def getGeneActivityDictionary(genefile):
    genepositions = pd.read_csv(genefile, header=None, comment="#", sep="\t").values
    
    genelabels = []
    rv = dict()
    geneidx = 0
    for curgene in genepositions:
        genename = curgene[0]
        curchr = curgene[1]
        curstart = curgene[2]
        curend = curgene[3]
        rv[genename] = [geneidx, curchr, curstart, curend]
        genelabels.append(genename)
        geneidx += 1

    return genelabels, rv

def getBarcodeDictionary(barcodefile, barcodecol='barcode',validreadcol='passed_filters', iscellcol='is__cell_barcode', sep=","):
    barcodes = pd.read_csv(barcodefile)
    
    selectedbarcodes = barcodes[barcodecol].values
    if iscellcol:
        selectedbarcodes = selectedbarcodes[barcodes[iscellcol] == 1]
    
    validreads=barcodes[validreadcol][barcodes[iscellcol] == 1]
    barcodelabels = []
    barcodedict = dict()
    bcidx = 0
    for curbarcode in selectedbarcodes:
        barcodelabels.append(curbarcode)
        barcodedict[curbarcode] = bcidx
        bcidx += 1
        
    return barcodelabels,barcodedict,validreads

def getCountsByBarcode(records, barcodes):
    """
    @records -  The record iterator that provides a list of chr, start, end, barcode for each entry
    @barcode - Dictionary of column indices. 
        
    @return - List of [colindex, count]
    """
    countdict = dict()
    for currecord in records:
        curbc = currecord[3]
        if curbc in barcodes:
            colidx = barcodes[curbc]
            if colidx not in countdict:
                countdict[colidx] = 0
            countdict[colidx] = countdict[colidx]+1
            
    colsortedcounts = sorted(countdict.keys())
    rv = []
    for curcol in colsortedcounts:
        rv.append([curcol, countdict[curcol]])
    return rv


#Running the code
genelabels, geneactivitydict = getGeneActivityDictionary(genefile)

barcodelabels,barcodedict, validreads = getBarcodeDictionary(singlecellcsv, barcodecol=barcodecol, iscellcol=iscellcol)
tb = tabix.open(fragmentfile)

rows = []
cols = []
data = []


for curgene in geneactivitydict:
    curgeneinfo = geneactivitydict[curgene]
    colidx = curgeneinfo[0]
    curchr = curgeneinfo[1]
    curstart = curgeneinfo[2]
    curend = curgeneinfo[3]

    try:
        colcounts = getCountsByBarcode(tb.query(curchr, curstart, curend), barcodedict)

        for curcount in colcounts:
            cols.append(colidx)
            rows.append(curcount[0])
            data.append(curcount[1])
    except Exception as e:
        print(e)
        print(curgene)
        print(curchr)
        print(curstart)
        print(curend)
        pass
        
counts = csr_matrix((data, (rows, cols)), shape=(len(barcodelabels), len(genelabels)))

adata = ad.AnnData(counts, dtype='float32')
adata.obs_names = barcodelabels
adata.var_names = genelabels

scalefactor = np.median(validreads)
rowsums = np.sum(adata.X, axis=1)

normX = csr_matrix(adata.X / rowsums * scalefactor)
adata.X = normX
sc.pp.log1p(adata)

print(adata)
adata.write(outputfile)
