import numpy as np
import pandas as pd
import tabix
from scipy.sparse import csr_matrix
from scipy.sparse import diags
import anndata as ad
import argparse

#Arguments:
#1. fragmentfile
#2. Single cell CSV File
#3. Peaks
#4. output file

parser = argparse.ArgumentParser(description="Create adata object.")
parser.add_argument("fragmentfile")
parser.add_argument("singlecellcsv")
parser.add_argument("peakfile")
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
peakfile = args.peakfile
outputfile = args.outputfile

barcodecol = args.barcodecol
iscellcol = args.iscellcol

#Functions
def getPeakDictionary(peakfile, chrfile=None):
    peaks = pd.read_csv(peakfile, header=None, comment="#", sep="\t").values
    if chrfile:
        chromosomes = list(pd.read_csv(chrfile, header=None).values)
    else:
        chromosomes = list(np.unique(peaks[:,0]))
        
    peaklabels = []
    peakdict = dict()
    peakidx = 0
    for curpeak in peaks:
        chrom = curpeak[0]
        if chrom in chromosomes:
            start = int(curpeak[1])
            end = int(curpeak[2])
            peaklabel = chrom+"_"+str(start)+"_"+str(end)
            peaklabels.append(peaklabel)
            peakdict[peaklabel] = [peakidx, chrom, start, end]
            peakidx += 1  
    return peaklabels, peakdict

def getBarcodeDictionary(barcodefile, barcodecol='barcode', iscellcol='is__cell_barcode', sep=","):
    barcodes = pd.read_csv(barcodefile)
    selectedbarcodes = barcodes[barcodecol].values
    if iscellcol:
        selectedbarcodes = selectedbarcodes[barcodes[iscellcol] == 1]
    
    barcodelabels = []
    barcodedict = dict()
    bcidx = 0
    for curbarcode in selectedbarcodes:
        barcodelabels.append(curbarcode)
        barcodedict[curbarcode] = bcidx
        bcidx += 1
        
    return barcodelabels,barcodedict

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

def addMetaData(adata, singlecellcsv, barcodecol='barcode', iscellcol='is__cell_barcode'):
    singlecelldf = pd.read_csv(singlecellcsv)
    selecteddata = singlecelldf[singlecelldf[iscellcol] == 1]
    for curcol in singlecelldf.columns:
        if curcol != barcodecol and curcol != iscellcol:
            columndata =  pd.Series(list(selecteddata[curcol].values), dtype=int)
            columndata.index = selecteddata[barcodecol]
            adata.obs[curcol] = columndata

#Running the code
peaklabels,peakdict = getPeakDictionary(peakfile)
barcodelabels,barcodedict = getBarcodeDictionary(singlecellcsv, barcodecol=barcodecol, iscellcol=iscellcol)
tb = tabix.open(fragmentfile)

rows = []
cols = []
data = []

for curpeak in peaklabels:
    peakinfo = peakdict[curpeak]
    colidx = peakinfo[0]
    curchr = peakinfo[1]
    curstart = peakinfo[2]
    curend = peakinfo[3]

    try:
        colcounts = getCountsByBarcode(tb.query(curchr, curstart, curend), barcodedict)

        for curcount in colcounts:
            cols.append(colidx)
            rows.append(curcount[0])
            data.append(curcount[1])
    except:
        pass
        
counts = csr_matrix((data, (rows, cols)), shape=(len(barcodelabels), len(peaklabels)))

adata = ad.AnnData(counts, dtype=int)
adata.obs_names = barcodelabels
adata.var_names = peaklabels
addMetaData(adata, singlecellcsv, barcodecol=barcodecol, iscellcol=iscellcol)
adata.write(outputfile)
