import numpy as np
import pandas as pd
import argparse
import peakoverlap as po


parser = argparse.ArgumentParser(description="Create summary and csv file with multiplets removes")
parser.add_argument("filelist")
parser.add_argument("outfile")

args = parser.parse_args()

infile = args.filelist
outfile = args.outfile

#read input file
inputdf = pd.read_csv(infile, sep="\t", header=None).values

cellranger = inputdf[:, 0]

allpeaks = []
for i in range(len(cellranger)):
    curpeaks = pd.read_csv(cellranger[i], sep="\t", header=None, comment="#")
    allpeaks.append(curpeaks.values[:,0:3])
    
union = po.getUnionPeaks(allpeaks)
peaksizes = union[:,2]-union[:,1]
union_sizefiltered = union[(peaksizes > 20) & (peaksizes < 10000),:]
union_noy = union_sizefiltered[(union_sizefiltered[:,0] != "chrY"),:]

chrunion = []
for i in range(len(union_noy)):
    if(union_noy[i][0].startswith("chr")):
        chrunion.append(union_noy[i])
        
union_chrfiltered = np.array(chrunion)

ocount,_ = po.getOverlapCount(union_chrfiltered,allpeaks)

union_final = union_chrfiltered[ocount > 3,:]

pd.DataFrame(union_final).to_csv(outfile, sep="\t", index=None, header=None)