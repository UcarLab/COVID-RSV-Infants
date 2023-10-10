library(Seurat)
library(Azimuth)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

infile=args[1]
outfile=args[2]

mydata <- readRDS(infile)
mydata <- RunAzimuth(mydata, reference = "pbmcref")
saveRDS(mydata, outfile)

