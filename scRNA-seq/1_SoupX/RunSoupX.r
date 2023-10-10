library(SoupX)
library(Seurat)
library(DropletUtils)

args <- commandArgs(trailingOnly = TRUE)
infile=args[1]

sc = load10X(infile)
sc = autoEstCont(sc)
out = adjustCounts(sc)

outdir = paste0(infile,"/strainedCounts")
DropletUtils:::write10xCounts(outdir, out)

write.table(sc$soupProfile, paste0(outdir, "/soupProfile.txt"), sep = "\t")
