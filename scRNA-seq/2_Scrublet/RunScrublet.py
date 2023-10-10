import scrublet
import scrublet as scr
import scipy.io
import numpy as np
import os
import glob
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Scrublet Help Script.')
parser.add_argument("matrix")
parser.add_argument("barcodes")
parser.add_argument("outputdirectory")

parser.add_argument('--dr', dest='doubletrate', type=float, default=0.06,
                    help='Doublet rate (Default: 0.06)')

parser.add_argument('--mincounts', dest='mincounts', type=int, default=2,
                    help='Minimum counts (Default: 2)')

parser.add_argument('--mincells', dest='mincells', type=int, default=3,
                    help='Minimum cells (Default: 3)')

parser.add_argument('--mingenevar', dest='mingenevar', type=int, default=85,
                    help='Min Gene Variability PCTL (Default: 85)')

parser.add_argument('--npcs', dest='npcs', type=int, default=30,
                    help='Number of prinicpal components (Default: 30)')

args = parser.parse_args()

inputmatrixfile = args.matrix
inputbarcodesfile = args.barcodes
outdir = args.outputdirectory
doubletrate = args.doubletrate
mincounts = args.mincounts
mincells = args.mincells
mingenevar = args.mingenevar
npcs = args.npcs

parameters = [["Input Matrix:", inputmatrixfile],
              ["Input Barcodes:", inputbarcodesfile],
              ["Output Directory:", outdir],
              ["Doublet Rate:", doubletrate],
              ["Min Counts:", mincounts],
              ["Min Cells:", mincells],
              ["Min Gene Variability PCTL:", mingenevar],
              ["Number of PCs:", npcs]]

pd.DataFrame(np.array(parameters)).to_csv(outdir+"/RunParameters.txt", sep="\t", index=None, header=None)

cm = scipy.io.mmread(inputmatrixfile).T.tocsc()
barcodes = pd.read_csv(inputbarcodesfile, header=None).values[:,0]

scrub = scr.Scrublet(cm, expected_doublet_rate=doubletrate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=mincounts,min_cells=mincells,min_gene_variability_pctl=mingenevar, n_prin_comps=npcs)

output = np.concatenate((barcodes.reshape(-1,1), doublet_scores.reshape(-1,1), predicted_doublets.reshape(-1,1)), axis=1)

pd.DataFrame(output, columns=["Cell Barcode", "Doublet Score", "Default Doublet Prediction"]).to_csv(outdir+"/DoubletSummary.txt", sep="\t", index=None)

fig_hist,_ = scrub.plot_histogram();
fig_hist.savefig(outdir+"/DoubletHist.pdf", bbox_inches='tight')

scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
fig_UMAP,_ = scrub.plot_embedding('UMAP', order_points=True);
fig_UMAP.savefig(outdir+"/DoubletUMAP.pdf", bbox_inches='tight')

scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9))
fig_tSNE,_ = scrub.plot_embedding('tSNE', order_points=True);
fig_tSNE.savefig(outdir+"/DoubletTSNE.pdf", bbox_inches='tight')
