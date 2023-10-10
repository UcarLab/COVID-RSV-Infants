#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=72:00:00
#SBATCH --mem=120G
#SBATCH --job-name=CellRangerATAC

export PATH=/projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/cellranger-atac-2.1.0:$PATH

sampleid=$1 #The name that you want to provide for the run, this will create a directory in the output directory

fastqpath=$2 #Directory of the fastq files

sample=$3 #Sample id found in front of the fastq files (Useful if multiple samples are present in the fastq path)

referencepath=$4 #/projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

cd $5 #Cange path to the outputdirectory

    cellranger-atac count --id=${sampleid} \
                   --reference=${referencepath} \
                   --fastqs=${fastqpath} \
                   --sample=${sample} \
                   --localcores=8 \
                   --localmem=64