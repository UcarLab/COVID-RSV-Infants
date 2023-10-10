#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=128G
#SBATCH --job-name=CellRanger6_1_2

export PATH=/PathTo/CellRanger_v6.1/cellranger-6.1.2:$PATH

sampleid=$1 #The name that you want to provide for the run, this will create a directory in the output directory
fastqpath=$2 #Directory of the fastq files
sample=$3 #Sample id found in front of the fastq files (Useful if multiple samples are present in the fastq path)
referencepath=$4 #PathTo/CellRanger_v6.1/References/refdata-gex-GRCh38-2020-A
cd $5 #Cange path to the outputdirectory

cellranger count --id=${sampleid} --transcriptome=${referencepath} --fastqs=${fastqpath} --sample=${sample}

