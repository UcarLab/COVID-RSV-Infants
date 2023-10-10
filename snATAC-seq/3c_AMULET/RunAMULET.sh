#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --job-name=AMULET

module load singularity

if [ -n $SLURM_JOB_ID ] ; then
  SCRIPTDIR=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -f1 -d" ")
else
  SCRIPTDIR=$(realpath $0)
fi

SCRIPTDIR=$(dirname "$SCRIPTDIR")
echo "$SCRIPTDIR"

fragmentfile=$1
filteredcsv=$2
curout=$3


singularity exec SingleCellPython.sif /projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/AMULET/AMULET-v1.1/AMULET.sh "${fragmentfile}" "${filteredcsv}" /projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/AMULET/AMULET-v1.1/human_autosomes.txt /projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed ${curout} /projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/AMULET/AMULET-v1.1/
    
singularity exec SingleCellPython.sif python3 /projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/AMULET/AMULET-v1.1/AMULET.py --q 0.05 --rfilter /projects/ucar-lab/USERS/athib/Scripts/snATAC-seq/AMULET/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed  ${curout}/Overlaps.txt ${curout}/OverlapSummary.txt ${curout}


