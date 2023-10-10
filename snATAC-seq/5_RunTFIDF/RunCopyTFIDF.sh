#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --mem=512G


module load singularity

if [ -n $SLURM_JOB_ID ] ; then
  SCRIPTDIR=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}' | cut -f1 -d" ")
else
  SCRIPTDIR=$(realpath $0)
fi

SCRIPTDIR=$(dirname "$SCRIPTDIR")
echo "$SCRIPTDIR"
singularity exec ${SCRIPTDIR}/SingleCellPython.sif python3 ${SCRIPTDIR}/CopyTFIDF_Step2.py "$@"



