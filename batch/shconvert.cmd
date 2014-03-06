#!/bin/bash
#SBATCH --output=/home/horbach/slurm-%j.out
#SBATCH --workdir=/DATA/
#SBATCH --job-name=SHconvert
#SBATCH --get-user-env
#SBATCH --partition=BIG
#SBATCH --ntasks=512
#SBATCH --mail-type=all
#SBATCH --export=ALL
#SBATCH --time=0-08:00:00

CASENUM=373

# create folder for final results
# mkdir -p /SCRATCH/${USER}/$CASENUM
# mkdir -p /SCRATCH/${USER}/$CASENUM/c-files

# distribute data to compute nodes
# pdsh -R ssh -f ${SLURM_JOB_NUM_NODES} -w ${SLURM_JOB_NODELIST} "mkdir -p /DATA/${USER}/$CASENUM && cp -ur /SCRATCH/${USER}/platemaps_mt512/ /DATA/${USER}/ ; cp -r /SCRATCH/${USER}/TERRA1024_mt512/* /DATA/${USER}/$CASENUM/ ; "

# run compute job
cd /SCRATCH/${USER}/conversion
mpirun.openmpi SHconvert

# remove data from compute nodes
#pdsh -R ssh -f ${SLURM_JOB_NUM_NODES} -w ${SLURM_JOB_NODELIST} "rm -rf /DATA/${USER}/$CASENUM "
