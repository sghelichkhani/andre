#!/bin/bash
#SBATCH --output=/home/horbach/slurm-%j.out
#SBATCH --workdir=/DATA/
#SBATCH --job-name=terra-adj
#SBATCH --get-user-env
#SBATCH --partition=BIG
#SBATCH --ntasks=512
#SBATCH --ntasks-per-node=24
#SBATCH --mail-type=all
#SBATCH --export=ALL
#SBATCH --time=60-00:00:00

CASENUM=357
PROC=0512
MT=512


# create folder for final results
 mkdir -p /SCRATCH/${USER}/$CASENUM
 mkdir -p /SCRATCH/${USER}/$CASENUM/c-files
 mkdir -p /SCRATCH/${USER}/$CASENUM/gmt

# distribute data to compute nodes
 pdsh -R ssh -f ${SLURM_JOB_NUM_NODES} -w ${SLURM_JOB_NODELIST} "mkdir -p /DATA/${USER}/$CASENUM && cp -ur /SCRATCH/${USER}/platemaps_mt$MT/ /DATA/${USER}/ ; mkdir -p /DATA/${USER}/$CASENUM/output; cp /SCRATCH/${USER}/TERRA$PROC''_mt$MT/*ter* /DATA/${USER}/$CASENUM/ ; "

# run compute job
cd /DATA/${USER}/$CASENUM
mpirun.openmpi terra

# remove data from compute nodes
pdsh -R ssh -f ${SLURM_JOB_NUM_NODES} -w ${SLURM_JOB_NODELIST} "rm -rf /DATA/${USER}/$CASENUM "
