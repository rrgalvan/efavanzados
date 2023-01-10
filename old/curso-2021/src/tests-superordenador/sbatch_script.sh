#!/bin/bash

#------------------------------------
# Global variables
MPI_PROCESSES=1 #32
export PROGRAM=diffusion-simple.py

## no funciona >>>>>>>>>
export NODES=4
#export PROGRAM=advect_diffus.py
export JOB_NAME="adv_dif"
export OUTPUT_DIR="./out"
mkdir -p $OUTPUT_DIR

#SBATCH --job-name=$JOBNAME
#SBATCH --partition=cn
#SBATCH --nodes=$NODES
#SBATCH --error=$OUTPUT_DIR/job.%J.err
#SBATCH --output=$OUTPUT_DIR/job.%J.out
#SBATCH --mail-user=rafael.rodriguez@uca.es
#SBATCH --mail-type=ALL
#SBATCH --requeue
# <<<<<<<<<<<<<<<<<<<<<<

#---------------------------------------------------------------------------
# Environment configuration

# Exported variables so that they are available in all children jobs
export WORKDIR=$SLURM_SUBMIT_DIR

export PYTHON=python

# Load MPI module
#module load mpi/mpich-3.0-x86_64
# module load mpi/mpich-x86_64
# module load  mpi/openmpi3-x86_64
# module load petsc/3.8.3
# module load anaconda/3.7

# Before conda activate, see https://fenicsproject.discourse.group/t/fenics-from-conda-doesnt-import/3502/6
CONDA_DIR=/apps/anaconda3-python3.7
source $CONDA_DIR/etc/profile.d/conda.sh

conda activate fenicsproject

#---------------------------------------------------------------------------
# Run pogram
time mpirun -n $MPI_PROCESSES $PYTHON $WORKDIR/$PROGRAM
#time srun --mpi=none $PYTHON $WORKDIR/$PROGRAM
RESULT=$?


#---------------------------------------------------------------------------
# Save results, remove tmp files and exit
# mv $SCRATCH1 $DATAEND/$SLURM_JOB_ID
# rm -rf $SCRATCH1
exit $RESULT
