##!/bin/bash
#SBATCH --job-name=diffusion-simple.py
#SBATCH --partition=cn
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=1
#SBATCH --error=./out/job.%J.err
#SBATCH --output=./out/job.%J.out
#SBATCH --mail-user=daniel.acosta@uca.es
#SBATCH --mail-type=ALL
#SBATCH --requeue

#------------------------------------
# Global variables
export PROGRAM=./diffusion-simple.py # Change this variable in #SBATCH --job-name too
export MPI_PROCESSES=10 # Change this variable in #SBATCH --nodes too

# REMEMBER TO MAKE A DIRECTORY CALLED out INSIDE THE WORK DIRECTORY TO SAVE THE OUTPUTS

#---------------------------------------------------------------------------
# Environment configuration

# Exported variables so that they are available in all children jobs
export WORKDIR=$SLURM_SUBMIT_DIR
export PYTHON=python

# module load anaconda/3.7

# Before conda activate, see https://fenicsproject.discourse.group/t/fenics-from-conda-doesnt-import/3502/6
CONDA_DIR=/apps/anaconda3-python3.7
source $CONDA_DIR/etc/profile.d/conda.sh

conda activate fenicsproject

#---------------------------------------------------------------------------
# Run pogram
time mpirun -n $MPI_PROCESSES $PYTHON $WORKDIR/$PROGRAM
RESULT=$?

#---------------------------------------------------------------------------
# Save results, remove tmp files and exit
# mv $SCRATCH1 $DATAEND/$SLURM_JOB_ID
# rm -rf $SCRATCH1
exit $RESULT
