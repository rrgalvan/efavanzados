#!/bin/bash

# Before conda activate, see https://fenicsproject.discourse.group/t/fenics-from-conda-doesnt-import/3502/6
CONDA_DIR=/apps/anaconda3-python3.7
source $CONDA_DIR/etc/profile.d/conda.sh
conda activate fenicsproject

time mpirun -n 16 python diffusion-simple.py
