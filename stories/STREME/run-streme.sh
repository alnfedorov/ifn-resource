#!/bin/bash

#SBATCH --job-name=streme
#SBATCH --partition=short
#SBATCH --account=$(whoami)
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-4:00:00
#SBATCH --output=streme-%j.out
#SBATCH --error=streme-%j.err

# Ensure tmp directory exists
TMPDIR="${PWD}/tmp"
mkdir -p "$TMPDIR"

# Pull the Singularity image if not present
test -f memesuite_5.5.5.sif || singularity pull docker://memesuite/memesuite:5.5.5

# Run the Python script inside the Singularity container
singularity exec --bind "$(pwd)/":/project memesuite_5.5.5.sif \
    bash -c 'export TMPDIR="/project/tmp"; export OMPI_MCA_rmaps_base_oversubscribe=1; export LD_LIBRARY_PATH="/usr/lib:$LD_LIBRARY_PATH"; cd /project/ && python3 run-streme.py'
