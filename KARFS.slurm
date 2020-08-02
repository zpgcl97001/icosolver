#!/bin/bash

#SBATCH --job-name=Re100
#SBATCH --partition=small
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -n 1 

ulimit -s unlimited
ulimit -l unlimited

srun --mpi=pmi2 CouetteFlow >log.txt
