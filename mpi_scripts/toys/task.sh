#!/bin/sh
#SBATCH --time 00:10:00
#SBATCH --partition mix
#SBATCH -n64
#SBATCH --job-name ^2

mpiexec python parabola.py
