#!/bin/bash
#SBATCH -N 6
#SBATCH --ntasks 96
#SBATCH -t 200
mpirun ./mpi-a25 --m 2688 --n 4096 --epsilon 0.001 --max-iterations 1000