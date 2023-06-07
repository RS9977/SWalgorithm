#!/bin/bash -l

#$ -l h_rt=00:15:00 # Specify the hard time limit for the job
#$ -l cpu_arch=skylake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 1


module load gcc/9.3.0

gcc -O3 -mavx512f -lgomp SWalgo_V3.c -o SWalgo_V3

./SWalgo_V3 1000 100000