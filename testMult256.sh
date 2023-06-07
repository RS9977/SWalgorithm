#!/bin/bash -l

#$ -l h_rt=00:15:00 # Specify the hard time limit for the job
#$ -l cpu_arch=skylake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 28


module load gcc/9.3.0

gcc -mavx512f -O3 SWalgo_V8.c -lpthread -lgomp -o SWalgo_V8

./SWalgo_V8 1000 100000
