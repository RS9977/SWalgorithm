#!/bin/bash -l

#$ -l h_rt=00:10:00 # Specify the hard time limit for the job
#$ -l cpu_arch=skylake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 28


module load gcc/12

gcc -mavx512bw -O3 SWalgo_V11.c -lpthread -lgomp -o SWalgo_V11

./SWalgo_V11 2048 100000
