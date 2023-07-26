#!/bin/bash -l

#$ -l h_rt=01:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=icelake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 32

module load gcc/12

rm SWalgo_OQ

gcc -O3 -mavx2 SWalgo_OQ.c -lpthread -lgomp -o SWalgo_OQ


./SWalgo_OQ 256 96000000 32

