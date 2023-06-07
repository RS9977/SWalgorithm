#!/bin/bash -l

#$ -l h_rt=03:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=skylake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 1


module load gcc/9.3.0

gcc -O1 -lgomp SWalgo.c -o SWalgo

gcc -O1 -lgomp SWalgo_V2.c -o SWalgo_V2

gcc -O1 -mavx512f -lgomp SWalgo_V3.c -o SWalgo_V3

gcc -O1 -mavx2 -lgomp SWalgo_V4.c -o SWalgo_V4

./SWalgo 1000 100000
./SWalgo_V2 1000 100000
./SWalgo_V3 1000 100000
./SWalgo_V4 1000 100000