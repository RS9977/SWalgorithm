#!/bin/bash -l

#$ -l h_rt=01:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=cascadelake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 32

lscpu

module load gcc/12

gcc -mavx2 -O3 SWalgo_V9.c -lgomp -lpthread -o SWalgo_V9

./SWalgo_V9 114688 140000 32
./SWalgo_V9 114688 140000 16
./SWalgo_V9 114688 140000 8
./SWalgo_V9 114688 140000 4
./SWalgo_V9 114688 140000 2
./SWalgo_V9 114688 140000 1