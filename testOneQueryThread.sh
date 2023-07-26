#!/bin/bash -l

#$ -l h_rt=01:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=broadwell
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 28

gcc -O3 -mavx2 SWalgo_V5a.c -lpthread -lgomp -o SWalgo_V5a


./SWalgo_V5a 100 1000 1

./SWalgo_V5a 100 1000 2

./SWalgo_V5a 100 1000 4

./SWalgo_V5a 100 1000 8

./SWalgo_V5a 100 1000 16

./SWalgo_V5a 100 1000 28