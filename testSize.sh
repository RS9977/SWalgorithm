#!/bin/bash -l

#$ -l h_rt=01:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=cascadelake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 32

lscpu

module load gcc/12


gcc -mavx2 -O3 SWalgo_V9.c -lgomp -lpthread -o SWalgo_V9

gcc -mavx2 -O3 SWalgo_V8.c -lgomp -lpthread -o SWalgo_V8

module load gcc/9

gcc -mavx512f -O3 SWalgo_V7.c -lpthread -lgomp -o SWalgo_V7

gcc -mavx512bw -O3 SWalgo_V11.c -lpthread -lgomp -o SWalgo_V11

echo "AVX2 1B"

./SWalgo_V9 4096 131200 32

./SWalgo_V9 1024 131200 32

./SWalgo_V9 512 131200 32

./SWalgo_V9 256 131200 32

./SWalgo_V9 128 131200 32

./SWalgo_V9 64 131200 32

./SWalgo_V9 32 131200 32

echo "AVX2 4B"

./SWalgo_V8 4096 131200 32

./SWalgo_V8 1024 131200 32

./SWalgo_V8 512 131200 32

./SWalgo_V8 256 131200 32

./SWalgo_V8 128 131200 32

./SWalgo_V8 64 131200 32

./SWalgo_V8 32 131200 32

echo "AVX512 1B"

./SWalgo_V11 4096 131200 32

./SWalgo_V11 1024 131200 32

./SWalgo_V11 512 131200 32

./SWalgo_V11 256 131200 32

./SWalgo_V11 128 131200 32

./SWalgo_V11 64 131200 32

./SWalgo_V11 32 131200 32

echo "AVX512 4B"

./SWalgo_V7 4096 131200 32

./SWalgo_V7 1024 131200 32

./SWalgo_V7 512 131200 32

./SWalgo_V7 256 131200 32

./SWalgo_V7 128 131200 32

./SWalgo_V7 64 131200 32

./SWalgo_V7 32 131200 32