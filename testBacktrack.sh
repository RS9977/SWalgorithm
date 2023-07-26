#!/bin/bash -l

#$ -l h_rt=01:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=cascadelake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 1

lscpu

module load gcc/12

gcc -mavx2 -O3 SWalgo_V4.c -lgomp -o SWalgo_V4

gcc -mavx512f -O3 SWalgo_V3.c -lgomp -o SWalgo_V3

gcc -mavx2 -O3 SWalgo_V4b.c -lgomp -o SWalgo_V4b

gcc -mavx512f -O3 SWalgo_V3b.c -lgomp -o SWalgo_V3b

echo "AVX512 -----------------------------------------------------------"

./SWalgo_V3 16 10000
./SWalgo_V3b 16 10000

./SWalgo_V3 32 10000
./SWalgo_V3b 32 10000

./SWalgo_V3 64 10000
./SWalgo_V3b 64 10000

./SWalgo_V3 128 10000
./SWalgo_V3b 128 10000

./SWalgo_V3 256 10000
./SWalgo_V3b 256 10000

./SWalgo_V3 1024 10000
./SWalgo_V3b 1024 10000

./SWalgo_V3 4096 10000
./SWalgo_V3b 4096 10000

echo "AVX2 -----------------------------------------------------------"

./SWalgo_V4 16 10000
./SWalgo_V4b 16 10000

./SWalgo_V4 32 10000
./SWalgo_V4b 32 10000

./SWalgo_V4 64 10000
./SWalgo_V4b 64 10000

./SWalgo_V4 128 10000
./SWalgo_V4b 128 10000

./SWalgo_V4 256 10000
./SWalgo_V4b 256 10000

./SWalgo_V4 1024 10000
./SWalgo_V4b 1024 10000

./SWalgo_V4 4096 10000
./SWalgo_V4b 4096 10000