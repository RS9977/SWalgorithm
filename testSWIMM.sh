#!/bin/bash -l

core=28

#$ -l h_rt=00:15:00 # Specify the hard time limit for the job
#$ -l cpu_arch=icelake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 32

lscpu

module load intel/2018

module load gcc/12

make clean

make AVX2

./swimm2 -S preprocess -i uniprot_sprot_varsplic.fasta -o out

./swimm2 -S search -q uniprot_sprot.fasta -d out

./swimm2 -S preprocess -i test3.fasta -o out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

./swimm2 -S search -q test2.fasta -d out

gcc -O3 -mavx2 SWalgo_OQ.c -lpthread -lgomp -o SWalgo_OQ

./SWalgo_OQ 256 96000000 32
