#!/bin/bash -l

#$ -l h_rt=12:00:00 # Specify the hard time limit for the job
#$ -l cpu_arch=skylake
#$ -N SWalgo # Give job a name
#$ -P caad 
#$ -pe omp 28

module load intel/2018

./swimm2 -S preprocess -i uniprot_sprot_varsplic.fasta -o out

./swimm2 -S search -q uniprot_sprot.fasta -d out

