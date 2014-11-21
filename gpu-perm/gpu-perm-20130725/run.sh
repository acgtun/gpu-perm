#! /bin/csh
#PBS -l walltime=00:10:00
#PBS -l nodes=1:ppn=2:gpus=1
#PBS -d .

source /usr/usc/cuda/5.0/setup.csh

./gpu-perm NC_008253.fasta e_coli_10000snp.fa 


