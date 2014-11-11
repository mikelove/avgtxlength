#!/bin/bash 
#SBATCH -n 4  #Number of cores 
#SBATCH -t 600  #Runtime in minutes 
#SBATCH -p general
#SBATCH --mem-per-cpu=5000 
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=  #Email to which notifications will be sent

module load centos6/parallel-20130922_gcc-4.4.7
module load bio/sratoolkit.2.3.3-4
ls sra/* | parallel -j 4 fastq-dump --split-files {}
