#!/bin/bash 
#SBATCH -n 20  #Number of cores 
#SBATCH -t 1800  #Runtime in minutes 
#SBATCH -p general  #Partition to submit to 
#SBATCH --mem-per-cpu=5000 
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=  #Email to which notifications will be sent

module load centos6/tophat-2.0.11.Linux_x86_64
tophat2 -o tophat_out/SRR577587 -p 20 /path/to/GRCh38 fastq/SRR577587_1.fastq fastq/SRR577587_2.fastq

