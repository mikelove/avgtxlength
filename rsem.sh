#!/bin/bash
#SBATCH -n 16 #Number of cores
#SBATCH -t 1800  #Runtime in minutes
#SBATCH -p general #Partition to submit to
#SBATCH --mem-per-cpu=5000
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=  #Email to which notifications will be sent

module load centos6/rsem-1.2.11
rsem-calculate-expression --paired-end -p 16 fastq/SRR577587_1.fastq fastq/SRR577587_2.fastq /path/to/GRCh38 SRR577587
