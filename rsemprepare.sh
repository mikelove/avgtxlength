#!/bin/bash
#SBATCH -n 2 #Number of cores
#SBATCH -t 1200  #Runtime in minutes
#SBATCH -p general #Partition to submit to
#SBATCH --mem-per-cpu=20000
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=  #Email to which notifications will be sent

module load centos6/rsem-1.2.11
rsem-prepare-reference --gtf /path/to/Homo_sapiens.GRCh38.77_filtered.gtf /path/to/rsem_GRCh38 GRCh38
