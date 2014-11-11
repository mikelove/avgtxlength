#!/bin/bash 
#SBATCH -n 10  #Number of cores 
#SBATCH -t 2880  #Runtime in minutes 
#SBATCH -p general  #Partition to submit to 
#SBATCH --mem-per-cpu=5000 
#SBATCH --mail-type=ALL      #Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user= #Email to which notifications will be sent

module load centos6/cufflinks-2.2.1.Linux_x86_64 
cufflinks -o cufflinks_out/SRR577587 -p 10 --GTF /path/to/Homo_sapiens.GRCh38.77_filtered.gtf --frag-bias-correct /path/to/GRCh38.fa tophat_out/SRR577587/accepted_hits.bam
