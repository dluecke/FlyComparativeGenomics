#!/bin/bash
#SBATCH -A vpgru

#SBATCH -J SRAfetch
#SBATCH -o "sra.%j.%N"
#SBATCH -e "sra.%j.%N"
#SBATCH -c 8
#SBATCH -p ceres
#SBATCH -t 2-00:00:00
#SBATCH -N 1

SRR_LIST=$1

date

module load sratoolkit

while read SRR; do
	echo $SRR
	prefetch $SRR --max-size 80G
	# fasterq-dump standard options 
	fasterq-dump --outdir $SRR --skip-technical --threads 8 --split-files ./$SRR
	pigz $SRR/*.fastq # multicore gzip
done < $SRR_LIST

date
