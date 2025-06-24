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
	prefetch $SRR
	# fastq-dump options explained: https://edwards.flinders.edu.au/fastq-dump/
	fastq-dump --outdir $SRR --skip-technical --read-filter pass --readids --dumpbase --clip --gzip --split-3 $SRR/$SRR.sra
done < $SRR_LIST

date
