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
	# fasterq-dump options from NCBI EGAPx (but modified for fastq): 
    #       https://github.com/ncbi/egapx/tree/main?tab=readme-ov-file#offline-mode
	fasterq-dump --outdir $SRR --skip-technical --threads 8 --split-files \
      --seq-defline "@\$ac.\$si.\$ri" --qual-defline "+\$ac.\$si.\$ri"  ./$SRR
done < $SRR_LIST

date
