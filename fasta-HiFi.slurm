#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J fasta-HiFi
#SBATCH -o "fq.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "fq.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH -t 6:00:00
#SBATCH -N 1

# script to convert BAM to fasta, intended for HiFi data (not paired)

date
module load samtools/1.17

BAM=$1
FA_OUT=$(basename ${BAM%.*}.fa.gz)

echo -e "\nUsing samtools fastq to convert ${BAM} to $FA_OUT"
echo -e "\nCMD: samtools fasta --threads 3 -o $FA_OUT ${BAM}"

samtools fasta --threads 3 -o $FA_OUT ${BAM}

date
samtools --version
