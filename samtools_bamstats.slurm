#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J samtools_bamstats
#SBATCH -o "sbs.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "sbs.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p short
#SBATCH -N 1
date
module load samtools/1.17

BAM=$1

date

echo -e "\nGenerating stats file for $BAM"
echo -e "\nUsing Samtools versions:"
samtools --version
echo -e "\nCommand:\n samtools stats $BAM > $BAM.stats"
samtools stats $BAM > $BAM.stats

echo
date