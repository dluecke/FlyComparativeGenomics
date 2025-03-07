#!/bin/bash
#SBATCH -J samtools_depth
#SBATCH -o "sd.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "sd.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p short
#SBATCH -N 1

# script to convert BAM to fastq, intended for HiFi data (not paired)

date
module load samtools/1.17

BAM=$1

echo -e "\nUsing samtools fastq to convert ${BAM} to ${BAM%.*}.fastq"
echo -e "\nCommand:\n samtools fastq --threads 3 -o ${BAM%.*}.fastq ${BAM}"
samtools fastq --threads 3 -o ${BAM%.*}.fastq ${BAM}

date
samtools --version
