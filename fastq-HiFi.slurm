#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J fastq-HiFi
#SBATCH -o "fq.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "fq.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH -t 6:00:00
#SBATCH -N 1

# script to convert BAM to fastq, intended for HiFi data (not paired)

date
module load samtools/1.17

BAM=$1
FQ_OUT=$(basename ${BAM%.*}.fastq.gz)

echo -e "\nUsing samtools fastq to convert ${BAM} to $FQ_OUT"
echo -e "\nCommand:\n samtools fastq --threads 3 ${BAM} | gzip > $FQ_OUT"

samtools fastq --threads 3 ${BAM} | gzip > $FQ_OUT

date
samtools --version
