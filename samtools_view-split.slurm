#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="split.bam"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p ceres
#SBATCH -t 2:00:00
#SBATCH -o "split.%j.%N.stdout"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "split.%j.%N.stderr"           #optional but it prints our standard error

# writes BAM file for all reads in input BAM $1 with a split mapping anywhere
# uses 2-pass strategy to get both primary and supplementary (split 0x800) alignments

ml samtools

BAM=$1

OUT_BAM=${BAM%.*}.split.bam

# list of reads with split mapping
samtools view -f 0x800 $BAM | cut -f1 | sort -u > split_reads-${BAM%.*}.txt

# all alignments from split reads
samtools view -N split_reads-${BAM%.*}.txt -b -o $OUT_BAM $BAM

# index output for IGV viewing
samtools index $OUT_BAM
