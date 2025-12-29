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

date
ml samtools

BAM=$1

OUT_BAM=${BAM%.*}.split.bam

echo "Extracting all reads with supplementary alignments (split mappings) from $BAM into $OUT_BAM"

# list of reads with split mapping
echo "Write all read IDs with supplementary alignments to split_reads-${BAM%.*}.txt"
echo "CMD: samtools view -f 0x800 $BAM | cut -f1 | sort -u > split_reads-${BAM%.*}.txt"
samtools view -f 0x800 $BAM | cut -f1 | sort -u > split_reads-${BAM%.*}.txt

# all alignments from split reads
echo "Write bam file $OUT_BAM with primary and supplementary alignments for split reads"
echo "CMD: samtools view -N split_reads-${BAM%.*}.txt -b -o $OUT_BAM $BAM"
samtools view -N split_reads-${BAM%.*}.txt -b -o $OUT_BAM $BAM

# index output for IGV viewing
echo "indexing $OUT_BAM"
echo "CMD: samtools index $OUT_BAM"
samtools index $OUT_BAM

samtools --version
date
