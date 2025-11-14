#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J bwa_hic
#SBATCH -o "bwahic.%j.%N.stdout"  # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bwahic.%j.%N.stderr"
#SBATCH -c 48
#SBATCH -p ceres
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -C INTEL

# USAGE: sbatch path/to/bwa-hic.slurm ASSEMBLY.fa HIC_R1.fq HIC_R2.fq

# script for generating BAM alignment for HiC reads
# map HiC reads to assembly contigs, filter based on HiC methods (same as OTB)
#   using ceres modules bwa_mem2 and samtools, 
#   and samblaster in $SOFTWARE_DIR

# directory with any non-module installs (samblaster)
SOFTWARE_DIR=/project/vpgru/software
THREADS=48

module load bwa_mem2/2.2.1
module load samtools/1.17

ASM=$1
HIC_R1=$2
HIC_R2=$3

ASM_FN=$(basename $ASM)
HIC_BAM=${ASM_FN%.*}.HiC.bam

date

echo -e "\nMapping HiC reads $HIC_R1 and $HIC_R2 onto $ASM using bwa-mem2:"
which bwa-mem2
echo -e "and filtering using samtools and samblaster:"
which samtools
$SOFTWARE_DIR/samblaster/samblaster --version
echo -e "\nOutput file: $HIC_BAM"
echo -e "\ncommands:\n bwa-mem2 index -p ${ASM%.*} $ASM"
echo -e " bwa-mem2 mem -5SP -t $THREADS ${ASM%.*} $HIC_R1 $HIC_R2 |\
 samtools view -h -F 2316 |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $THREADS -o $HIC_BAM"
bwa-mem2 index -p ${ASM%.*} $ASM
bwa-mem2 mem -5SP -t $THREADS ${ASM%.*} $HIC_R1 $HIC_R2 |\
 samtools view -h -F 2316 |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $THREADS -o $HIC_BAM

echo -e "Generating BAM stats file. Command:\n samtools stats $HIC_BAM > $HIC_BAM.stats"
samtools stats $HIC_BAM > $HIC_BAM.stats

echo "finished"
date
