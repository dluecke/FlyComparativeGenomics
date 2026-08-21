#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J bwa_hic_subset
#SBATCH -o "bwahicsub.%j.%N.stdout"  # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bwahicsub.%j.%N.stderr"
#SBATCH -c 48
#SBATCH -p ceres
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -C INTEL

# USAGE: sbatch path/to/yahs_subset.slurm FULL_ASSEMBLY.fa HIC_R1.fq HIC_R2.fq path/to/SUBSET.list YAHS_OUT_ID
# Starts YaHS scaffolds for contigs in SUBSET.list, using HiC mapped to full assembly
# script for generating BAM alignment for HiC reads only on SUBSET of assembly contigs
# map HiC reads to all assembly contigs, filter based SUBSET list then on HiC methods (same as OTB)
#   using ceres modules bwa_mem2 and samtools, 
#   and samblaster in $SOFTWARE_DIR

# directory with any non-module installs (samblaster)
SOFTWARE_DIR=/project/vpgru/software
GIT_REPO=~/FlyComparativeGenomics
THREADS=48

module load bwa_mem2/2.2.1
module load samtools/1.17

ASM=$1
HIC_R1=$2
HIC_R2=$3
SUBSET=$4
YAHS_OUT_ID=$5

ASM_FN=$(basename $ASM)
SUBSET_FN=$(basename $SUBSET)
SUBSET_FA=${ASM_FN%.*}.${SUBSET_FN%.*}.fa
HIC_FULL_BAM=${ASM_FN%.*}.HiC.bam
HIC_SUBSET_BAM=${ASM_FN%.*}.${SUBSET_FN%.*}.HiC.bam

date

echo -e "\nMapping HiC reads $HIC_R1 and $HIC_R2 onto $ASM using bwa-mem2:"
which bwa-mem2
echo -e "and filtering using samtools and samblaster:"
which samtools
$SOFTWARE_DIR/samblaster/samblaster --version
echo -e "\nOutput file: $HIC_SUBSET_BAM"
echo -e "\nCMD: bwa-mem2 index -p ${ASM%.*} $ASM"
echo -e "CMD: bwa-mem2 mem -5SP -t $THREADS ${ASM%.*} $HIC_R1 $HIC_R2 -o $HIC_FULL_BAM"
echo -e "CMD: samtools view -N $SUBSET -h -F 2316 $HIC_FULL_BAM |\\
 $SOFTWARE_DIR/samblaster/samblaster |\\
 samtools sort -n -@ $THREADS -o $HIC_SUBSET_BAM"
bwa-mem2 index -p ${ASM%.*} $ASM
bwa-mem2 mem -5SP -t $THREADS ${ASM%.*} $HIC_R1 $HIC_R2 -o $HIC_FULL_BAM
samtools view -N $SUBSET -h -F 2316 $HIC_FULL_BAM |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $THREADS -o $HIC_SUBSET_BAM
echo "finished bwa"
date

echo -e "\nStarting YaHS scaffolding for subset of contigs in $SUBSET"
echo -e "Writing fasta with sequences in $SUBSET\nCMD: samtools faidx $ASM -r $SUBSET > $SUBSET_FA"
samtools faidx $ASM -r $SUBSET > $SUBSET_FA

echo -e "\nSubmitting YaHS scaffolding\nCMD: sbatch $GIT_REPO/yahs_from_bam.slurm $YAHS_OUT_ID $SUBSET_FA $HIC_SUBSET_BAM"
sbatch $GIT_REPO/yahs_from_bam.slurm $YAHS_OUT_ID $SUBSET_FA $HIC_SUBSET_BAM

echo -e "Generating BAM stats files."
echo -e "CMD: samtools stats $HIC_FULL_BAM > $HIC_FULL_BAM.stats"
samtools stats $HIC_FULL_BAM > $HIC_FULL_BAM.stats
echo -e "CMD: samtools stats $HIC_SUBSET_BAM > $HIC_SUBSET_BAM.stats"
samtools stats $HIC_SUBSET_BAM > $HIC_SUBSET_BAM.stats

echo "finished"
date
