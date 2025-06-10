#!/bin/bash

# self_align-ChrsVsUnplaced.sh takes assembly ASM.fa and number of chromosomes N_CHR
# splits the largest N_CHR scaffolds into ASM-Chrs.fa
# splits the remainder into ASM-unplaced.fa
# runs mummer alignment of unplaced onto chromosome scaffolds via RunMummerAlignment-VPGRU.sh

# REQUIRES
# FlyComparativeGenomics/RunMummerAlignment-VPGRU.sh
# samtools faidx
# seqtk subseq

# USAGE
# $ self_align-ChrsVsUnplaced.sh ASM.fa N_CHR

# Containing directory for FlyComparativeGenomics/ repo
GIT_REPOS=~

module load samtools/1.17
module load seqtk/1.3

ASM=$1
N_CHR=$2

# split ASM
# index fasta for seq lengths
samtools faidx $ASM
# largest N_CHR sequences into Chrs.fa
seqtk subseq $ASM <(sort -nr -k2,2 $ASM.fai | head -n${N_CHR} | awk {'print $1'}) \
    > ${ASM%.*}-Chrs.fa
# remaining seqs into unplaced.fa
seqtk subseq $ASM <(sort -nr -k2,2 $ASM.fai | tail -n+$((N_CHR+1)) | awk {'print $1'}) \
    > ${ASM%.*}-unplaced.fa

# submit alignment job to SLURM
$GIT_REPOS/FlyComparativeGenomics/RunMummerAlignment-VPGRU.sh \
    -o ${ASM%.*}-Chrs-vs-unplaced \
    ${ASM%.*}-Chrs.fa \
    ${ASM%.*}-unplaced.fa
