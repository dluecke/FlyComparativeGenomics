#!/bin/bash

# stv-region_bam.sh quick shell to write BAM file for specified region
# takes BAM (or SAM) file with headers, and region in samtools view format seq:c1-c2
# writes BAM "inbamname-seq_c1toc2.bam"
# requires samtools, written with Ceres module 1.17
# run from salloc node with samtools loaded

# USAGE: stv-region_bam.sh $IN_BAM "seq:c1-c2"

IN_ALN=$1
REGION=$2

[[ -n $REGION ]] || \
  { echo "USAGE: sbatch samtools_view-region.slurm IN.SAM|IN.BAM seq:c1-c2"; exit; }

# convert region syntax seq_N:0-100 to filename syntax seqN_0to100
REG_TAG=$(echo $REGION | tr -d '_' | tr ':' '_' | sed -E 's/-([0-9])/to\1/')

# output BAM
OUT_BAM=${IN_ALN%.*}-${REG_TAG}.bam

# check for index of alignment file, use -c to avoid max seq length issue
if [ ! -f ${IN_ALN}.csi ]; then
    samtools index -c $IN_ALN
fi

# write output bam file, skip secondary alignments 
samtools view -F 256 -b -o $OUT_BAM $IN_ALN $REGION

# index output bam
samtools index $OUT_BAM
