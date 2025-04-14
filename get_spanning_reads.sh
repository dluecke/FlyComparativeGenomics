#!/bin/bash

# get_spanning_reads.sh writes SAM file with reads spanning a target sequence/coordinate combination
# intended to extract HiFi reads covering genomic regions of interest
# takes positional variables:
#  MAPPING.bam file with reads mapped to sequence
#  target SEQuence (contig/scaffold/chromosome)
#  target COORDinate (bp on target sequence)
#  SPAN length (distance from target coordinate required to cover on both sides, default 1kb)
# writes:
#  MAPPING-SEQ_COORDspanSPANbp.sam file with all reads spanning target coordinate +/- span distance

# REQUIRES: samtools sort (written with samtools 1.17)

# positional variables
IN_BAM=$1
TARGET_SEQ=$2
TARGET_COORD=$3
# default span lenth to 1kb if not provided
[[ -n $4 ]] && SPAN_LENGTH=$4 || SPAN_LENGTH=1000

# output filename
OUT_SAM="${IN_BAM%.*}-${TARGET_SEQ}_${TARGET_COORD}span${SPAN_LENGTH}bp.sam"

# samtools view to extract reads covering target seq/coord
samtools view $IN_BAM ${TARGET_SEQ}:${TARGET_COORD}:${TARGET_COORD} |\
# awk to filter for reads covering full span
awk -v t=$TARGET_COORD -v s=$SPAN_LENGTH '$4 <= t-s && $4+$9 >= t+s' > $OUT_SAM
