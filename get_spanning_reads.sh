#!/bin/bash

# get_spanning_reads.sh writes SAM file with reads spanning a target sequence/coordinate combination
# intended to extract HiFi reads covering genomic regions of interest
# takes positional variables:
#  MAPPING.bam file with reads mapped to sequence
#  target SEQuence (contig/scaffold/chromosome)
#  target COORDinate (bp on target sequence)
#  SPAN length (distance from target coordinate required to cover on both sides, default 5kb)
# writes:
#  MAPPING-SEQ_COORDspanSPANbp.sam file with all reads spanning target coordinate +/- span distance

# REQUIRES: samtools sort (written with samtools 1.17)

# positional variables
IN_BAM=$1
TARGET_SEQ=$2
TARGET_COORD=$3
# default span lenth to 5kb if not provided
[[ -n $4 ]] && SPAN_LENGTH=$4 || SPAN_LENGTH=5000

# output filename
OUT_SAM="${IN_BAM%.*}-${TARGET_SEQ}_${TARGET_COORD}span${SPAN_LENGTH}bp.sam"

# screen output
echo -e "Writing SAM file with reads aligned across specific coordinate"
echo -e " target sequence name and coordinate: $TARGET_SEQ:$TARGET_COORD"
echo -e " read bp span in both directions: $SPAN_LENGTH"
echo -e " input alignment file: $IN_BAM"
echo -e " output alignment file: $OUT_SAM\n"
echo -e "Using samtools, version info:"
samtools --version
echo -e "\nCommand:\n samtools view $IN_BAM ${TARGET_SEQ}:${TARGET_COORD}-${TARGET_COORD} |\
 awk -v t=$TARGET_COORD -v s=$SPAN_LENGTH '\$4 <= t-s && \$4+\$9 >= t+s' > $OUT_SAM"

# samtools view to extract reads covering target seq/coord
samtools view $IN_BAM ${TARGET_SEQ}:${TARGET_COORD}-${TARGET_COORD} |\
# awk to filter for reads covering full span
awk -v t=$TARGET_COORD -v s=$SPAN_LENGTH '$4 <= t-s && $4+$9 >= t+s' > $OUT_SAM
