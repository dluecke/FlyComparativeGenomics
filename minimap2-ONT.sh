#!/bin/bash

# minimap2-ONT.sh maps Q20 ONT reads to reference

module load minimap2/2.24

REF=$1
ONT=$2

OUT_NAME=${REF%.*}-${ONT%.*}-mm2

minimap2 -a -k19 -w19 -U50,500 -g10k $REF $ONT > $OUT_NAME.sam

paftools.js sam2paf ${OUT_NAME}.sam | gzip -c - > ${OUT_NAME}.paf.gz
