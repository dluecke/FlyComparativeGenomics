#!/bin/bash

# minimap2-ONT.sh maps Q20 ONT reads to reference

module load minimap2/2.24

REF=$1
ONT=$2

OUT_SAM=${REF%.*}-${ONT%.*}-mm2.sam

minimap2 -a -k19 -w19 -U50,500 -g10k $REF $ONT > $OUT_SAM
