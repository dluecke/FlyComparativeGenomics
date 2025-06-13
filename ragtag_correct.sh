#!/bin/bash

# ragtag_correct.sh runs ragtag correct with reference reads
# needs at least 32 threads

# USAGE: ragtag_correct.sh REFERENCE.fa QUERY.fa QUERY_READS.fa|fq

module load ragtag/2.1.0

REF=$1
QRY=$2
READS=$3

# split the scaffolds first, prevent trailing Ns in output from correct
ragtag.py splitasm $QRY > ${QRY%.*}.split.fa

# run correct to break contigs if alignment and reads support
ragtag.py correct -R $READS -T corr -t 32 $REF ${QRY%.*}.split.fa
