#!/bin/bash

# ragtag.sh re-scaffolds QRY based on alignment to REF and read mapping
# runs ragtag correct with QRY read mapping and REF alignment to break misjoins
# runs ragtag scaffold on corrected-QRY to reorder based on REF alignment
# needs at least 32 threads

# USAGE: ragtag.sh REFERENCE.fa QUERY.fa QUERY_READS.fa|fq

module load ragtag/2.1.0

REF=$1
QRY=$2
READS=$3

# run correct to break contigs if alignment and reads support
# -u to tag unmodified sequences (better for AGP)
# -T corr for high quality long reads
# output file (default): ragtag_output/ragtag.correct.fasta
ragtag.py correct -u -R $READS -T corr -t 32 $REF $QRY

# run scaffold to reorder corrected QRY scaffolds
# -r to estimate gap size, -u to tag unmodified sequences (better for AGP)
ragtag.py scaffold -r -t 32 $REF ragtag_output/ragtag.correct.fasta
