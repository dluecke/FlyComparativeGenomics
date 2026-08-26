#!/bin/bash

# ragtag_scaffold.sh runs ragtag scaffold

# USAGE: ragtag_scaffold.sh REFERENCE.fa QUERY.fa

module load ragtag/2.1.0

echo "Starting RagTag scaffolding"
date
echo -e "\nRagTag loaded: $(ragtag.py --version)\n"

REF=$1
QRY=$2
[ -z $3 ] && THREADS=32 || THREADS=$3

# run correct to break contigs if alignment and reads support
# -r to estimate gap size, -u to tag unmodified sequences (better for AGP)
echo "CMD: ragtag.py scaffold -r -u -t $THREADS $REF $QRY"
ragtag.py scaffold -r -u -t $THREADS $REF $QRY

