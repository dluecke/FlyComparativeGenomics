#!/bin/bash

# ragtag_scaffold.sh runs ragtag scaffold

# USAGE: ragtag_scaffold.sh REFERENCE.fa QUERY.fa

module load ragtag/2.1.0

REF=$1
QRY=$2

# run correct to break contigs if alignment and reads support
# -r to estimate gap size, -u to tag unmodified sequences (better for AGP)
ragtag.py scaffold -r -t 32 $REF $QRY
