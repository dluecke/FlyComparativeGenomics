#!/bin/bash

# ragtag_correct.sh runs ragtag correct with reference reads

# USAGE: ragtag_correct.sh REFERENCE.fa QUERY.fa QUERY_READS.f?

module load ragtag/2.1.0

REF=$1
QRY=$2
READS=$3

ragtag.py correct -R $READS -T corr $REF $QRY
