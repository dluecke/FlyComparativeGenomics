#!/bin/bash

# get number of variant sites per window from a vcf file

WINDOW=100000
VCFFILE=$1


grep -a -v ^# $VCFFILE | awk -v win=$WINDOW -v OFS='\t' '{w=int($2/win)} {print $1"-"w}' | \
    sort | uniq -c | awk -v OFS='\t' '{print $2, $1}' | \
    awk -v OFS='\t' '{split($1, arr, "-")} {print length(arr[1]), arr[1], arr[2], $0}' | \
    sort -k1,1n -k2,2 -k3,3n | cut -f2-5 > $VCFFILE.nvar_by_windows.tsv