#!/bin/bash

# get_variants_by_sequence.sh takes FILE.vcf, 
# writes FILE.nvar_by_seq.csv with seqname, n_variants, seq_length

# USAGE: get_variants_by_sequence.sh FILE.vcf

VCF=$1
CSV=${VCF%.*}.nvar_by_seq.csv

while read -r -a line; do 
    seq=${line[1]}
    n_var=${line[0]}
    len=$(grep -m1 $seq $VCF | cut -d'=' -f4 | tr -d '>')
    echo $seq,$n_var,$len
done < <(grep -a -v ^# $VCF | cut -f1 | sort | uniq -c) \
 > $CSV
