#!/bin/bash

# combine_haps.sh takes hap1 and hap2 assemblies, adds hap tag to seqeunce headers, and prints both to stdout

# USAGE: combine_haps.sh hap1.fa hap2.fa > diploid.fa

[ $# -eq 2 ] || { echo "USAGE: combine_haps.sh hap1.fa hap2.fa > diploid.fa"; exit; }

sed 's/>/>h1./' $1

sed 's/>/>h2./' $2
