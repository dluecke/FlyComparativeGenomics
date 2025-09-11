#!/bin/bash

# order_chromosomes.sh converts chr-scale scaffolds to properly oriented chromosomes
# takes SCAFFOLDS.fa and chr_assignment.tsv, outputs SCAFFOLDS-Chromosomes.fa
# chr_assignment.tsv format is Chromosome/Scaffold/Orientation:
#   Chromosome1 scaffold_i  F
#   Chromosome2 scaffold_j  R
# partial scaffold name match OK

# REQUIRES: samtools, written with version 1.17
# USAGE: order_chromosomes.sh scaffolds.fa chr_assignment.tsv

usage() {
    echo -e "\nUSAGE: order_chromosomes.sh SCAFFOLDS.fa CHR_ASSIGNMENT.tsv"
    echo "CHR_ASSIGNMENT.tsv format (ChromosomeN/scaffold_i/FRorientation) eg:"
    echo "  Chromosome1 scaffold_j  R"
    echo "  Chromosome2 scaffold_k  F"
    echo "Will take longest partial match of provided scaffold ID"
    echo -e "Outputs assembly file:\n SCAFFOLD-Chromosomes.fa"
    echo -e "and exact scaffold-to-chromosome relationships:\n SCAFFOLDS-chr_assignment.tsv"
    echo "requires samtools faidx"
    exit
}
[[ $# == 2 ]] || usage

# inputs
SCAFFOLDS=$1
ORDERING=$2

# output files
# new assembly
CHROMOSOMES=${SCAFFOLDS%.*}-Chromosomes.fa
# assignments with full scaffold names
ORDER_OUT=${SCAFFOLDS%.*}-chr_assignment.tsv

# index input fasta
samtools faidx $SCAFFOLDS || usage

echo -e "\nConverting scaffolds assembly:\n $SCAFFOLDS"
echo -e "into chromosome assembly:\n $CHROMOSOMES"
echo -e "using renaming and orienting information in:\n $ORDERING"

# sort by sequence length 
sort -nr -k2,2 $SCAFFOLDS.fai > $SCAFFOLDS.fai.sort


### Write chromosomes to new assembly

# loop through ORDERING, each line as tab sep array
while IFS=$'\t' read -r -a arr_ORDERING; do

    # Chromosome name
    CHRi=${arr_ORDERING[0]}

    # check orientation
    ORIENTi=${arr_ORDERING[2]}
    if [ "$ORIENTi" = "F" ]; then
        RC_TAG=""
    elif [ "$ORIENTi" = "R" ]; then
        RC_TAG="--reverse-complement "
    else
        echo -e "\nDon't recognize orientation, should be F or R"
        echo -e "Problem line in $ORDERING:\n$line"
        usage
    fi

    # get full scaffold name, longest sequence starting with input scaffold_i
    SCAFi=$(grep -m1 ${arr_ORDERING[1]} $SCAFFOLDS.fai.sort | cut -f1)

    # write sequence to output assembly
    echo -e "\nWriting sequence $SCAFi as $CHRi in file $CHROMOSOMES"
    echo 'CMD: samtools faidx --mark-strand no '${RC_TAG}${SCAFFOLDS} $SCAFi' | sed "s/'$SCAFi/$CHRi'/" | sed "s@/rc@@" >> '$CHROMOSOMES
    # extract SCAFi
    samtools faidx --mark-strand no ${RC_TAG}${SCAFFOLDS} $SCAFi | \
    # change name to ChromosomeN and append to new assembly file
        sed "s/$SCAFi/$CHRi/" >> $CHROMOSOMES

    # write line for full scaffold name chr_assignment.tsv output
    # both a good record and will be used to extract non-chromosome scaffolds
    echo -e "$CHRi\t$SCAFi\t$ORIENTi" >> $ORDER_OUT

done < $ORDERING


### Write remaining scaffolds to new assembly

# list of all scaffolds, will remove those already extracted
cut -f1 $SCAFFOLDS.fai.sort > non-chr-scafs.list

# remove scaffolds already extracted as recorded in ORDER_OUT
while read rm_scaf; do
    sed -i "/$rm_scaf/d" non-chr-scafs.list
done < <(cut -f2 $ORDER_OUT)

# append all remaining scaffolds to CHROMOSOMES assembly
echo -e "\nWriting remaining scaffolds to assembly"
echo "CMD: samtools faidx -r non-chr-scafs.list $SCAFFOLDS >> $CHROMOSOMES"
samtools faidx -r non-chr-scafs.list $SCAFFOLDS >> $CHROMOSOMES

echo -e "\nChromosome assembly $CHROMOSOMES finished.\nSee $ORDER_OUT for scaffold-to-chromosome relationships"
