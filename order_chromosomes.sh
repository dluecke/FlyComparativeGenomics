#!/bin/bash

# order_chromosomes.sh converts chr-scale scaffolds to properly oriented chromosomes
# also removes leading/trailing Ns and wraps line lengths to 80
# takes SCAFFOLDS.fa and chr_assignment.tsv, outputs SCAFFOLDS-Chromosomes.fa
# chr_assignment.tsv format is Chromosome/Scaffold/Orientation:
#   Chromosome1 scaffold_i  F
#   Chromosome2 scaffold_j  R
# partial scaffold name match OK

# REQUIRES: samtools, written with version 1.17
# USAGE: order_chromosomes.sh scaffolds.fa chr_assignment.tsv

usage() {
    echo -e "USAGE: order_chromosomes.sh SCAFFOLDS.fa CHR_ASSIGNMENT.tsv [ > out.log]"
    echo "CHR_ASSIGNMENT.tsv format (ChromosomeN/scaffold_i/FRorientation) eg:"
    echo "  Chromosome1 scaffold_j  R"
    echo "  Chromosome2 scaffold_k  F"
    echo "Will take longest partial match of provided scaffold ID"
    echo -e "Outputs assembly files:\n SCAFFOLD-Chromosomes.all.fa\n SCAFFOLD-Chromosomes.fa (purge sequences <1kb)"
    echo -e "and exact scaffold-to-chromosome relationships:\n SCAFFOLDS-chr_assignment.tsv"
    echo "NOTE: existing copies of these files will be overwritten"
    echo "requires samtools faidx"
    exit
}
[[ $# == 2 ]] || usage

# threshold for purging final sequences
PURGE_BELOW_BP=999

# inputs
SCAFFOLDS=$1
ORDERING=$2

# output files
# new assembly
CHROMOSOMES=${SCAFFOLDS%.*}-Chromosomes
# assignments with full scaffold names
ORDER_OUT=${SCAFFOLDS%.*}-chr_assignment.tsv
# list of non-chromosome scaffolds
NON_SCAF_LIST=${SCAFFOLDS%.*}-non_chr-scafs.list

# index input fasta
samtools faidx $SCAFFOLDS || usage

echo -e "\nConverting scaffolds assembly:\n $SCAFFOLDS"
echo -e "into chromosome assembly:\n $CHROMOSOMES.fa"
echo -e "using renaming and orienting information in:\n $ORDERING"

# sort by sequence length 
sort -nr -k2,2 $SCAFFOLDS.fai > $SCAFFOLDS.fai.sort


### Write chromosomes to new assembly

# delete prior output file (all output is appended)
> $CHROMOSOMES.all.fa
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
    echo -e "\nWriting sequence $SCAFi as $CHRi in file $CHROMOSOMES.all.fa"
    echo 'CMD: samtools faidx --mark-strand no '${RC_TAG}${SCAFFOLDS} $SCAFi' | sed "s/'$SCAFi/$CHRi'/" >> '$CHROMOSOMES.all.fa
    # extract SCAFi
    samtools faidx --mark-strand no ${RC_TAG}${SCAFFOLDS} $SCAFi | \
    # change name to ChromosomeN and append to new assembly file
        sed "s/$SCAFi/$CHRi/" >> $CHROMOSOMES.all.fa

    # write line for full scaffold name chr_assignment.tsv output
    # both a good record and will be used to extract non-chromosome scaffolds
    echo -e "$CHRi\t$SCAFi\t$ORIENTi" >> $ORDER_OUT

done < $ORDERING


### Write remaining scaffolds to new assembly

# list of all scaffolds, will remove those already extracted
cut -f1 $SCAFFOLDS.fai.sort > $NON_SCAF_LIST

# remove scaffolds already extracted as recorded in ORDER_OUT
while read rm_scaf; do
    sed -i "/$rm_scaf/d" $NON_SCAF_LIST
done < <(cut -f2 $ORDER_OUT)

# append all remaining scaffolds to CHROMOSOMES assembly
echo -e "\nWriting remaining scaffolds to assembly, see $NON_SCAF_LIST for sequence IDs"
echo "CMD: samtools faidx -r $NON_SCAF_LIST $SCAFFOLDS >> $CHROMOSOMES.all.fa"
samtools faidx -r $NON_SCAF_LIST $SCAFFOLDS >> $CHROMOSOMES.all.fa


### Clean up sequences: remove leading/trailing Ns and wrap to 80 bp per line, purge small sequences
# code cribbed from https://www.biostars.org/p/412636/#9494761
echo -e "\nRemoving any leading or trailing N/n characters and wrapping line length"
echo "CMD: awk '/^>/ "'{printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}'"'\
 $CHROMOSOMES.all.fa "'| tr "\t" "\n" |'" sed -r '/^>/! s/[Nn]+$|^[Nn]+//g' | fold -w 80 > $CHROMOSOMES.clean.tmp; mv $CHROMOSOMES.clean $CHROMOSOMES.all.fa"

# Linearize sequences to single line
awk '
  # Header lines
  /^>/ {
    # Print the header line followed by a tab, with leading newline if not first record
    printf("%s%s\t", (N>0?"\n":""), $0)
    N++
    next # skip next block for header lines
  }
  { printf("%s", $0) }    # Print the sequence line without newline  
  END { printf("\n") }    # Final newline after last record   
' $CHROMOSOMES.all.fa | tr "\t" "\n" | \
# Remove leading or trailing N or n from non-header lines
  sed -r '/^>/! s/[Nn]+$|^[Nn]+//g' | \
 # Wrap to 80 character lines (will also wrap header lines, but header will only be "ChromosomeN" or scaffold ID so should be fine)
  fold -w 80 > $CHROMOSOMES.clean.tmp
mv $CHROMOSOMES.clean.tmp $CHROMOSOMES.all.fa

echo -e "\nWriting final assembly with cleaned sequences >$PURGE_BELOW_BP bps"
echo "CMD: samtools faidx $CHROMOSOMES.all.fa; awk -v bpthreshold=$PURGE_BELOW_BP '"'$2 > bpthreshold {print $1}'"'\
 $CHROMOSOMES.all.fa.fai > $CHROMOSOMES.final_seqs.txt;\
 samtools faidx -r $CHROMOSOMES.final_seqs.txt $CHROMOSOMES.all.fa > $CHROMOSOMES.fa"
# index clean assembly
samtools faidx $CHROMOSOMES.all.fa

# get sequences longer than 1kb
awk -v bpthreshold=$PURGE_BELOW_BP '$2 > bpthreshold {print $1}' $CHROMOSOMES.all.fa.fai > $CHROMOSOMES.final_seqs.txt

# write sequences above threshold to final assembly
samtools faidx -r $CHROMOSOMES.final_seqs.txt $CHROMOSOMES.all.fa > $CHROMOSOMES.fa

# finished
echo -e "\nChromosome assembly $CHROMOSOMES.fa finished.\nSee $ORDER_OUT for scaffold-to-chromosome relationships and $CHROMOSOMES.all.fa for unpurged assembly"
date
