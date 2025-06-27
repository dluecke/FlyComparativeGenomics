#!/bin/bash

# splitasm_non_autosomes.sh takes scaffolded assembly fasta, splits non-autosome scaffolds to contigs
#   then can redo HiC scaffolding to improve sex chromosomes/microchromosomes/unplaced scaffold regions
# outputs to stdout
# requires: samtools, seqtk, ragtag

# preserves scaffolding for N_CHR largest scaffolds
N_AUTOSOME=5

# input assembly fasta
IN_ASM=$1

# output fasta 

# output named 

samtools faidx $IN_ASM

sort $IN_ASM.fai -nr -k2,2 | \
    cut -f1 | \
    head -n $N_AUTOSOME \
    > scaffolds-autosomes.list

seqtk subseq $IN_ASM scaffolds-autosomes.list > ${IN_ASM%.*}.autosomes.fa

sort $IN_ASM.fai -nr -k2,2 | \
    cut -f1 | \
    tail -n $((N_AUTOSOME+1)) \
    > scaffolds-rest.list

seqtk subseq $IN_ASM scaffolds-rest.list > ${IN_ASM%.*}.rest.fa

ragtag.py splitasm ${IN_ASM%.*}.rest.fa > ${IN_ASM%.*}.rest.split.fa

cat ${IN_ASM%.*}.autosomes.fa
cat ${IN_ASM%.*}.rest.split.fa

