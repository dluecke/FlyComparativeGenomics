#!/bin/bash

# splitasm_skip_list.sh takes scaffolded assembly fasta and list of scaffold names, 
# splits unlisted scaffolds to contigs
#   then can redo HiC scaffolding to improve sex chromosomes/microchromosomes/unplaced scaffold regions
#   or input for reciprocal_ragtag
#   or run through rephasing pipeline
# outputs to stdout
# requires: samtools, seqtk, ragtag

# USAGE: splitasm_skip_list.sh IN_ASM.fa skip.list > NONSKIP_SPLIT_ASM.fa

# submit to CERES:
# sbatch splitasm_skip_list.slurm IN_ASM.fa SKIP.list
#  (outputs to IN_ASM-nonSKIP_split.fa)

# input assembly fasta
IN_ASM=$1
SKIP_LIST=$2

samtools faidx $IN_ASM

IN_ASM_FN=$(basename $IN_ASM)
SKIP_LIST_FN=$(basename $SKIP_LIST)

sort $IN_ASM.fai -nr -k2,2 | \
    cut -f1 | \
    grep -F -w -f $SKIP_LIST \
    > ${IN_ASM%.*}.${SKIP_LIST_FN%.*}_scaffolds.list

seqtk subseq $IN_ASM ${IN_ASM_FN%.*}.${SKIP_LIST_FN%.*}_scaffolds.list > ${IN_ASM_FN%.*}.scaffolds.fa

sort $IN_ASM.fai -nr -k2,2 | \
    cut -f1 | \
    grep -v -w -F -f $SKIP_LIST \
    > ${IN_ASM_FN%.*}.non${SKIP_LIST_FN%.*}_scaf_to_split.list

seqtk subseq $IN_ASM ${IN_ASM_FN%.*}.non${SKIP_LIST_FN%.*}_scaf_to_split.list \
  > ${IN_ASM_FN%.*}.non${SKIP_LIST_FN%.*}_scaf_to_split.fa

ragtag.py splitasm -o ${IN_ASM_FN%.*}.${SKIP_LIST_FN%.*}_contigs.agp \
  ${IN_ASM_FN%.*}.non${SKIP_LIST_FN%.*}_scaf_to_split.fa \
  > ${IN_ASM_FN%.*}.${SKIP_LIST_FN%.*}_contigs.fa

cat ${IN_ASM%.*}.scaffolds.fa
cat ${IN_ASM%.*}.${SKIP_LIST_FN%.*}_contigs.fa
