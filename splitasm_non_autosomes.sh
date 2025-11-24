#!/bin/bash

# splitasm_non_autosomes.sh takes scaffolded assembly fasta, splits non-autosome scaffolds to contigs
#   then can redo HiC scaffolding to improve sex chromosomes/microchromosomes/unplaced scaffold regions
#   or input for reciprocal_ragtag
# outputs to stdout
# requires: samtools, seqtk, ragtag

# USAGE: splitasm_non_autosomes.sh IN_ASM.fa [N_CHR] > NONCHR_SPLIT_ASM.fa

# submit to CERES:
# sbatch splitasm_nonchr.slurm IN_ASM.fa [N_CHR]
#  (outputs to IN_ASM-nonchr_split.fa)

# preserves scaffolding for N_CHR largest scaffolds
# if no $2 provided defaults to 5 (Muscoid standard)
[[ -n $2 ]] && N_CHR=$2 || N_CHR=5

# input assembly fasta
IN_ASM=$1

samtools faidx $IN_ASM

sort $IN_ASM.fai -nr -k2,2 | \
    cut -f1 | \
    head -n $N_CHR \
    > scaffolds-autosomes.list

seqtk subseq $IN_ASM scaffolds-autosomes.list > ${IN_ASM%.*}.autosomes.fa

sort $IN_ASM.fai -nr -k2,2 | \
    cut -f1 | \
    tail -n +$((N_CHR+1)) \
    > scaffolds-rest.list

seqtk subseq $IN_ASM scaffolds-rest.list > ${IN_ASM%.*}.rest.fa

ragtag.py splitasm ${IN_ASM%.*}.rest.fa > ${IN_ASM%.*}.rest.split.fa

cat ${IN_ASM%.*}.autosomes.fa
cat ${IN_ASM%.*}.rest.split.fa

