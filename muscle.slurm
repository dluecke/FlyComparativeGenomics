#!/bin/bash
#SBATCH --job-name="muscle_align"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p short
#SBATCH -o "ma.stdout.%j.%N"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "ma.stderr.%j.%N"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load muscle/5.1.0

IN_FASTA=$1

# output files for different formats
OUT_FASTA=${IN_FASTA%.*}.MSA.afa
OUT_CLUSTAL=${IN_FASTA%.*}.MSA.aln
OUT_PHYLIP=${IN_FASTA%.*}.MSA.phy

muscle -in $IN_FASTA -threads 32 -fastaout $OUT_FASTA -clwout $OUT_CLUSTAL -phyiout $OUT_PHYLIP

