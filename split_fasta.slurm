#!/bin/bash

#SBATCH --time=2:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="split_fasta"
#SBATCH --mail-user=david.luecke@usda.gov   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output="std/split_fasta-%j-%N.out" # job standard output file (%j replaced by job id)
#SBATCH --error="std/split_fasta-%j-%N.err" # job standard error file (%j replaced by job id)

# USAGE: sbatch path/to/split_fasta.slurm ASSEMBLY.fasta

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load samtools
module load seqtk

ASSEMBLY=$1

# length threshold for single scaffold fastas
BIG_SCAF=10000000

# build fasta index if not present
if [ ! -f $ASSEMBLY.fai ]; then
	samtools faidx $ASSEMBLY
fi

# count scaffolds longer than threshold
N_BIG_SCAF=$(sort -nr -k2 $ASSEMBLY.fai | awk -v big_scaf="$BIG_SCAF" '$2 > big_scaf' | wc -l)

mkdir split_fastas

# split the N_BIG_SCAF largest scaffolds into individual fastas
while read s; do
	seqtk subseq $ASSEMBLY <(echo $s) > split_fastas/${ASSEMBLY%.f*a}-$s.fasta
done < <(sort -nr -k2 $ASSEMBLY.fai | head -n $N_BIG_SCAF | awk '{print $1}')

# write multifasta with the remaining smaller scaffolds
seqtk subseq $ASSEMBLY <(sort -nr -k2 $ASSEMBLY.fai | tail -n +$((N_BIG_SCAF+1)) | awk '{print $1}') > split_fastas/${ASSEMBLY%.f*a}-SmallScafs.fasta

