#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="ULI-lima_pbmarkdups"
#SBATCH -N 1
#SBATCH -n 17
#SBATCH -p ceres
#SBATCH -t 6:00:00
#SBATCH -o "ULI.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "ULI.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load smrtlink/25.1.0.257715 
echo "smrtlink loaded"
which lima
which pbmarkdup

# input BAM, need to trim ULI adaptors
RAW_BAM=$1

# ULI adaptor sequence
echo -e ">ULI_adaptor\nAAGCAGTGGTATCAACGCAGAGTACT" > ULI_adaptor.fa

# intermediate/output file names
TMP_BAM=$(echo $RAW_BAM | sed 's/.bam/.tmptrim.bam/')
OUT_BAM=$(echo $RAW_BAM | sed 's/.bam/.ULItrim_markdups.bam/')

# run lima to trim adaptors
echo "running lima to trim ULI adaptors, command:"
echo " lima --num-threads 16 $RAW_BAM ULI_adaptor.fa $TMP_BAM"
lima --num-threads 16 $RAW_BAM ULI_adaptor.fa $TMP_BAM

# run pbmarkdups to remove PCR duplicates
echo "running pbmarkdups to remove PCR duplicates, command:"
echo " pbmarkdup $TMP_BAM $OUT_BAM"
pbmarkdup $TMP_BAM $OUT_BAM 

# rm $TMP_BAM
