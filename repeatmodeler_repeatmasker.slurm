#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J repeatModMask
#SBATCH -o "rm-stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "rm-stderr.%j.%N"
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -N 1

date
module load repeatmodeler/2.0.4
module load repeatmasker/4.1.5

ASM=$1
DB_NAME=$2

echo "Running repeatmodeler on assembly $1 with output $2"
which RepeatModeler
echo "command:"
echo "BuildDatabase -name $DB_NAME $ASM"
BuildDatabase -name $DB_NAME $ASM

echo "command:"
echo "RepeatModeler -threads 32 -database $DB_NAME"
RepeatModeler -threads 32 -database $DB_NAME
date

echo "Running RepeatMasker on $ASM with library RM_*/consensi.fa.classified"
which RepeatMasker

ln -s RM_*/consensi.fa.classified .
echo "command:"
echo "RepeatMasker -pa 32 -xsmall -dir repmask_RepModLib -e rmblast -lib consensi.fa.classified $ASM"
RepeatMasker -pa 32 -xsmall -dir repmask_RepModLib -e rmblast -lib consensi.fa.classified $ASM

date
