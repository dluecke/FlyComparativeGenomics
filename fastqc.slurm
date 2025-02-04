#!/bin/bash
#SBATCH --job-name="fastqc"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p short
#SBATCH -o "fqc.stdout.%j.%N"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "fqc.stderr.%j.%N"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load fastqc

fastqc $1
