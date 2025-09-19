#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="fastqc"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p ceres
#SBATCH -t 4:00:00
#SBATCH -o "fqc.%j.%N.stdout"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "fqc.%j.%N.stderr"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load fastqc/0.12.1
fastqc --version
fastqc -t 32 $1
