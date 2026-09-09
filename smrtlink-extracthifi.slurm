#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="smrtlink-extract"
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p ceres
#SBATCH -t 24:00:00
#SBATCH -o "extract.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "extract.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date
module load smrtlink/26.1.0.284828

CCSREADS=$1

echo -e "\nRunning SMRTlink extracthifi on CCS reads: $CCSREADS"

MOVIE=$(echo $CCSREADS | cut -d'.' -f1)

echo -e "\nCMD: extracthifi $CCSREADS $MOVIE.hifi.bam\n"
extracthifi $CCSREADS $MOVIE.hifi.bam

echo -e "extracthifi finished at $(date)\n\nextracthifi version:"
extracthifi --version