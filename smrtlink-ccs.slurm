#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="smrtlink-ccs"
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p ceres
#SBATCH -t 24:00:00
#SBATCH -o "ccs.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "ccs.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date
module load smrtlink/25.3.0.273777

SUBREADS=$1

echo -e "\nRunning SMRTlink CCS on subreads: $SUBREADS\n"

MOVIE=$(echo $SUBREADS | cut -d'.' -f1)

echo -e "\nCMD: ccs $SUBREADS $MOVIE.ccs.bam\n"
ccs $SUBREADS $MOVIE.ccs.bam

echo -e "CCS finished at $(date)\n\nCCS version:"
ccs --version