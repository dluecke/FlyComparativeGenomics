#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="smrtlink-pbi"
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p ceres
#SBATCH -t 24:00:00
#SBATCH -o "pbi.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "pbi.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date
module load smrtlink/26.1.0.284828

PB_BAM=$1

echo -e "\nRunning SMRTlink pbindex on $PB_BAM"
echo -e "CMD: pbindex $PB_BAM\n"
pbindex $PB_BAM

echo -e "Extracting index information with pbindexdump"
echo "CMD: pbindexdump ${PB_BAM}.pbi > ${PB_BAM}.dump.json"
pbindexdump ${PB_BAM}.pbi > ${PB_BAM}.dump.json

date
echo -e "\nFinished. versions:"
pbindex --version | head -n1
pbindexdump --version