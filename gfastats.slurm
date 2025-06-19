#!/bin/bash
#SBATCH --job-name="gfastats"
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p ceres
#SBATCH -t 1:00:00
#SBATCH -o "gfastats.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "gfastats.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load apptainer

ASM=$1

singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img gfastats $ASM > $ASM.stats
