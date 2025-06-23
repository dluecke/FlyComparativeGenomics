#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="MitoHiFi"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p ceres
#SBATCH -t 24:00:00
#SBATCH -o "mhf.%j.%N.stdout"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "mhf.%j.%N.stderr"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# Run MitoHiFi from a singularity image
# Usage: sbatch ~/FlyComparativeGenomics mitohifi.slurm CONTIGS.fa REFMT.fa REFMT.gb
# -o option for code, 5 = invertebrate mitochondrial

module load apptainer

singularity run /project/vpgru/software/OtherSingularityImages/mitohifi_master.sif \
    mitohifi.py \
    -c $1 \
    -f $2 \
    -g $3 \
    -o 5 \
    -t 32
