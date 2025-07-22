#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="st_view-reg"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=16G   # replacing mem queue
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 8                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "str.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "str.%j.%N.stderr"     #optional, prints our standard error

# takes BAM file and runs samtools view to extract all reads from ALN $1 in region $2

# USAGE: sbatch samtools_view-region.slurm IN.BAM|IN.SAM "seq:c1-c2"

module load samtools/1.17

IN_ALN=$1
REGION=$2
# default map quality 50
[[ -n $REGION ]] || \
  { echo "USAGE: sbatch samtools_view-region.slurm IN.SAM|IN.BAM seq:c1-c2"; exit; }

OUT_SAM=${IN_ALN%.*}-${REGION}.sam

echo "Extracting reads mapped to region"
echo "Input alignment: $IN_ALN"
echo "Region: $REGION"
echo "Output SAM: $OUT_SAM"
date

# samtools view options:
# -@ threads used in compression
echo -e "\nCMD: samtools view -@ 7 -o $OUT_SAM $BAM $REGION\n"
samtools view -@ 7 -o $OUT_SAM $BAM $REGION

samtools --version
date

