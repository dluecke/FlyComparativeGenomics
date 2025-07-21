#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="st_view-q"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 4:00:00
#SBATCH --mem-per-cpu=16G   # replacing mem queue
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 8                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "stq.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "stq.%j.%N.stderr"     #optional, prints our standard error

# takes BAM file and runs samtools view to extract all reads from BAM $1 with map scores above $2

# USAGE: sbatch samtools_view-q.slurm IN.BAM [MAPQ default 50]

module load samtools/1.17

BAM=$1
MAPQ=$2
# default map quality 50
[[ -n $MAPQ ]] || MAPQ=50

OUT_SAM=${BAM%.*}-Q${MAPQ}.sam

echo "Filtering reads by mapping quality"
echo "Input BAM: $BAM"
echo "Map quality: $MAPQ"
echo "Output SAM: $OUT_SAM"
date

# samtools view options:
# -h keep header to make subsequent samtools view calls easier (eg extract region)
# -@ threads used in compression
echo -e "\nCMD: samtools view -h -@ 7 -o $OUT_SAM $BAM\n"
samtools view -h -@ 7 -o $OUT_SAM $BAM

samtools --version
date

