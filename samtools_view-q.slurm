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

# takes BAM file $1 and runs samtools view to extract all read with map scores above $2
# writes and indexes sorted BAM with headers

# USAGE: sbatch samtools_view-q.slurm IN.BAM [MAPQ default 50]

module load samtools/1.17

IN_BAM=$1
MAPQ=$2
# default map quality 50
[[ -n $MAPQ ]] || MAPQ=50

OUT_BAM=${IN_BAM%.*}-Q${MAPQ}.bam

echo "Filtering reads by mapping quality"
echo "Input BAM: $IN_BAM"
echo "Map quality: $MAPQ"
echo "Output BAM: $OUT_BAM"
date

# samtools options:
# -h keep header to make subsequent samtools view calls easier (eg extract region)
# -u keeps uncompressed until final step (and sets BAM output format)
# -q minimum quality
# -@ threads used in compression
# -c use CSI indexing (BAI index breaks when sequence > ~512Mb)

echo -e "\nCMD: samtools view -h -u -q $MAPQ $IN_BAM | samtools sort -@ 7 -o $OUT_BAM\n"
samtools view -h -u -q $MAPQ $IN_BAM | samtools sort -@ 7 -o $OUT_BAM

echo -e "\nIndexing $OUT_BAM"
echo -e "CMD: samtools index -c $OUT_BAM\n"
samtools index -c $OUT_BAM

samtools --version
date

