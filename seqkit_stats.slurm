#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="seqkit_stats"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=3G   # replacing mem queue
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 4                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "sks.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "sks.%j.%N.stderr"     #optional, prints our standard error

# USAGE: sbatch seqkit_stats.slurm OUTNAME FAST*gz SEQS.fa FASTQ_FILE2.gz ...
# outputs OUTNAME.tsv with all seqkit stats (see https://bioinf.shenwei.me/seqkit/usage/#stats)

date

module load seqkit/2.4.0

OUTFILE="${1%.tsv}".tsv

echo -e "Calculating statistics using seqkit for files:"
echo ${@:2} | tr ' ' '\n'
echo -e "Output file: $OUTFILE\n"
echo -e "CMD: seqkit stats --all --basename -o $OUTFILE ${@:2}\n"
seqkit stats --all --basename -o "$OUTFILE" ${@:2}

seqkit --help | head
