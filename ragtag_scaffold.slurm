#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="rt_scaffold"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 2-00:00:00
#SBATCH --mem-per-cpu=16G   # replacing mem queue
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "rts.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "rts.%j.%N.stderr"     #optional, prints our standard error

# USAGE: sbatch ragtag_scaffold.slurm REFERENCE.fa QUERY.fa

module load ragtag/2.1.0

REF=$1
QRY=$2

# run correct to break contigs if alignment and reads support
# -r to estimate gap size, -u to tag unmodified sequences (better for AGP)
ragtag.py scaffold -r -u -t 32 $REF $QRY
