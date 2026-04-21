#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="seqshannon"
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p ceres
#SBATCH -t 1-00:00:00
#SBATCH -o "sqsh.%j.%N.stdout"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "sqsh.%j.%N.stderr"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# USAGE: sbatch seqshannon.slurm SEQS.fa

PY_ENVS=/project/vpgru/software/python_envs

IN_SEQS=$1
OUT_SHA=${IN_SEQS%.*}-shannon.txt

module load python_3

source $PY_ENVS/SeqShannon/bin/activate

seqshannon -i $IN_SEQS > $OUT_SHA
