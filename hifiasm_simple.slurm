#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="hifiasm_simple"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p ceres
#SBATCH -t 1-00:00:00
#SBATCH -o "hfa.%j.%N.stdout"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "hfa.%j.%N.stderr"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# USAGE: sbatch hifiasm_simple.slurm OUTNAME "READS_BLOB*[bam|fq.gz]"
#   quote blobs for reads

module load hifiasm/0.25.0

OUTNAME=$1
READS=$2

hifiasm -o $OUTNAME -t 32 $READS

# convert to fasta contigs https://hifiasm.readthedocs.io/en/latest/faq.html#how-do-i-get-contigs-in-fasta
awk '/^S/{print ">"$2;print $3}' $OUTNAME.bp.p_ctg.gfa > $OUTNAME.bp.p_ctg.fa
