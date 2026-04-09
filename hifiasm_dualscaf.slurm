#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="hifiasm_dualscaf"
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -p ceres
#SBATCH -t 3-00:00:00
#SBATCH -o "hfa.%j.%N.stdout"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "hfa.%j.%N.stderr"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# USAGE: sbatch hifiasm_dualscaf.slurm OUTNAME "READS_BLOB*.fq.gz"
#   quote blobs for reads

date
module load hifiasm/0.25.0

OUTNAME=$1
READS=$2

# dual-scaf option tries to scaffold haps against eachother
echo -e "Starting HiFiASM.\nCMD: hifiasm --dual-scaf -o $OUTNAME -t 48 $READS\n"
hifiasm --dual-scaf -o $OUTNAME -t 48 $READS
date

# convert to fasta contigs https://hifiasm.readthedocs.io/en/latest/faq.html#how-do-i-get-contigs-in-fasta
echo "HiFiASM finished. Converting GFA output to fasta."
awk '/^S/{print ">"$2;print $3}' $OUTNAME.bp.p_ctg.gfa > $OUTNAME.bp.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' $OUTNAME.bp.hap1.p_ctg.gfa > $OUTNAME.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' $OUTNAME.bp.hap2.p_ctg.gfa > $OUTNAME.bp.hap2.p_ctg.fa

echo "Finished"
date
