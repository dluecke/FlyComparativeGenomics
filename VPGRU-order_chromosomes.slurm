#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J order_chrs
#SBATCH -o "oc.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "oc.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH -t 2:00:00
#SBATCH -N 1

# location for order_chromosomes.sh
FCG=~/FlyComparativeGenomics

# usage: sbatch VPGRU-order_chromosomes.slurm CHR_ASSIGNMENT.tsv
# Will run order_chromosomes.sh on all fastas (*.fa or *.fasta) in current directory using CHR_ASSIGNMENT.tsv

date
module load samtools/1.17

CHRS=$1

# array with fastas (if none will include '*.fa' etc)
FASTAS=( *.fa *.fasta )

for f in ${FASTAS[@]}; do
    # check that file exists to avoid running with '*.fa' or '*.fasta'
    if [ -e "$f" ]; then
        echo -e "\n_________________________________________________\n"
        echo "Starting chromosome ordering for scaffolds $f using assignments $1"
        $FCG/order_chromosomes.sh $f $1
    fi
done

echo -e "\n------------------------------------------------\n"
samtools --version
