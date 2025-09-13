#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J order_chrs
#SBATCH -o "oc.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "oc.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH -t 4:00:00
#SBATCH -N 1

# location for order_chromosomes.sh
FCG=~/FlyComparativeGenomics

# usage: sbatch VPGRU-order_chromosomes.slurm CHR_ASSIGNMENT.tsv
# Will run order_chromosomes.sh on all fastas (*.fa or *.fasta) in current directory using CHR_ASSIGNMENT.tsv

date
module load samtools/1.17

CHR_ASSIGN=$1

# array with fasta file candidates (if no match this will include '*.fa' etc)
FASTA_SCAN=( *.fa *.fasta )

# check each candidate file exists, build FASTAS filename array
unset FASTAS
for f in "${FASTA_SCAN[@]}"; do
    if [ -e "$f" ]; then
        FASTAS+=("$f")
    fi
done

echo -e "\n================================================"
echo -e "\nRE-ORDERING AND ORIENTING CHROMOSOMES"
echo -e "\nScaffold assemblies:"
printf '  %s\n' "${FASTAS[@]}"
echo -e "\nChromosome assignments:"
awk '{print "  " $0}' $CHR_ASSIGN

for f in ${FASTAS[@]}; do
    echo -e "\n------------------------------------------------\n"
    echo "Starting chromosome ordering for scaffolds $f using assignments $CHR_ASSIGN"
    $FCG/order_chromosomes.sh $f $CHR_ASSIGN
done

echo -e "\n================================================\n"
samtools --version
