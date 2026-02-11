#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J samtools_coverage
#SBATCH -o "sc.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "sc.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH --mem-per-cpu=16G
#SBATCH -N 1
date
module load samtools/1.17

BAM=$1

# check if BAM sorted
SORTCHECK=$(samtools stats --threads 3 $BAM | grep "is sorted:" | awk '{print $NF}') # 0 if unsorted, 1 if sorted
echo "SORTCHECK=$SORTCHECK"
if [[ $SORTCHECK == 0 ]]; then
    echo -e "\nBAM is unsorted, sorting now."
    echo -e "\nCMD: samtools sort --threads 3 $BAM > ${BAM%.*}.sort.bam"
    samtools sort --threads 3 $BAM > ${BAM%.*}.sort.bam
    echo -e "\nUsing samtools to write ${BAM}.coverage.tsv"
    echo -e "\nCMD: samtools coverage ${BAM%.*}.sort.bam > ${BAM}.coverage.tsv"
    samtools coverage ${BAM%.*}.sort.bam > ${BAM}.coverage.tsv
else
    echo -e "\nUsing samtools to write ${BAM}.coverage.tsv"
    echo -e "\nCMD: samtools coverage ${BAM} > ${BAM}.coverage.tsv"
    samtools coverage ${BAM} > ${BAM}.coverage.tsv
fi
