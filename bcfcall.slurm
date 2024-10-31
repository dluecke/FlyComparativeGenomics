#!/bin/bash
#SBATCH -J bcfcall
#SBATCH -o "bcf.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bcf.%j.%N.stderr"
#SBATCH -c 9
#SBATCH -p mem
#SBATCH -N 1

module load samtools/1.17
module load bcftools/1.20

BAM=$1
ASSEMBLY=$2

samtools index -M $BAM

bcftools mpileup --threads 8 -Ou -f $ASSEMBLY $BAM | bcftools call --threads 8 -cv -Ov -o $BAM.vcf

# write csv counting variant lines per scaffold
grep -v ^# $BAM.vcf | cut -f1 | sort | uniq -c | head | awk '{print $1","$2}' > $BAM.nseg_by_scaffold.csv

