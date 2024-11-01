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

date

echo "calling variants for $BAM against $ASSEMBLY with:"
which bcftools
which samtools

echo "indexing $BAM"
samtools index -M $BAM

echo "writing VCF file, command:"
echo "bcftools mpileup --threads 8 -Ou -f $ASSEMBLY $BAM | bcftools call --threads 8 -cv -Ov -o $BAM.vcf"
bcftools mpileup --threads 8 -Ou -f $ASSEMBLY $BAM | bcftools call --threads 8 -cv -Ov -o $BAM.vcf
echo "VCF file done"
date

# write csv counting variant lines per scaffold
echo "counting variants per scaffold, command:"
echo "grep -v ^# $BAM.vcf | cut -f1 | sort | uniq -c | awk '{print $1","$2}' > $BAM.nseg_by_scaffold.csv"
grep -v ^# $BAM.vcf | cut -f1 | sort | uniq -c | awk '{print $1","$2}' > $BAM.nseg_by_scaffold.csv

date
