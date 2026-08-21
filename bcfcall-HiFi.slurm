#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J bcfcall-HiFi
#SBATCH -o "bcfHF.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bcfHF.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH --mem-per-cpu=16G
#SBATCH -t 3-00:00:00
#SBATCH -N 1

# variant calling for HiFi reads, using bcftools mpileup and call

FCG=~/FlyComparativeGenomics
module load samtools/1.17
module load bcftools/1.20

BAM=$1
ASSEMBLY=$2

date

echo "calling variants for $BAM against $ASSEMBLY with:"
which bcftools
which samtools

echo "indexing $BAM"
samtools index -M -c $BAM

echo "writing VCF file, command:"
echo "bcftools mpileup --threads 4 -X pacbio-ccs -Ou -f $ASSEMBLY $BAM | bcftools call --threads 4 -cv -Ov -o $BAM.vcf"
bcftools mpileup --threads 4 -X pacbio-ccs -Ou -f $ASSEMBLY $BAM | bcftools call --threads 4 -cv -Ov -o $BAM.vcf
echo "VCF file done"
date

# write csv counting variant lines per scaffold
echo "counting variants per scaffold with $FCG/get_variants_by_sequence.sh"
$FCG/get_variants_by_sequence.sh $BAM.vcf

#echo "grep -a -v ^# $BAM.vcf | cut -f1 | sort | uniq -c | awk '{print \$2","\$1}' > $BAM.nseg_by_scaffold.csv"
#grep -a -v ^# $BAM.vcf | cut -f1 | sort | uniq -c | awk '{print $2","$1}' > $BAM.nseg_by_scaffold.csv

date
