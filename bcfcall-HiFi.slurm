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

echo "writing BCF file, command:"
echo "CMD: bcftools mpileup --threads 4 -X pacbio-ccs -Ou -f $ASSEMBLY $BAM | bcftools call --threads 4 -cv -Ob -o $BAM.bcf"
bcftools mpileup --threads 4 -X pacbio-ccs -Ou -f $ASSEMBLY $BAM | bcftools call --threads 4 -cv -Ob -o $BAM.bcf
echo "BCF file done"
date

echo -e "\nCalculating variant statistics for $BAM.bcf"
echo "CMD: bcftools stats --threads 4 $BAM.bcf > $BAM.bcf.stats"
bcftools stats --threads 4 $BAM.bcf > $BAM.bcf.stats
echo "bcftools stats done"
date

echo -e "\nConverting BCF to VCF"
echo "CMD: bcftools view --threads 4 -Ov -o $BAM.vcf $BAM.bcf"
bcftools view --threads 4 -Ov -o $BAM.vcf $BAM.bcf
echo "VCF file done"
date

# write csv counting variant lines per scaffold
echo "counting variants per scaffold with $FCG/get_variants_by_sequence.sh"
$FCG/get_variants_by_sequence.sh $BAM.vcf

date
