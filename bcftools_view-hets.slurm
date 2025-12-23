#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J bcfv_hets
#SBATCH -o "bcfv.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bcfv.%j.%N.stderr"
#SBATCH -c 17
#SBATCH -p ceres
#SBATCH --mem-per-cpu=16G
#SBATCH -t 1-00:00:00
#SBATCH -N 1

module load bcftools/1.20

VCF=$1

# select heterozygous high quality variants
echo "extracting heterozygous variants with quality > 100 from $VCF"
echo "CMD: bcftools view -i 'QUAL>100 && GT="het"' -Ov -o ${VCF%.*}.hets.vcf $VCF"
bcftools view -i 'QUAL>100 && GT="het"' -Ov -o ${VCF%.*}.hets.vcf $VCF

# write csv counting variant lines per scaffold
echo "counting variants per scaffold"
echo "CMD: grep -a -v ^# ${VCF%.*}.hets.vcf | cut -f1 | sort | uniq -c | awk '{print \$2","\$1}' > ${VCF%.*}.nhets_by_scaffold.csv"
grep -a -v ^# ${VCF%.*}.hets.vcf | cut -f1 | sort | uniq -c | awk '{print $2","$1}' > ${VCF%.*}.nhets_by_scaffold.csv

date
