#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J paftools_bam2paf
#SBATCH -o "b2p.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "b2p.%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p ceres
#SBATCH -t 4:00:00
#SBATCH -N 1

# portions of VPGRU-pbmm2_TEMPLATE.slurm triggered by RunPbmm2-VPGRU.sh -p

date

module load samtools/1.17
module load minimap2/2.24

BAM=$1
OUT_PREFIX=${BAM%.*}
[[ -f ${OUT_PREFIX}.bam ]] || \
    (echo "no BAM file"; exit; )


samtools --version
minimap2 --version

echo -e "\nUsing samtools to produce SAM file ${OUT_PREFIX}.sam"
echo -e "\nCommand:\n samtools view -h ${OUT_PREFIX}.bam > ${OUT_PREFIX}.sam\n"
samtools view -h ${OUT_PREFIX}.bam > ${OUT_PREFIX}.sam
echo -e "\n samtools view finished at:"
date

echo -e "\nUsing minimap2 paftools.js to produce ${OUT_PREFIX}.paf.gz"
echo -e "\nCommand:\n paftools.js sam2paf ${OUT_PREFIX}.sam | gzip -c - > ${OUT_PREFIX}.paf.gz"
paftools.js sam2paf ${OUT_PREFIX}.sam | gzip -c - > ${OUT_PREFIX}.paf.gz
echo -e "\n paftools.js sam2paf finshed at:"
date

