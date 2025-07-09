#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J paftools_sam2paf
#SBATCH -o "s2p.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "s2p.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p ceres
#SBATCH -t 4:00:00
#SBATCH -N 1

# final step in portion of VPGRU-pbmm2_TEMPLATE.slurm triggered by RunPbmm2-VPGRU.sh -p

date

module load minimap2/2.24

SAM=$1
OUT_PREFIX=${SAM%.*}
[[ -f ${OUT_PREFIX}.sam ]] || \
  { echo "no SAM file"; exit; }


minimap2 --version

echo -e "\nUsing minimap2 paftools.js to produce ${OUT_PREFIX}.paf.gz"
echo -e "\nCommand:\n paftools.js sam2paf ${OUT_PREFIX}.sam | gzip -c - > ${OUT_PREFIX}.paf.gz"
paftools.js sam2paf ${OUT_PREFIX}.sam | gzip -c - > ${OUT_PREFIX}.paf.gz
echo -e "\n paftools.js sam2paf finshed at:"
date

