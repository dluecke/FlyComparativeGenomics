#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="splitasm_skiplist"
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --time=2:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4   # 20 processor core(s) per node X 2 threads per core
#SBATCH --mail-user=david.luecke@usda.gov   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output="sl-%j-%N.stdout" # job standard output file (%j replaced by job id)
#SBATCH --error="sl-%j-%N.stderr" # job standard error file (%j replaced by job id)

# USAGE: sbatch splitasm_skip_list.slurm IN_ASM.fa SKIP.LIST

# runs splitasm_skip_list.sh to write IN_ASM-nonSKIP_split.fa
# assembly with all except listed scaffolds split into contigs

# uses ragtag splitasm
# also requires samtools, seqtk

date

ml ragtag/2.1.0 seqtk/1.3 samtools/1.17

# git repo location
GITLOC=~

IN_ASM=$1
SKIP_LIST=$2

# output to working directory, filename based on input filename
IN_ASM_FN=$(basename $IN_ASM)
SKIP_LIST_FN=$(basename $SKIP_LIST)
OUT_ASM=${IN_ASM_FN%.*}-non${SKIP_LIST_FN%.*}_split.fa

echo -e "\nRunning splitasm_skip_list.sh on $IN_ASM with skip list $SKIP_LIST"
echo "Preserves scaffolds in $SKIP_LIST"
echo "Splits remaining scaffolds to contigs at N characters with ragtag.py splitasm"
echo "writes new assembly $OUT_ASM\n"
echo -e "\nragtag splitasm version:"
ragtag.py --version
samtools --version
echo -e "seqtk subseq version:"
seqtk 2>&1 | grep Version

echo -e "\nCMD: $GITLOC/FlyComparativeGenomics/splitasm_skip_list.sh $IN_ASM $SKIP_LIST > $OUT_ASM"
$GITLOC/FlyComparativeGenomics/splitasm_skip_list.sh $IN_ASM $SKIP_LIST > $OUT_ASM

date