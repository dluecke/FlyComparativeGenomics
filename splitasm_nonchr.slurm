#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="splitasm_nonchr"
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --time=2:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=4   # 20 processor core(s) per node X 2 threads per core
#SBATCH --mail-user=david.luecke@usda.gov   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output="sn-%j-%N.stdout" # job standard output file (%j replaced by job id)
#SBATCH --error="sn-%j-%N.stderr" # job standard error file (%j replaced by job id)

# USAGE: sbatch splitasm_nonchr.slurm IN_ASM.fa [N_CHR]

# runs splitasm_non_autosomes.sh to write IN_ASM-nonchr_split.fa
# assembly with all but largest N_CHR scaffolds split into contigs

# uses ragtag splitasm
# also requires samtools, seqtk

date

ml ragtag/2.1.0 seqtk/1.3 samtools/1.17

# git repo location
GITLOC=~

IN_ASM=$1

# preserves scaffolding for N_CHR largest scaffolds
# if no $2 provided defaults to 5 (Muscoid standard)
[[ -n $2 ]] && N_CHR=$2 || N_CHR=5

# output to working directory, filename based on input filename
IN_ASM_FN=$(basename $IN_ASM)
OUT_ASM=${IN_ASM_FN%.*}-nonchr_split.fa

echo -e "\nRunning splitasm_non_autosomes.sh on $IN_ASM"
echo "Preserves largest $N_CHR scaffolds"
echo "Splits remaining smaller scaffolds to contigs at N characters"
echo "writes new assembly $OUT_ASM\n"
echo -e "\nragtag splitasm version:"
ragtag.py --version
samtools --version
echo -e "seqtk subseq version:"
seqtk 2>&1 | grep Version

echo -e "\nCMD: $GITLOC/FlyComparativeGenomics/splitasm_non_autosomes.sh $IN_ASM $N_CHR > $OUT_ASM"
$GITLOC/FlyComparativeGenomics/splitasm_non_autosomes.sh $IN_ASM $N_CHR > $OUT_ASM

date