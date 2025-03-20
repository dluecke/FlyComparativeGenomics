#!/bin/bash
#SBATCH -J windows_vnts
#SBATCH -o "win_var-%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "win_var-%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p short
#SBATCH -N 1

date

VCFFILE=$1
[[ -n $2 ]] && WIN_KB=$2 || WIN_KB=100

echo "Counting variants in ${WIN_KB}kb windows of $VCFFILE"
echo "output: $VCFFILE.nvar_by_windows${WIN_KB}kb.tsv"
echo "command:"
echo "~/FlyComparativeGenomics/get_variants_by_window.sh $VCFFILE $WIN_KB"
~/FlyComparativeGenomics/get_variants_by_window.sh $VCFFILE $WIN_KB
