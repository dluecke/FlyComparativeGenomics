#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J windows_nmasked
#SBATCH -o "win_nmask-%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "win_nmask-%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p short
#SBATCH -N 1

date
module load seqtk/1.3

HMASM=$1
[[ -n $2 ]] && WIN_KB=$2 || WIN_KB=100

echo "Counting number masked bp in ${WIN_KB}kb windows of $HMASM using"
which seqtk
echo "output: $HMASM.windows${WIN_KB}kb_nmasked.csv"
echo "command:"
echo "~/FlyComparativeGenomics/get_nmasked_by_window.sh $HMASM $WIN_KB"
~/FlyComparativeGenomics/get_nmasked_by_window.sh $HMASM $WIN_KB
