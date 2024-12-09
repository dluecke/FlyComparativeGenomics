#!/bin/bash
#SBATCH -J windows_nmasked
#SBATCH -o "win_nmask-%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "win_nmask-%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p short
#SBATCH -N 1

date
module load seqtk/1.3

HMASM=$1

echo "Counting number masked bp in 100kb windows of $HMASM using"
which seqtk
echo "output: $HMASM.windows_nmasked.csv"
echo "command:"
echo "~/FlyComparativeGenomics/get_nmasked_by_window.sh $HMASM"
~/FlyComparativeGenomics/get_nmasked_by_window.sh $HMASM
