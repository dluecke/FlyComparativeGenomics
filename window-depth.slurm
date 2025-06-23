#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J windows_depth
#SBATCH -o "win_depth-%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "win_depth-%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p short
#SBATCH -N 1

date

DEPTHFILE=$1
[[ -n $2 ]] && WIN_KB=$2 || WIN_KB=100

echo "Counting average depth in ${WIN_KB}kb windows of $DEPTHFILE"
echo "output: $(echo $DEPTHFILE | sed "s/depth/depth_by_windows${WIN_KB}kb/")"
echo "command:"
echo "~/FlyComparativeGenomics/get_depth_by_window.sh $DEPTHFILE $WIN_KB"
~/FlyComparativeGenomics/get_depth_by_window.sh $DEPTHFILE $WIN_KB
