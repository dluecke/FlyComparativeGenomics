#!/bin/bash
#SBATCH -J samtools_depth
#SBATCH -o "sd.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "sd.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p mem
#SBATCH -N 1
date
module load samtools/1.17

BAM=$1

echo -e "\nUsing samtools to write ${BAM}.depth.tsv"
echo -e "\nCommand:\n samtools depth -a ${BAM} > ${BAM}.depth.tsv"
samtools depth -a ${BAM} > ${BAM}.depth.tsv
# get average per scaffold, asorti() to sort by scaffold name
awk '{ 
        sum[$1] += $3
        count[$1]++
    } 
    END {
        n=asorti(count, indices); 
        for (i=1; i<=n; i++) {
            print indices[i], "\t", sum[indices[i]]/count[indices[i]]
        }
    }' ${BAM}.depth.tsv > ${BAM}.depth_by_scaffold.tsv

date
