#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J samtools_depth
#SBATCH -o "sd.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "sd.%j.%N.stderr"
#SBATCH -c 4
#SBATCH -p mem
#SBATCH -N 1
date
module load samtools/1.17

BAM=$1

# check if BAM sorted
SORTCHECK=$(samtools stats --threads 3 $BAM | grep "is sorted:" | awk '{print $NF}') # 0 if unsorted, 1 if sorted
echo "SORTCHECK=$SORTCHECK"
if [[ $SORTCHECK == 0 ]]; then
    echo -e "\nBAM is unsorted, sorting now."
    echo -e "\nCommand:\n samtools sort --threads 3 $BAM > ${BAM%.*}.sort.bam"
    samtools sort --threads 3 $BAM > ${BAM%.*}.sort.bam
    echo -e "\nUsing samtools to write ${BAM}.depth.tsv"
    echo -e "\nCommand:\n samtools depth --threads 3 -a ${BAM%.*}.sort.bam > ${BAM}.depth.tsv"
    samtools depth --threads 3 -a ${BAM%.*}.sort.bam > ${BAM}.depth.tsv
else
    echo -e "\nUsing samtools to write ${BAM}.depth.tsv"
    echo -e "\nCommand:\n samtools depth --threads 3 -a ${BAM} > ${BAM}.depth.tsv"
    samtools depth --threads 3 -a ${BAM} > ${BAM}.depth.tsv
fi

# get average per scaffold, asorti() to sort by scaffold name
awk -v OFS='\t' '{ 
        sum[$1] += $3
        count[$1]++
    } 
    END {
        n=asorti(count, indices); 
        for (i=1; i<=n; i++) {
            print indices[i], count[indices[i]], sum[indices[i]]/count[indices[i]]
        }
    }' ${BAM}.depth.tsv > ${BAM}.depth_by_scaffold.tsv

date
