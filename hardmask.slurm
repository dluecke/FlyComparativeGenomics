#!/bin/bash
#SBATCH -J hardmask
#SBATCH -o "hm.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "hm.%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p short
#SBATCH -N 1

# USAGE: sbatch path/to/this.slurm assembly.masked.fa 
# takes softmasked assembly assembly.masked.fa
# writes hardmasked assembly assembly.hardmasked.fa
# writes assembly.hardmasked.fa.nmasked.csv file with
#    total length and number masked bases per scaffold

date 

module load seqtk
which seqtk

ASM=$1
HMASM=$(echo $ASM | sed 's/masked/hardmasked/')

# make sure rename worked, otherwise will delete $ASM file
if [[ "$ASM" == "$HMASM" ]]; then
    HMASM="${ASM%.*}-hm.fa"
fi

echo "writing hardmasked verison of $ASM: $HMASM"

tr '[actg]' 'N' < $ASM | sed 's/sNNffold/scaffold/' > $HMASM

echo "writing CSV with scaffold lengths and number masked bp: $HMASM.nmasked.csv"

while read SCAF; do 
    LENGTH=$(seqtk subseq $HMASM <(echo $SCAF) | tail -n1 | tr -d '\n' | wc -c)
    MASKED=$(seqtk subseq $HMASM <(echo $SCAF) | tr -c -d 'N' | wc -c)
    echo $SCAF,$LENGTH,$MASKED
done < <(grep ">" $HMASM | tr -d '>') > $HMASM.nmasked.csv

date
