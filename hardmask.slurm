#!/bin/bash
#SBATCH -A vpgru
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
# Also useful to get number of Ns in scaffold gaps for any assembly
#   (ignore "-hm" tag on nmask.csv file, "hm.fa" will be copy of unmasked)

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

# if line doesn't start with > replace any of a, c, t, g, or n with N
sed '/^>/!s/[actgn]/N/g' $ASM > $HMASM

echo "writing CSV with scaffold lengths and number masked bp: $HMASM.nmasked.csv"

while read SCAF; do 
    # writing tmp file for both LENGTH and MASKED bc I was paranoid about +/-1 errors
    # have to do seqtk subseq to get count of Ns anyway 
    seqtk subseq $HMASM <(echo $SCAF) > tmp.Scaf-$HMASM-$SCAF.fa
    LENGTH=$(tail -n1 tmp.Scaf-$HMASM-$SCAF.fa | tr -d '\n' | wc -c)
    MASKED=$(tail -n1 tmp.Scaf-$HMASM-$SCAF.fa | tr -c -d 'N' | wc -c)
    rm tmp.Scaf-$HMASM-$SCAF.fa
    echo $SCAF,$LENGTH,$MASKED
done < <(grep ">" $HMASM | awk '{print $1}' | tr -d '>') > $HMASM.nmasked.csv

date
