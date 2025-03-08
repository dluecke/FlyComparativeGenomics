#!/bin/bash
#SBATCH -J trimfq
#SBATCH -o "trimfq.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "trimfq.%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p short
#SBATCH -N 1

# trim fastq sequences to 0.001 base error rate

date
module load seqtk/1.3 
which seqtk

FQ_IN=$1
FQ_OUT=$(basename ${FQ_IN%.*}.q001.fq)

echo "trimming sequences in $FQ_IN to quality 0.001"
echo -e "Command:\n seqtk trimfq -q 0.001 $FQ_IN > $FQ_OUT"

seqtk trimfq -q 0.001 $FQ_IN > $FQ_OUT

date