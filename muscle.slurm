#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="muscle_align"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p short
#SBATCH -o "ma.stdout.%j.%N"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "ma.stderr.%j.%N"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load muscle/5.1.0

IN_FASTA=$1

# output files for different formats
OUT_FASTA=${IN_FASTA%.*}.MSA.afa
# couldn't get this version of muscle to auto-generate formats, I'll homebrew an interleaved format
OUT_INTERLEAVED=${IN_FASTA%.*}.MSA.iafa

# Running multiple sequence alignment
echo -e "Running muscle aligner on sequences in $IN_FASTA, producing $OUT_FASTA\nUsing:"
muscle --version
echo -e "Command:\n muscle -align $IN_FASTA -threads 32 -output $OUT_FASTA\n"
muscle -align $IN_FASTA -threads 32 -output $OUT_FASTA

# Write interleaved aligned fasta
# line/seq counts for the aligned fasta
N_TOTAL=$(wc -l $OUT_FASTA | awk '{print $1}')
N_SEQS=$(grep -c ">" $OUT_FASTA)
N_PER=$(echo "$N_TOTAL / $N_SEQS" | bc)
# check numbers divide evenly (bc gives float so won't be reversible if not evenly divisible)
if [[ $N_SEQS*$N_PER -ne $N_TOTAL ]]; then 
    echo "Number sequences and number lines in $OUT_FASTA not evenly divisible!"
    exit
fi
# check each sequence is same number of lines
for i in $(seq 1 $N_SEQS); do 
    HEADERLINEi=$(grep -n -m${i} ">" $OUT_FASTA | tail -n1 | cut -d':' -f1)
    if [[ $((1+(i-1)*N_PER)) -ne $HEADERLINEi ]]; then 
        echo "Number of lines per sequence not equal!"
        exit
    fi
done

echo "Writing interleaved aligned fasta file $OUT_INTERLEAVED"
# number characters per aligned sequence line to track alignment coordinates
N_PERLINE=$(awk 'NR==2{print length($1)}' $OUT_FASTA)
# outer loop per interleaved block
for i in $(seq 1 $N_PER); do 
    # inner loop per sequence in block
    for j in $(seq 0 $((N_SEQS-1)) ); do
        # print line i of sequence j
        LINEij=$((i+j*N_PER))
        sed -n "${LINEij}p" $OUT_FASTA
    done
    # starting coordinate for next block (block 1 is headers), skip last
    [[ $i -ne $N_PER ]] && echo $((1+N_PERLINE*(i-1))) || exit 0
done > $OUT_INTERLEAVED
