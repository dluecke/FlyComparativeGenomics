#!/bin/bash

# get n masked bases per window
# requires seqtk, written with v1.3

HMASM=$1

[[ -n $2 ]] && WIN_KB=$2 || WIN_KB=100

WINDOW=$((WIN_KB*1000))

while read SCAF; do 
    LENGTH=$(seqtk subseq $HMASM <(echo $SCAF) | tail -n1 | wc -c)
    N_WINDOWS=$((LENGTH/WINDOW))
    for i in $(seq 0 $N_WINDOWS); do 
        BEG=$((i*WINDOW))
        END0=$(((i+1)*WINDOW))
        END=$((END0 < LENGTH ? END0 : LENGTH))
        WINDOW_LENGTH=$(seqtk subseq $HMASM \
            <(echo -e "$SCAF\t$BEG\t$END") | \
            tail -n1 | tr -d '\n' | wc -c)
        WINDOW_MASKED=$(seqtk subseq $HMASM \
            <(echo -e "$SCAF\t$BEG\t$END") | \
            tail -n1 | tr -c -d 'N' | wc -c)
        echo $SCAF-$i, $WINDOW_LENGTH, $WINDOW_MASKED
    done
done < <(grep ">" $HMASM | tr -d '>' | awk '{print $1}') > $HMASM.windows${WIN_KB}kb_nmasked.csv
