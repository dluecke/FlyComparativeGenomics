#!/bin/bash

# get n masked bases per window


HMASM=$1
WINDOW=100000

while read SCAF; do 
    LENGTH=$(seqtk subseq $HMASM <(echo $SCAF) | tail -n1 | wc -c)
    N_WINDOWS=$((LENGTH/WINDOW))
    for i in $(seq 0 $N_WINDOWS); do 
        BEG=$((i*WINDOW))
        END0=$(((i+1)*WINDOW-1))
        END=$((END0 < LENGTH ? END0 : LENGTH))
        WINDOW_LENGTH=$((END-BEG+1))
        WINDOW_MASKED=$(seqtk subseq Cmac_v2h.fasta \
            <(echo -e "$SCAF\t$BEG\t$END") | \
            tr -c -d 'N' | wc -c)
        echo $SCAF-$i, $WINDOW_LENGTH, $WINDOW_MASKED
    done
done < <(grep ">" $HMASM | tr -d '>') > $HMASM.windows_nmasked.csv
