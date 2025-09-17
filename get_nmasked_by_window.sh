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
        # writing tmp file for both LENGTH and MASKED bc I was paranoid about LENGTH=END-BEG giving +/-1 errors
        # have to do seqtk subseq to get count of Ns anyway 
        seqtk subseq $HMASM <(echo -e "$SCAF\t$BEG\t$END") > tmp.ScafWindow-$HMASM-$SCAF-$BEG-$END.fa
        WINDOW_LENGTH=$(tail -n1 tmp.ScafWindow-$HMASM-$SCAF-$BEG-$END.fa | tr -d '\n' | wc -c)
        WINDOW_MASKED=$(tail -n1 tmp.ScafWindow-$HMASM-$SCAF-$BEG-$END.fa | tr -c -d 'N' | wc -c)
        rm tmp.ScafWindow-$HMASM-$SCAF-$BEG-$END.fa
        echo ${SCAF}:$i, $WINDOW_LENGTH, $WINDOW_MASKED
    done
done < <(grep ">" $HMASM | tr -d '>' | awk '{print $1}') > $HMASM.windows${WIN_KB}kb_nmasked.csv
