#!/bin/bash

# ReciprocalBLAST.sh takes:
#   list of protein query IDs (one per line)
#   target and query protein blast DBs 
# returns:
#   top hit blastp outputs 
#     blastp_out6-query_on_target.tsv
#     blastp_out6-target_on_query.tsv
#   results table with reciprocal blast calls ReciprocalBlastResults.tsv
#   and sequence file with reciprocal hits pairs ReciprocalHitPairs.faa

module load seqtk/1.3
module load blast+/2.15.0

QUERYLIST=$1
TARGETDB=$2
QUERYDB=$3

while read QRY; do 
    seqtk subseq $QUERYDB <(echo $QRY) | \
      blastp -db $TARGETDB -outfmt 6 -evalue 0.1 -max_target_seqs 1 \
      2>/dev/null >> blastp_out6-query_on_target.tsv
    TOPHIT=$(tail -n1 blastp_out6-query_on_target.tsv | cut -f2)
    seqtk subseq $TARGETDB <(echo $TOPHIT) | \
      blastp -db $QUERYDB -outfmt 6 -evalue 0.1 -max_target_seqs 1 \
      2>/dev/null >> blastp_out6-target_on_query.tsv
    RECIP=$(tail -n1 blastp_out6-target_on_query.tsv | cut -f2)
    if [ "$QRY" == "$RECIP" ]; then 
        CHECK=yes
        seqtk subseq $QUERYDB <(echo $QRY) >> ReciprocalHitPairs.faa
        seqtk subseq $TARGETDB <(echo $TOPHIT) >> ReciprocalHitPairs.faa
        echo >> ReciprocalHitPairs.faa
    else 
        CHECK=no
    fi
    echo -e "$QRY\t$TOPHIT\t$RECIP\t$CHECK"
done < $QUERYLIST > ReciprocalBlastResults.tsv

