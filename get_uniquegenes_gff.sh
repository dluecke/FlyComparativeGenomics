#!/bin/bash

# get_uniquegenes_gff.sh takes 2 gffs of same assembly A.gff B.gff
# produces gff gene lines for unique regions in each A.uniquegenes.gff B.uniquegenes.gff

# just use bedtools subtract, way faster 

GFF_A=$1
GFF_B=$2

OUT_A="${GFF_A%.gff}.uniquegenes.gff"
OUT_B="${GFF_B%.gff}.uniquegenes.gff"

GFF_A_BN=$(basename $GFF_A)
GFF_B_BN=$(basename $GFF_B)

TMPCOMM_A="comm-$(basename $GFF_A).tmp"
TMPCOMM_B="comm-$(basename $GFF_B).tmp"
TMPAWK="awk-lab-comm-cat${GFF_A_BN%.gff}${GFF_B_BN%.gff}.gff.tsv.tmp"

IFS=$'\n'

while read GFFLOC; do 
    grep "$GFFLOC" $GFF_A
done < <(comm -23 <(awk '$3 == "gene"' $GFF_A | cut -f1-4 | sort) \
                  <(awk '$3 == "gene"' $GFF_B | cut -f1-4 | sort)) \
     > $TMPCOMM_A

while read GFFLOC; do 
    grep "$GFFLOC" $GFF_B
done < <(comm -13 <(awk '$3 == "gene"' $GFF_A | cut -f1-4 | sort) \
                  <(awk '$3 == "gene"' $GFF_B | cut -f1-4 | sort)) \
     > $TMPCOMM_B

awk 'prevBEG>prev2END && prevEND<$5 {print prevLINE} 
    {prev2END=prevEND; prevBEG=$5; prevEND=$6; prevLINE=$0}' \
    <(cat <(awk -v OFS='\t' '{print "A", $0}' $TMPCOMM_A) \
          <(awk -v OFS='\t' '{print "B", $0}' $TMPCOMM_B) | sort -k2,5 ) \
    > $TMPAWK

awk '$1 == "A"' $TMPAWK | cut -f2-10 > $OUT_A
awk '$1 == "B"' $TMPAWK | cut -f2-10 > $OUT_B

unset IFS
rm *.tmp
