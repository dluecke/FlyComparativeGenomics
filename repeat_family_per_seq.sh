#!/bin/bash
# repeat_family_per_seq.sh
# parses stockholm alignment files to count the number of sequences per repeat family
# prepared by David Luecke 7/31/2026 using GPT-4o
# usage: repeat_family_per_seq.sh <input_file.stk> > output_file.tsv

INPUT_FILE=$1

awk '
BEGIN {
    FS="[ \t:]";  # split on space, tab, or colon
}
# Capture family name from header
/^#=GF ID/ {
    family = $3
    next
}
# Skip comments and end markers
/^#/ || /^\/\// || NF == 0 { next }

{
    seqPrefix = $1  # e.g., Chromosome5
    counts[family, seqPrefix]++
    famSeen[family] = 1
    seqSeen[seqPrefix] = 1
}
END {
    # Print header
    printf "%-25s", "Family"
    for (seq in seqSeen) {
        seqList[++seqCount] = seq
    }
    # Sort sequence IDs alphabetically
    n = asort(seqList)
    for (i = 1; i <= n; i++) {
        printf "\t%s", seqList[i]
    }
    printf "\n"

    # Print counts per family
    for (fam in famSeen) {
        printf "%-25s", fam
        for (i = 1; i <= n; i++) {
            val = counts[fam, seqList[i]]
            if (val == "") val = 0
            printf "\t%d", val
        }
        printf "\n"
    }
}
' $INPUT_FILE
