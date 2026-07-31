#!/bin/bash
# repeat_family_per_seq.sh
# parses stockholm alignment files to count the number of sequences per repeat family
# prepared by David Luecke 7/31/2026 using GPT-4o, Claude Haiku 4.5, and GNU awk 5.2.1
# usage: repeat_family_per_seq.sh <input_file.stk> > output_file.tsv

INPUT_FILE=$1

awk '
BEGIN {
    # Set field separator to split on spaces, tabs, or colons
    FS="[ \t:]+"
}
# Capture the repeat family name from header lines
/^#=GF ID/ {
    family = $3
    next
}
# Skip comment lines, end markers, or blank lines
/^#/ || /^\/\// || NF == 0 { next }

{
    # Only count sequence rows after a family header has been seen
    if (family == "") next

    seqPrefix = $1
    # Count occurrences of each family/sequence combination
    counts[family, seqPrefix]++
    # Track which families and sequence prefixes have been encountered
    famSeen[family] = 1
    seqSeen[seqPrefix] = 1
}
END {
    # Build sorted list of sequence prefixes
    seqCount = 0
    for (seq in seqSeen) {
        seqCount++
        seqList[seqCount] = seq
    }
    nSeq = asort(seqList)

    # Build sorted list of families
    famCount = 0
    for (fam in famSeen) {
        famCount++
        famList[famCount] = fam
    }
    nFam = asort(famList)

    # Print header row: Family followed by sequence prefixes
    printf "%-25s", "Family"
    for (i = 1; i <= nSeq; i++) {
        printf "\t%s", seqList[i]
    }
    printf "\n"

    # Print one row per family with counts for each sequence prefix
    for (f = 1; f <= nFam; f++) {
        fam = famList[f]
        printf "%-25s", fam
        for (i = 1; i <= nSeq; i++) {
            val = counts[fam, seqList[i]]
            if (val == "") val = 0
            printf "\t%d", val
        }
        printf "\n"
    }
}
' "$INPUT_FILE"
