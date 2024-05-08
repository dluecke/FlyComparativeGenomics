#!/bin/bash

# collect_FCS_taxa_summaries.sh generates taxon summary text files from FCS-genome output
# Taxa summaries based on first match (of up to 4) in taxonomy report
# Expects a *.taxonomy.rpt file in SEQID/ results directory
# Produces in current directory:
#  taxon_summary-SEQID.txt lists all reference species and primary divs, sorted by prevalence
#  taxonBroad_summary-SEQID.txt lists all primary divs, sorted by prevalence 
#    input for R script plot_FC-genome_taxa.R which compares taxa in reads vs assembly
#  NonInsectSeqs-SEQID.txt is all taxonomy.rpt lines with only non-insect matches

usage() {
    echo "USAGE: $0 FCS_OUTPUT_DIR/"
    echo " expects file FCS_OUTPUT_DIR/*.taxonomy.rpt"
    echo ' output files named for "FCS_OUTPUT_DIR"'
    exit 0
}

[ -d $1 ] || usage

TAXA_REPORT=$(ls $1/*taxonomy.rpt)
[ -f "$TAXA_REPORT" ] || usage 

SEQ_ID=$(basename $1)

tail -n+3 $TAXA_REPORT | cut -f6,8 | sort -k 2 -t $'\t' | \
    uniq -c > taxonomy_summary-${SEQ_ID}.txt

tail -n+3 $TAXA_REPORT | cut -f8 | sort | \
    uniq -c > taxonomyBroad_summary-${SEQ_ID}.txt

grep : $TAXA_REPORT | grep -v anml:insects > NonInsectSeqs-${SEQ_ID}.txt

