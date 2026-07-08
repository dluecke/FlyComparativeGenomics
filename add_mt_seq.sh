#!/bin/bash

# Script: add_mt_seq.sh
# Usage: ./add_mt_seq.sh SCAFFOLDS.fa MITO.fa SCAF_MITO_SEQS.txt
# Appends MITO.fa to SCAFFOLDS.fa and removes sequences listed in SCAF_MITO_SEQS.txt
# Output: SCAFFOLDS-w_mt.fa
# Requires: samtools, written with version 1.17

# Prepared 7/8/2026 by David Luecke using Claude Haiku 4.5

usage() {
    cat << EOF
Usage: $(basename "$0") SCAFFOLDS.fa MITO.fa SCAF_MITO_SEQS.txt

Description:
    Appends MITO.fa to SCAFFOLDS.fa and removes sequences listed in SCAF_MITO_SEQS.txt

Arguments:
    SCAFFOLDS.fa           Path to scaffold sequences FASTA file
    MITO.fa                Path to mitochondrial sequence FASTA file
    SCAF_MITO_SEQS.txt    Path to file containing sequence IDs to remove (one per line)

Output:
    {SCAFFOLDS_basename}-w_mt.fa

Requirements:
    samtools

Example:
    $(basename "$0") input_scaffolds.fa mitochondria.fa sequences_to_remove.txt

EOF
    exit 1
}

# Check arguments
if [[ $# -ne 3 ]]; then
    echo "Error: Incorrect number of arguments" >&2
    usage
fi

SCAFFOLDS="$1"
MITO="$2"
REMOVE_LIST="$3"
LOG_FILE="add_mt_seq.log"

# Initialize log file
> "$LOG_FILE"

log_msg() {
    local msg="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $msg" | tee -a "$LOG_FILE"
}

if [[ ! -f "$SCAFFOLDS" ]]; then
    log_msg "Error: SCAFFOLDS file '$SCAFFOLDS' not found"
    exit 1
fi

if [[ ! -f "$MITO" ]]; then
    log_msg "Error: MITO file '$MITO' not found"
    exit 1
fi

if [[ ! -f "$REMOVE_LIST" ]]; then
    log_msg "Error: Remove list file '$REMOVE_LIST' not found"
    exit 1
fi

log_msg "Starting add_mt_seq.sh"
log_msg "Input SCAFFOLDS: $SCAFFOLDS"
log_msg "Input MITO: $MITO"
log_msg "Remove list: $REMOVE_LIST"

# Extract base name without extension for output file (in current working directory)
BASE_NAME="$(basename "${SCAFFOLDS%.*}")"
OUTPUT_FILE="${BASE_NAME}-w_mt.fa"
log_msg "Output file: $OUTPUT_FILE"

# Create temporary combined file
TEMP_FILE=$(mktemp)
log_msg "Creating temporary combined file: $TEMP_FILE"
cat "$SCAFFOLDS" "$MITO" > "$TEMP_FILE"
log_msg "Concatenated $SCAFFOLDS and $MITO"

# Index the combined file
log_msg "Indexing combined file with samtools faidx"
samtools faidx "$TEMP_FILE" 2>&1 | tee -a "$LOG_FILE"

# Create list of sequences to keep
KEEP_LIST=$(mktemp)
log_msg "Creating keep list, filtering sequences from remove list"
samtools faidx "$TEMP_FILE" | cut -f1 | while IFS= read -r seq_id; do
    if ! grep -qx "$seq_id" "$REMOVE_LIST"; then
        echo "$seq_id"
    fi
done > "$KEEP_LIST"
KEEP_COUNT=$(wc -l < "$KEEP_LIST")
log_msg "Sequences to keep: $KEEP_COUNT"

# Extract sequences to keep using samtools faidx
log_msg "Extracting sequences to keep using samtools faidx"
samtools faidx "$TEMP_FILE" $(cat "$KEEP_LIST") > "$OUTPUT_FILE" 2>&1 | tee -a "$LOG_FILE"

log_msg "Cleaning up temporary files"
rm "$TEMP_FILE" "$TEMP_FILE.fai" "$KEEP_LIST"

log_msg "Completed. Output written to $OUTPUT_FILE"
