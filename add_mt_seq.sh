#!/bin/bash

# Script: add_mt_seq.sh
# Usage: ./add_mt_seq.sh SCAFFOLDS.fa MITO.fa SCAF_MITO_SEQS.txt
# Appends MITO.fa to SCAFFOLDS.fa and removes sequences listed in SCAF_MITO_SEQS.txt
# Output: SCAFFOLDS-w_mt.fa
# Requires: samtools, written with version 1.17

# Prepared 7/8/2026 by David Luecke using Claude Haiku 4.5

usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] SCAFFOLDS.fa MITO.fa SCAF_MITO_SEQS.txt

Description:
    Appends MITO.fa to SCAFFOLDS.fa and removes sequences listed in SCAF_MITO_SEQS.txt

Options:
    -h, --help          Show this help message and exit
    -k, --keep-temp     Do not delete temporary files after completion

Arguments:
    SCAFFOLDS.fa           Path to scaffold sequences FASTA file
    MITO.fa                Path to mitochondrial sequence FASTA file
    SCAF_MITO_SEQS.txt    Path to file containing sequence IDs to remove (one per line)

Output:
    ./{SCAFFOLDS_basename}-w_mt.fa

Requirements:
    samtools

Example:
    $(basename "$0") --keep-temp input_scaffolds.fa mitochondria.fa sequences_to_remove.txt

EOF
    exit 1
}

KEEP_TEMP=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -k|--keep-temp)
            KEEP_TEMP=true
            shift
            ;;
        --)
            shift
            break
            ;;
        -*)
            echo "Error: Unknown option '$1'" >&2
            usage
            ;;
        *)
            break
            ;;
    esac
done

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

# Log command invocations explicitly
run_and_log() {
    local cmd="$1"
    log_msg "+ $cmd"
    eval "$cmd"
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

# Create temporary files in current working directory
log_msg "Creating temporary files"
TEMP_FILE=$(mktemp -p . "${BASE_NAME}-tmp.XXXXXX.fa")
log_msg "+ mktemp -p . \"${BASE_NAME}-tmp.XXXXXX.fa\" => $TEMP_FILE"
KEEP_LIST=$(mktemp -p . "${BASE_NAME}-keep.XXXXXX.txt")
log_msg "+ mktemp -p . \"${BASE_NAME}-keep.XXXXXX.txt\" => $KEEP_LIST"
log_msg "Creating temporary combined file: $TEMP_FILE"
run_and_log "cat \"$SCAFFOLDS\" \"$MITO\" > \"$TEMP_FILE\""
log_msg "Concatenated $SCAFFOLDS and $MITO"

# Index the combined file
log_msg "Indexing combined file with samtools faidx"
run_and_log "samtools faidx \"$TEMP_FILE\" 2>&1 | tee -a \"$LOG_FILE\""

# Create list of sequences to keep
log_msg "Creating keep list, filtering sequences from remove list"
# Ensure index exists (samtools faidx writes to ${TEMP_FILE}.fai)
if [[ ! -f "${TEMP_FILE}.fai" ]]; then
    log_msg "Index ${TEMP_FILE}.fai not found; creating index with samtools faidx"
    run_and_log "samtools faidx \"$TEMP_FILE\" 2>&1 | tee -a \"$LOG_FILE\""
fi
# Read sequence IDs from the .fai index file and exclude any IDs listed in remove list
run_and_log "cut -f1 \"${TEMP_FILE}.fai\" | grep -v -F -f \"$REMOVE_LIST\" > \"$KEEP_LIST\""
KEEP_COUNT=$(wc -l < "$KEEP_LIST")
log_msg "Sequences to keep: $KEEP_COUNT"

# Extract sequences to keep using samtools faidx
log_msg "Extracting sequences to keep using samtools faidx"
run_and_log "samtools faidx \"$TEMP_FILE\" $(cat \"$KEEP_LIST\") > \"$OUTPUT_FILE\" 2>&1 | tee -a \"$LOG_FILE\""

log_msg "Completed. Output written to $OUTPUT_FILE"

if [[ "$KEEP_TEMP" == false ]]; then
    log_msg "Cleaning up temporary files"
    run_and_log "rm -f \"$TEMP_FILE\" \"$TEMP_FILE.fai\" \"$KEEP_LIST\""
else
    log_msg "Temporary files retained: $TEMP_FILE, $TEMP_FILE.fai, $KEEP_LIST"
fi
