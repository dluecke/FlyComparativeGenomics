#!/bin/bash

# summarize_assembly.sh collects assembly summary and qc stats,
# writes a single CSV output for manuscript use
# for now assumes primary assembly

# INPUT: JSON file with following keys:
#   ASM_TAG: STRING for output name, empty string defaults to input filename stem
#   GFASTATS_ASM: FILE ASM.gfastats (gfastats output for full assembly)
#   GFASTATS_CHR_DIR: DIRECTORY gfastats_by_chr-ASM/ (output from gfastats-by_chr.slurm)
#   PCT_CHRS: FILE ASM-pct_chrs.txt (output from get_pct_chrs.sh)
#   DEPTH_BY_SCAF: FILE ASM.ALN-depth_by_scaffold.tsv (output from get_depth_by_scaffold.sh on ASM-HiFi depth.tsv)
#   HIFI_ALN_STATS: FILE ASM.ALN.bam.stats (output from samtools stats on ASM-HiFi BAM alignment)
#   MERQURY_DIR: DIRECTORY merqury/ (with output from merqury on ASM with HiFi kmers: ASM*.qv and ASM.completeness.stats)
#   BUSCO_DIRS: [
#       DIRECTORY run_LINEAGE1/ (output from busco on ASM with LINEAGE1 eg run_arthropoda_odb10/),
#       DIRECTORY run_LINEAGE2/ (output from busco on ASM with LINEAGE2 eg run_insecta_odb10/), ...
#     ]

# Prepared 9/10/2026 by David Luecke using MAI-Code-1.1-Flash and GitHub Copilot agent in VS Code

usage() {
  cat <<'EOF'
Usage: summarize_assembly.sh <input.json>

Input JSON keys:
{
  "ASM_TAG": "STRING (optional; empty defaults to input filename stem)",
  "GFASTATS_ASM": "FILE ASM.gfastats (gfastats output for full assembly)",
  "GFASTATS_CHR_DIR": "DIRECTORY gfastats_by_chr-ASM/ (output from gfastats-by_chr.slurm)",
  "PCT_CHRS": "FILE ASM-pct_chrs.txt (output from get_pct_chrs.sh)",
  "DEPTH_BY_SCAF": "FILE ASM.ALN-depth_by_scaffold.tsv (output from get_depth_by_scaffold.sh on ASM-HiFi depth.tsv)",
  "HIFI_ALN_STATS": "FILE ASM.ALN.bam.stats (output from samtools stats on ASM-HiFi BAM alignment)",
  "MERQURY_DIR": "DIRECTORY merqury/ (with output from merqury on ASM with HiFi kmers: ASM*.qv and ASM.completeness.stats)",
  "BUSCO_DIRS": [
    "DIRECTORY run_LINEAGE1/ (output from busco on ASM with LINEAGE1 eg run_arthropoda_odb10/)",
    "DIRECTORY run_LINEAGE2/ (output from busco on ASM with LINEAGE2 eg run_insecta_odb10/)",
    "..."
  ]
}
EOF
  exit 1
}

# empty globs will be used for file checks
shopt -s nullglob

# INITIAL CHECKS FOR CORRECT INPUTS

# check for correct number of arguments
if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  usage
fi

input_json="$1"

# check for jq tool for parsing JSON
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required to parse $input_json" >&2
  usage
fi

# check for required keys in JSON
if ! jq -e '
  type == "object" and
  (.GFASTATS_ASM | type == "string") and
  (.GFASTATS_CHR_DIR | type == "string") and
  (.PCT_CHRS | type == "string") and
  (.DEPTH_BY_SCAF | type == "string") and
  (.HIFI_ALN_STATS | type == "string") and
  (.MERQURY_DIR | type == "string") and
  (.BUSCO_DIRS | type == "array" and all(.[]; type == "string")) and
  ((.ASM_TAG // "") | type == "string")
' "$input_json" >/dev/null 2>&1; then
  usage
fi

# define variables from JSON key fields
ASM_TAG=$(jq -r '.ASM_TAG // ""' "$input_json")
GFASTATS_ASM=$(jq -r '.GFASTATS_ASM' "$input_json")
GFASTATS_CHR_DIR=$(jq -r '.GFASTATS_CHR_DIR' "$input_json")
PCT_CHRS=$(jq -r '.PCT_CHRS' "$input_json")
DEPTH_BY_SCAF=$(jq -r '.DEPTH_BY_SCAF' "$input_json")
HIFI_ALN_STATS=$(jq -r '.HIFI_ALN_STATS' "$input_json")
MERQURY_DIR=$(jq -r '.MERQURY_DIR' "$input_json")
mapfile -t BUSCO_DIRS < <(jq -r '.BUSCO_DIRS[]' "$input_json")

# check that required files and directories exist
for var_name in \
  GFASTATS_ASM \
  PCT_CHRS \
  DEPTH_BY_SCAF \
  HIFI_ALN_STATS; do
  value="${!var_name}"
  if [ ! -f "$value" ]; then
    echo "ERROR: $var_name is not a file: $value" >&2
    usage
  fi
done

for var_name in \
  GFASTATS_CHR_DIR \
  MERQURY_DIR; do
  value="${!var_name}"
  if [ ! -d "$value" ]; then
    echo "ERROR: $var_name is not a directory: $value" >&2
    usage
  fi
done

for i in "${!BUSCO_DIRS[@]}"; do
  value="${BUSCO_DIRS[$i]}"
  if [ ! -d "$value" ]; then
    echo "ERROR: BUSCO_DIRS[$i] is not a directory: $value" >&2
    usage
  fi
done

# Filename arrays from provided directories (missing files return empty array due to nullglob)

# all by-chr gfastats files, default order (check these after writing full assembly summary)
GFASTATS_CHR_FILES=("$GFASTATS_CHR_DIR"/*.gfastats)

# array for all merqury qv files, ordered by wc -l
mapfile -t MERQURY_QV_FILES < <(
    for f in "$MERQURY_DIR"/*.qv; do
        wc -l "$f"
    done | sort -n | awk '{print $2}'
)
# check if MERQURY_QV_FILES array is empty
if [ ${#MERQURY_QV_FILES[@]} -eq 0 ]; then
    echo "ERROR: No .qv files found in $MERQURY_DIR" >&2
    exit 1
fi
if [ ! -f ${MERQURY_QV_FILES[0]} ]; then
    echo "ERROR: First .qv file not found: ${MERQURY_QV_FILES[0]}" >&2
    exit 1
fi

# completeness stats file, only proceed if single file found
if [ "$(compgen -G $MERQURY_DIR/*.completeness.stats | wc -l)" -eq 1 ]; then
    MERQURY_COMPLETENESS="$MERQURY_DIR"/*.completeness.stats
else
    echo "ERROR: Expected single completeness.stats file in $MERQURY_DIR, found $(compgen -G "$MERQURY_DIR"/*.completeness.stats | wc -l)" >&2
    exit 1
fi


# Default output tag: if ASM_TAG is empty, use input filename stem
if [ -z "$ASM_TAG" ]; then
  ASM_TAG=$(basename "${input_json%.*}")
fi

# Output files for full assembly summary and per-chromosome stats
OUTFILE_FULL="${ASM_TAG}-FULL_assembly_summary.csv"
OUTFILE_CHRS="${ASM_TAG}-CHRS_assembly_summary.csv"


# FULL ASSEMBLY SUMMARY

# Extract full assembly stats
# Searches based on gfastats output syntax
TOTAL_BP=$(grep -m1 "Total scaffold length" "$GFASTATS_ASM" | awk '{print $NF}')
MASKED_BP=$(grep -m1 "soft-masked bases" "$GFASTATS_ASM" | awk '{print $NF}')
UNMASKED_BP=$(echo $TOTAL_BP - $MASKED_BP | bc)
N_SCAFFOLDS=$(grep -m1 "scaffolds" "$GFASTATS_ASM" | awk '{print $NF}')
SCAFFOLD_N50=$(grep -m1 "Scaffold N50" "$GFASTATS_ASM" | awk '{print $NF}')
SCAFFOLD_L50=$(grep -m1 "Scaffold L50" "$GFASTATS_ASM" | awk '{print $NF}')
N_CONTIGS=$(grep -m1 "contigs" "$GFASTATS_ASM" | awk '{print $NF}')
CONTIG_N50=$(grep -m1 "Contig N50" "$GFASTATS_ASM" | awk '{print $NF}')
CONTIG_L50=$(grep -m1 "Contig L50" "$GFASTATS_ASM" | awk '{print $NF}')
N_GAPS=$(grep -m1 "gaps in scaffolds" "$GFASTATS_ASM" | awk '{print $NF}')
GAPS_BP=$(grep -m1 "Total gap length" "$GFASTATS_ASM" | awk '{print $NF}')
AVG_GAP_BP=$(grep -m1 "Average gap length" "$GFASTATS_ASM" | awk '{print $NF}')
GC_PCT=$(grep -m1 "GC content" "$GFASTATS_ASM" | awk '{print $NF}')
# pct_chrs.txt stat, convert to rounded percent same as GC_PCT
IN_CHR_DEC=$(grep "PctInChroms" "$PCT_CHRS" | awk '{print $NF}')
IN_CHR_PCT=$(printf "%.2f" $(echo "$IN_CHR_DEC * 100" | bc -l))
# Stats from merqury files
QV_FULL=$(awk '{print $4}' "${MERQURY_QV_FILES[0]}")
QV_COMPLETE=$(awk '{print $5}' "$MERQURY_COMPLETENESS")

# Build BUSCO results keyed by lineage dataset name and value by one-line summary.
# Each BUSCO_DIR must contain a short_summary.json file with:
#   .lineage_dataset.name
#   .results.one_line_summary
# The final associative array length must match BUSCO_DIRS length.
declare -A BUSCO_RESULTS=()
for i in "${!BUSCO_DIRS[@]}"; do
  busco_dir="${BUSCO_DIRS[$i]}"
  short_summary="$busco_dir/short_summary.json"

  if [ ! -f "$short_summary" ]; then
    echo "ERROR: BUSCO_DIRS[$i] is missing short_summary.json: $short_summary" >&2
    exit 1
  fi

  lineage_name=$(jq -er '.lineage_dataset.name' "$short_summary" 2>/dev/null || true)
  one_line_summary=$(jq -er '.results.one_line_summary' "$short_summary" 2>/dev/null || true)

  if [ -z "$lineage_name" ] || [ "$lineage_name" = "null" ] || [ -z "$one_line_summary" ] || [ "$one_line_summary" = "null" ]; then
    echo "ERROR: Missing required BUSCO keys in $short_summary" >&2
    echo "       Expected .lineage_dataset.name and .results.one_line_summary" >&2
    exit 1
  fi

  BUSCO_RESULTS["$lineage_name"]="$one_line_summary"
done

if [ "${#BUSCO_RESULTS[@]}" -ne "${#BUSCO_DIRS[@]}" ]; then
  echo "ERROR: BUSCO_RESULTS length (${#BUSCO_RESULTS[@]}) does not match BUSCO_DIRS length (${#BUSCO_DIRS[@]})" >&2
  exit 1
fi


