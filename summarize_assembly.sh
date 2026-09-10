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

if [ "$#" -ne 1 ] || [ ! -f "$1" ]; then
  usage
fi

input_json="$1"

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required to parse $input_json" >&2
  usage
fi

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

ASM_TAG=$(jq -r '.ASM_TAG // ""' "$input_json")
GFASTATS_ASM=$(jq -r '.GFASTATS_ASM' "$input_json")
GFASTATS_CHR_DIR=$(jq -r '.GFASTATS_CHR_DIR' "$input_json")
PCT_CHRS=$(jq -r '.PCT_CHRS' "$input_json")
DEPTH_BY_SCAF=$(jq -r '.DEPTH_BY_SCAF' "$input_json")
HIFI_ALN_STATS=$(jq -r '.HIFI_ALN_STATS' "$input_json")
MERQURY_DIR=$(jq -r '.MERQURY_DIR' "$input_json")
mapfile -t BUSCO_DIRS < <(jq -r '.BUSCO_DIRS[]' "$input_json")

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

