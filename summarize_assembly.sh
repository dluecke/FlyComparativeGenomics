#!/bin/bash

# summarize_assembly.sh collects assembly summary and qc stats,
# writes a single CSV output for manuscript use
# for now assumes primary assembly

# INPUT: JSON file with following keys:
#   ASM_TAG: STRING for output name, empty string defaults to FILENAME of input FILENAME.json
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
}

if [ "$#" -eq 0 ] || [ ! -f "$1" ]; then
  usage
  exit 1
fi

