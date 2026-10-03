#!/usr/bin/env bash

set -euo pipefail

usage() {
  echo "Usage: $0 <output.csv> <BUSCO_DIR> [BUSCO_DIR ...]" >&2
  echo "Each BUSCO directory must contain short_summary.json." >&2
  exit 2
}

if [ "$#" -lt 2 ]; then
  usage
fi

output_csv=$1
shift

if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required to parse BUSCO summaries" >&2
  exit 1
fi

for busco_dir in "$@"; do
  if [ ! -d "$busco_dir" ]; then
    echo "ERROR: BUSCO directory not found: $busco_dir" >&2
    exit 1
  fi
  if [ ! -f "$busco_dir/short_summary.json" ]; then
    echo "ERROR: Missing short_summary.json: $busco_dir/short_summary.json" >&2
    exit 1
  fi
done

{
  echo 'filename,lineage,C%,S%,D%,F%,M%,n'

  for busco_dir in "$@"; do
    summary="$busco_dir/short_summary.json"

    if ! jq -e '.lineage_dataset.name and .results.one_line_summary' "$summary" >/dev/null; then
      echo "ERROR: Missing lineage or result summary in $summary" >&2
      exit 1
    fi

    input_file=$(jq -r '.parameters.in // empty' "$summary")
    if [ -n "$input_file" ]; then
      filename=$(basename "$input_file")
    else
      filename=$(basename "$busco_dir")
    fi

    lineage=$(jq -r '.lineage_dataset.name' "$summary")
    one_line_summary=$(jq -r '.results.one_line_summary' "$summary")

    if [[ $one_line_summary =~ C:([^%]+)%\[S:([^%]+)%,D:([^%]+)%\],F:([^%]+)%,M:([^%]+)%,n:([0-9]+) ]]; then
      complete=${BASH_REMATCH[1]}
      single=${BASH_REMATCH[2]}
      duplicated=${BASH_REMATCH[3]}
      fragmented=${BASH_REMATCH[4]}
      missing=${BASH_REMATCH[5]}
      count=${BASH_REMATCH[6]}
    else
      echo "ERROR: Could not parse BUSCO summary in $summary" >&2
      exit 1
    fi

    jq -nr \
      --arg filename "$filename" \
      --arg lineage "$lineage" \
      --arg complete "$complete" \
      --arg single "$single" \
      --arg duplicated "$duplicated" \
      --arg fragmented "$fragmented" \
      --arg missing "$missing" \
      --arg count "$count" \
      '[$filename, $lineage, $complete, $single, $duplicated, $fragmented, $missing, $count] | join(",")'
  done
} > "$output_csv"