#!/bin/bash
# GetPurgedBuscoContigs.sh compares Busco output between assemblies 
#   pre- and post- purge_dups, across multiple lineages, 
#   to ID contigs with Busco genes lost in the purge_dups run
# takes 2 directories:
#   $1 with pre-purge_dups Busco results and 
#   $2 with post-purge_dups Busco results
#   REQUIRES SAME LINEAGES in both
# writes for each lineage list of purged Busco benchmarks:
#   BuscosPurged-genes-LINEAGE.txt with Busco IDs 
#   BuscosPurged-contigs-LINEAGE.txt with contig IDs 
#  also list of purged contigs containing lost Buscos from any lineage:
#   BuscosPurged-contigs_union.txt

usage() {
    echo "USAGE: $0 ORIGINAL_ASM_BUSCO/ PURGED_ASM_BUSCO/"
    echo "Directories with Busco runs BEFORE and AFTER contig purging"
    echo "Each Busco lineage must be represented in both"
    exit 0
}

# call usage if no args
[ $# -eq 0 ] && usage

# call usage if both inputs aren't directories
[[ ! -d $1 || ! -d $2 ]] && { echo "inputs not directories"; usage; }

DIR_ORIGINAL=$1
DIR_PURGED=$2

# make arrays files with missing buscos for all lineages in original/purged subdirectories
# NOTE: find output to array only works if no directories have space or special characters
MISSING_ORIGINAL=($(find $DIR_ORIGINAL -type f -name "missing_busco_list.tsv" -print | sort))
MISSING_PURGED=($(find $DIR_PURGED -type f -name "missing_busco_list.tsv" -print | sort))

# get lineage array from each file path, will use to check are making correct comps
LINEAGES=()
for MO in ${MISSING_ORIGINAL[*]}; do
    LINEAGES+=($(echo $MO | sed 's/.*run_\([a-z]*\)_odb[0-9]*.*/\1/'))
done

LINEAGES_PURGED=()
for MP in ${MISSING_PURGED[*]}; do
    LINEAGES_PURGED+=($(echo $MP | sed 's/.*run_\([a-z]*\)_odb[0-9]*.*/\1/'))
done

# check lineage lists the same, otherwise end
if [[ "${LINEAGES[*]}" == "${LINEAGES_PURGED[*]}" ]]; then
    echo -e "Finding purged Buscos from lineages:\n ${LINEAGES[*]}"
else
    echo "Lineages in $DIR_ORIGINAL and $DIR_PURGED do not match"
    echo -e "$DIR_ORIGINAL lineages:\n ${LINEAGES[*]}"
    echo -e "$DIR_PURGED lineages:\n ${LINEAGES_PURGED[*]}"
    usage
fi


