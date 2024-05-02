#!/bin/bash

# RunMummerAlignment-VPGRU.sh a wrapper to run nucmer whole genome alignment and plot.tsv conversion
# Can take resulting *.plot.tsv and *.breaks.tsv files and plot with alignment_region_dotplot.R
# Requires:
#    convert_gnuplot_to_tsv.sh from git@github.com:dluecke/annotation_tools.git
#    Ceres module mummer/4.0.0rc1

usage() { 
    echo "USAGE: $0 [-c|-o|-p|-t|-g|-h] REFERENCE.fa QUERY.fa"
    echo "  -c INT minimum match length, default=1000"
    echo "  -o STRING alignment name, default REFERENCE_vs_QUERY-cINT"
    echo "  -p STRING slurm partition, default short"
    echo "  -t INT threads, default 8"
    echo "  -g PATH to annotation_tools git repo, default ~/annotation_tools"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args
[ $# -eq 0 ] && usage

# last 2 args the reference and query seq files
REF_FASTA="${@: -2}"
QRY_FASTA="${@: -1}"
# default run parameters
C_VAL=1000
N_CORES=8
PARTITION="short"
ANNOTATIONTOOLSPATH=~/annotation_tools
# trim filenames for default run name
REFFILE=$(basename $REF_FASTA)
QRYFILE=$(basename $QRY_FASTA)
[[ "$REFFILE" == *".fa" ]] && REFNAME="${REFFILE%%.fa}"
[[ "$REFFILE" == *".fasta" ]] && REFNAME="${REFFILE%%.fasta}"
[[ "$QRYFILE" == *".fa" ]] && QRYNAME="${QRYFILE%%.fa}"
[[ "$QRYFILE" == *".fasta" ]] && QRYNAME="${QRYFILE%%.fasta}"
# default run name
RUN_ID="${REFNAME}_vs_${QRYNAME}-c${C_VAL}"

# get options, including call usage if -h flag
while getopts ":hc:o:p:t:g:" arg; do
    case $arg in
        c) # min match length, default 1000
            C_VAL=${OPTARG}
            ;;
        o) # name for RunID and output directory, default filename
            RUN_ID="${OPTARG}"
            ;;
        p) # partition, default short
            PARTITION="${OPTARG}"
            ;;
        t) # number of threads for SLURM submission
            N_CORES=${OPTARG}
            ;;
        g) # path to annotation_tools/ git repo
            ANNOTATIONTOOLSPATH=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

# call usage if not bam or fastq file
[[ "$REFFILE" == *".fasta" || "$REFFILE" == *".fa" ]] \
    && [[ "$QRYFILE" == *".fasta" || "$QRYFILE" == *".fa" ]] \
    || { echo "need fasta file"; usage; }

echo $ANNOTATIONTOOLSPATH
# call usage if no annotation_tools found
[[ -d $ANNOTATIONTOOLSPATH ]] || { echo "can't find annotation_tools/"; usage; }

# STARTING SCRIPT ACTIONS
echo "Script starting for ${RUN_ID}"