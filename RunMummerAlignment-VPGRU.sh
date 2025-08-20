#!/bin/bash

# RunMummerAlignment-VPGRU.sh a wrapper to run nucmer whole genome alignment and plot.tsv conversion
# Can take resulting *.plot.tsv and *.breaks.tsv files and plot with alignment_region_dotplot.R
# Requires:
#    convert_gnuplot_to_tsv.sh from git@github.com:dluecke/annotation_tools.git
#    Ceres module mummer/4.0.0rc1

usage() { 
    echo "USAGE: $0 [-c|-o|-p|-t|-h] REFERENCE.fa QUERY.fa"
    echo "  -c INT minimum match length, default=2000"
    echo "  -m FLAG run nucmer with --maxmatch, also sets --mem-per-cpu=16G (vs default 2G)"
    echo "  -o STRING alignment name, default REFERENCE-vs-QUERY (c value and maxmatch appended automatically)"
    echo "  -p STRING slurm partition, default ceres"
    echo "  -t INT threads, default 16"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h" 
[ $# -eq 0 ] && usage
[[ "$*" == *" -h"* || $1 == "-h" ]] && usage

# last 2 args the reference and query seq files
REF_FASTA="${@: -2:1}"
QRY_FASTA="${@: -1}"
# default run parameters
C_VAL=2000
MAXMATCH=""
MAXMATCH_TAG=""
N_CORES=16
MEM_PER_CORE="2G"
PARTITION="ceres"
# trim filenames for default run name
REFFILE=$(basename $REF_FASTA)
QRYFILE=$(basename $QRY_FASTA)
[[ "$REFFILE" == *".fa" ]] && REFNAME="${REFFILE%%.fa}"
[[ "$REFFILE" == *".fasta" ]] && REFNAME="${REFFILE%%.fasta}"
[[ "$QRYFILE" == *".fa" ]] && QRYNAME="${QRYFILE%%.fa}"
[[ "$QRYFILE" == *".fasta" ]] && QRYNAME="${QRYFILE%%.fasta}"
# default run name
RUN_ID="${REFNAME}-vs-${QRYNAME}"

# get options, including call usage if -h flag
while getopts ":hmc:o:p:t:" arg; do
    case $arg in
        c) # min match length, default 1000
            C_VAL=${OPTARG}
            ;;
        m) # maxmatch command passed to nucmer
            MAXMATCH="--maxmatch"
            MEM_PER_CORE="16G"
            MAXMATCH_TAG="_maxmatch"
            ;;
        o) # name for RunID and output directory, default filename
            RUN_ID="${OPTARG}"
            ;;
        p) # partition, default ceres
            PARTITION="${OPTARG}"
            ;;
        t) # number of threads for SLURM submission
            N_CORES=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

# call usage if not fasta file
[[ "$REFFILE" == *".fasta" || "$REFFILE" == *".fa" ]] \
    && [[ "$QRYFILE" == *".fasta" || "$QRYFILE" == *".fa" ]] \
    || { echo "need fasta file"; usage; }

# call usage if files aren't found
[[ -f $REF_FASTA && -f $QRY_FASTA ]] || { echo "can't find input files"; usage; }


# STARTING SCRIPT ACTIONS

# alignment name is Run_ID and C value, with maxmatch tag if used
ALIGNMENT_NAME="${RUN_ID}-c${C_VAL}${MAXMATCH_TAG}"

# screen output before SLURM submission
echo "Performing nucmer alignment ${ALIGNMENT_NAME}:"
echo "  Reference: ${REF_FASTA}"
echo "  Query: ${QRY_FASTA}"
echo "  Seed Match Length: ${C_VAL}"
echo "Submitting to SLURM with job name ${RUN_ID}-mummer:"
echo "  Partition:  ${PARTITION}"
echo "  Cores/Threads:  ${N_CORES}"

# location of this script and slurm template
GITLOC=$(dirname $0)

# launch slurm template with proper variables
sbatch --job-name="${RUN_ID}-mummer" \
    --mail-user="${USER}@usda.gov" \
    -p ${PARTITION} \
    -c ${N_CORES} \
    --mem-per-cpu=$MEM_PER_CORE \
    -o "Mummer-${RUN_ID}.stdout.%j.%N" \
    -e "Mummer-${RUN_ID}.stderr.%j.%N" \
    --export=ALL,REFERENCE=${REF_FASTA},QUERY=${QRY_FASTA},\
C_VAL=${C_VAL},ALIGNMENT_NAME=${ALIGNMENT_NAME},\
RUN_ID=${RUN_ID},PARTITION=${PARTITION},\
THREADS=${N_CORES},MAXMATCH=${MAXMATCH} \
    ${GITLOC}/VPGRU-MummerAlignment_TEMPLATE.slurm
