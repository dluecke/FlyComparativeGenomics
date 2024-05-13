#!/bin/bash

# RunFCS-VPGRU.sh runs FCS-genome screen and summarizes taxonomy report

usage() { 
    echo "USAGE: $0 [-o|-t|-g|-h] -x TAXON_ID SEQFILE.[fa|fq|bam]"
    echo " REQUIRED:"
    echo "  -x INT NCBI taxon ID"
    echo "  SEQFILE.[fa|fq|bam] file to screen"
    echo " OPTIONAL:"
    echo "  -o STRING FCS output directory, default SEQFILE"
    echo "  -t INT threads, default 32"
    echo "  -g PATH to FlyComparativeGenomics git repo, default ~/FlyComparativeGenomics"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h" 
[ $# -eq 0 ] && usage

# last 2 args the reference and query seq files
SEQFILE="${@: -1}"
SEQNAME=$(basename $SEQFILE)
OUT_DIR="${SEQNAME%%.[A-Za-z]*}"
# default run parameters
N_CORES=32
FCG_PATH=~/FlyComparativeGenomics

# get options, including call usage if -h flag
while getopts ":hc:o:p:t:g:" arg; do
    case $arg in
        x) # NCBI taxon ID for species
            TAXID=${OPTARG}
            ;;
        o) # name for RunID and output directory, default filename
            OUT_DIR="${OPTARG}"
            ;;
        t) # number of threads for SLURM submission
            N_CORES=${OPTARG}
            ;;
        g) # path to annotation_tools/ git repo
            FCG_PATH=${OPTARG}
            ;;
        h | *) # print help
            echo "Requested help"
            usage
            ;;
    esac
done

[[ -f $SEQFILE ]] || { echo "Can't find input file $SEQFILE"; usage; }
[[ ! -z $TAXID ]] || { echo "No taxon ID provided"; usage; }
[[ -d $FCG_PATH ]] || { echo "Can't find FlyComparativeGenomics repo at $FCG_PATH"; usage; }

[[ -d $OUT_DIR ]] || mkdir -p $OUT_DIR


echo -e "Running FCS-genome screen on file:\n ${SEQFILE}"
echo -e "using taxonomic ID:\n $TAXID"
echo -e "with output directed to:\n $OUTDIR"
echo "Will also write summaries of taxonomy report in $OUTDIR"

echo -e "\nSubmitting to the mem partition with $THREADS CPU"
echo "with job name FCS-${TAXID}-${OUT_DIR}"

# launch slurm template with proper variables
sbatch --job-name="FCS-${TAXID}-${OUT_DIR}" \
    --mail-user="${USER}@usda.gov" \
    -c ${N_CORES} \
    -o "FCS-${TAXID}-${OUT_DIR}.stdout.%j.%N" \
    -e "FCS-${TAXID}-${OUT_DIR}.stderr.%j.%N" \
    --export=ALL,IN_SEQFILE=${SEQFILE},OUTDIR=${OUT_DIR}\
THREADS=${N_CORES},FCG_REPO=${FCG_PATH} \
    ${FCG_PATH}/VPGRU-FCS_screen_TEMPLATE.slurm
