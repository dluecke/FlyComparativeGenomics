#!/bin/bash

# RunFCS-VPGRU.sh runs FCS-genome screen and summarizes taxonomy report

usage() { 
    echo "USAGE: $0 [-o|-t|-g|-h] -x TAXON_ID PATH/TO/SEQFILE.[fa|fq|bam]"
    echo " REQUIRED:"
    echo "  -x INT NCBI taxon ID"
    echo "  SEQFILE.[fa|fq|bam] file to screen"
    echo " OPTIONAL:"
    echo "  -o STRING FCS output directory, default SEQFILENAME-FCSout"
    echo "  -t INT threads, default 32"
    echo "  -g PATH to FlyComparativeGenomics git repo, default ~/FlyComparativeGenomics"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h" 
[ $# -eq 0 ] && usage

# last 2 args the reference and query seq files
SEQFILE="${@: -1}"
#OUT_DIR="$(basename $SEQFILE)-FCSout"
OUT_DIR="FCS-${SEQFILE%.*}"
# default run parameters
CLEAN="c"
N_THREAD=32
FCG_PATH=~/FlyComparativeGenomics

# get options, including call usage if -h flag
while getopts ":hx:o:t:g:" arg; do
    case $arg in
        x) # NCBI taxon ID for species
            TAXID=${OPTARG}
            ;;
        o) # name for RunID and output directory, default filename
            OUT_DIR="${OPTARG}"
            ;;
        t) # number of threads for SLURM submission
            N_THREAD=${OPTARG}
            ;;
        g) # path to annotation_tools/ git repo
            FCG_PATH=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

[[ -f $SEQFILE ]] || { echo "Can't find input file $SEQFILE"; usage; }
[[ ! -z $TAXID ]] || { echo "No taxon ID provided"; usage; }
[[ -d $FCG_PATH ]] || { echo "Can't find FlyComparativeGenomics repo at $FCG_PATH"; usage; }

if [[ "$SEQFILE" == *".bam" || \
      "$SEQFILE" == *".fq" || "$SEQFILE" == *".fastq" || \
      "$SEQFILE" == *".fq.gz" || "$SEQFILE" == *".fastq.gz" || \
      "$SEQFILE" == *".fa" || "$SEQFILE" == *".fasta" || \
      "$SEQFILE" == *".fa.gz" || "$SEQFILE" == *".fasta.gz" ]]; then
    echo "$SEQFILE extension recognized"
else
    echo "file extension not recognized: $SEQFILE"
    exit;
fi

[[ -d $OUT_DIR ]] || mkdir -p $OUT_DIR


echo -e "\nRunning FCS-genome screen on file:\n ${SEQFILE}"
echo -e "using taxonomic ID:\n $TAXID"
echo -e "with output directed to:\n $OUT_DIR"
echo "Will also write summaries of taxonomy report in $OUT_DIR"

if [[ -n $CLEAN ]]; then
    echo -e "\nRunning FCS clean genome with action report to be generated in $OUT_DIR"
fi

echo -e "\nSubmitting to the mem partition with $N_THREAD tasks"
echo "with job name FCS-${TAXID}-${OUT_DIR}"

# launch slurm template with proper variables
sbatch --job-name="FCS-${TAXID}-${OUT_DIR}" \
    --mail-user="${USER}@usda.gov" \
    -n ${N_THREAD} \
    -o "FCS-${TAXID}-${OUT_DIR}.stdout.%j.%N" \
    -e "FCS-${TAXID}-${OUT_DIR}.stderr.%j.%N" \
    --export=ALL,IN_SEQFILE=${SEQFILE},OUTDIR=${OUT_DIR},\
TAXID=${TAXID},CLEAN=${CLEAN},THREADS=${N_THREAD},FCG_REPO=${FCG_PATH} \
    ${FCG_PATH}/VPGRU-FCS_screen_TEMPLATE.slurm
