#!/bin/bash

# RunMerylMerqury-VPGRU.sh runs meryl kmer analysis and merqury reads vs assembly tool

usage() { 
    echo "USAGE: $0 [-o|-k|-t|-g|-h] READS.FASTQ ASSEMBLY.FASTA"
    echo "  -o STRING FCS output directory, default ASSEMBLYNAME"
    echo "  -k INT kmer length, default 21"
    echo "  -t INT threads, default 48"
    echo "  -g PATH to FlyComparativeGenomics git repo, default ~/FlyComparativeGenomics"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h" 
[ $# -eq 0 ] && usage

# last 2 args the reference and query seq files
READS_FASTQ="${@: -2:1}"
ASM_FASTA="${@: -1}"
ASM_FN="$(basename $ASM_FASTA)"
OUT_PREFIX="${ASM_FN%%.f*a}"
# default run parameters
K_LEN=21
N_THREAD=48
FCG_PATH=~/FlyComparativeGenomics

# get options, including call usage if -h flag
while getopts ":ho:k:t:g:" arg; do
    case $arg in
        o) # name for RunID and output directory, default filename
            OUT_PREFIX="${OPTARG}"
            ;;
        k) # kmer length for reads meryl
            K_LEN=${OPTARG}
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

[[ -f $READS_FASTQ ]] || { echo "Can't find reads file ${READS_FASTQ}"; usage; }
[[ -f $ASM_FASTA ]] || { echo "Can't find assembly file ${ASM_FASTA}"; usage; }
[[ -d $FCG_PATH ]] || { echo "Can't find FlyComparativeGenomics repo at $FCG_PATH"; usage; }

echo -e "Running meryl and merqury with k=${K_LEN} on reads:\n ${READS_FASTQ}"
echo -e "and assembly:\n ${ASM_FASTA}\nRun in :\n $PWD"
echo -e "with merqury output prefix ${OUT_PREFIX}"
echo "Will also write summaries of taxonomy report in ${OUT_DIR}"

echo -e "\nSubmitting to the mem partition with $N_THREAD tasks"
echo "with job name Meryl-${OUT_PREFIX}"

# launch slurm template with proper variables
sbatch --job-name="Meryl-${OUT_PREFIX}" \
    --mail-user="${USER}@usda.gov" \
    -n ${N_THREAD} \
    -o "Meryl-${OUT_PREFIX}.stdout.%j.%N" \
    -e "Meryl-${OUT_PREFIX}.stderr.%j.%N" \
    --export=ALL,READS_FASTQ=${READS_FASTQ},ASM_FASTA=${ASM_FASTA},K_LEN=${K_LEN}\
MERQURY_OUT=${OUT_PREFIX},THREADS=${N_THREAD},FCG_REPO=${FCG_PATH} \
    ${FCG_PATH}/VPGRU-meryl_merqury_TEMPLATE.slurm

