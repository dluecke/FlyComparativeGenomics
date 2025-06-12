#!/bin/bash

# RunPurgeDups-VPGRU.sh a wrapper to run purge_dups 
#    on provided primary contigs and PAF HiFi alignment
# Designed for USDA Ceres VPGRU project space:
#    Ceres modules for purge_dups and minimap2

usage() { 
    echo "USAGE: $0 [-o|-m|-t|-p|-h] CTG_ASM.fa[sta] HIFI_ASM_ALIGN.paf[.gz]"
    echo "  -o STRING prefix for purged output, default CTG_ASM-DeDup[_cutoffMVAL]"
    echo "  -m INT depth for manual cutoff, default NULL and auto-cutoff"
    echo "  -t INT threads, default 32"
    echo "  -p STRING slurm partition, default medium"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h"
[ $# -eq 0 ] && usage
[[ "$*" == *"-h "* ]] && usage

# default parameter values
# last arg is the seq file
PRI_ASM="${@: -2:1}"
PAF_IN="${@: -1}"
M_CUT=""
O_SET=""
N_THREAD=32
QUEUE="ceres"
ASM_NAME=$(basename $PRI_ASM)
[[ "$ASM_NAME" == *".fa" ]] && PD_OUT="${ASM_NAME%%.fa}-DeDup"
[[ "$ASM_NAME" == *".fasta" ]] && PD_OUT="${ASM_NAME%%.fasta}-DeDup"

# get options, including call usage if -h flag
while getopts ":hp:o:t:m:" arg; do
    case $arg in
        o) # name for RunID and output directory, default filename
            PD_OUT="${OPTARG}"
            O_SET="set"
            ;;
        m) # cutoff value for manual cutoff, flag to run manual cutoff
            M_CUT=${OPTARG}
            ;;
        t) # number of threads for SLURM submission
            N_THREAD=${OPTARG}
            ;;
        p) # slurm partition
            QUEUE="${OPTARG}"
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

# call usage if not PAF formated alignment file
[[ "$PAF_IN" == *".paf" || \
    "$PAF_IN" == *".paf.gz" ]] || { echo "need alignment in PAF format"; usage; }

# STARTING SCRIPT ACTIONS

# Print info to screen pre-submission
echo -e "Starting purge_dups on assembly:\n ${PRI_ASM}"
echo -e "using HiFi alignment:\n ${PAF_IN}\nRun in :\n $PWD"

if [[ -n $M_CUT ]]; then
    echo "Manual cutoff with -m $M_CUT"
    if [[ -z $O_SET ]]; then
        PD_OUT="${PD_OUT}_cutoff$M_CUT"
    fi
fi

echo -e "Output file prefix:\n ${PD_OUT}"

echo -e "using purge_dups and minimap2 modules\n"
echo -e "Submitting to $QUEUE partition with $N_THREAD threads\n"

# Launch SLURM script via template in Git Repo
# location of Git Repo for this script (and SLURM template)
GITLOC=$(dirname $0)
# launch slurm template with proper variables
sbatch --job-name="purge_dups-${PD_OUT}" \
    --mail-user="${USER}@usda.gov" \
    -p ${QUEUE} \
    -n ${N_THREAD} \
    -o "pd.${PD_OUT}.%j.%N.stdout" \
    -e "pd.${PD_OUT}.%j.%N.stderr" \
    --export=ALL,PRI_ASM=${PRI_ASM},PAF_IN=${PAF_IN},\
PD_OUT=${PD_OUT},N_THREAD=${N_THREAD},QUEUE=${QUEUE},M_CUT=${M_CUT} \
    ${GITLOC}/VPGRU-PurgeDups_TEMPLATE.slurm


