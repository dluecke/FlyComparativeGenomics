#!/bin/bash

# RunYahs.sh a wrapper to submit Yahs job to cluster with contigs and HiC reads

usage() {
    echo "USAGE: $0 [-o|-g|-h] CTG_ASM.fa[sta] HiC_R1.fq HiC_R2.fq"
    echo "REQUIRED:"
    echo "  CTG_ASM.fa - contig assembly"
    echo "  HiC_R[1,2].fq - HiC paired reads"
    echo "OPTIONAL:"
    echo "  -o STRING prefix for curated scaffold output, default CTG_ASM.yahs"
    echo "      output file STRING_scaffolds_final.fa"
    echo "  -g PATH to directory with git repo, default ~"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h"
[ $# -lt 3 ] && usage
[[ "$*" == *"-h "* ]] && usage

# positional variables
CTG_ASM="${@: -3:1}"
HIC_R1="${@: -2:1}"
HIC_R2="${@: -1}"

# call usage if inputs don't look right
[[ ! -f $CTG_ASM || "$CTG_ASM" != *".fa"* ]] && { echo "no FASTA file"; usage; }
[[ ! -f $HIC_R1 || "$HIC_R1" != *".f"*"q"* ]] && { echo "no HiC R1 file"; usage; }
[[ ! -f $HIC_R2 || "$HIC_R2" != *".f"*"q"* ]] && { echo "no HiC R2 file"; usage; }

# default parameter values
PREFIX=$(basename ${CTG_ASM%.*})".yahs"
GITLOC=~

# get options, including call usage if -h or unrecognized flag
while getopts ":hp:o:t:m:" arg; do
    case $arg in
        o) # output prefix
            PREFIX="${OPTARG}"
            ;;
        g) # location of git repo
            GITLOC=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

echo -e "submitting:\n sbatch $GITLOC/FlyComparativeGenomics/VPGRU-YaHS_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_R1 $HIC_R2"
sbatch $GITLOC/FlyComparativeGenomics/VPGRU-YaHS_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_R1 $HIC_R2