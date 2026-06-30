#!/bin/bash

# RunYahs.sh a wrapper to submit Yahs job to cluster with contigs and HiC reads

usage() {
    echo "USAGE: $0 [-o|-s|-g|-h] CTG_ASM.fa[sta] HiC_R1.fq HiC_R2.fq"
    echo "REQUIRED:"
    echo "  CTG_ASM.fa - contig assembly"
    echo "  HiC_R[1,2].fq - HiC paired reads"
    echo "OPTIONAL:"
    echo "  -o STRING prefix for curated scaffold output, default CTG_ASM.yahs"
    echo "      output file STRING_scaffolds_final.fa"
    echo "  -s FLAG split the assembly to contigs before running yahs"
    echo "      modifies -o STRING to CTG_ASM.split.yahs (default) or STRING.split (-o STRING provided)"
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
PREFIX_DEFAULT=true
SPLIT=false
GITLOC=~

# get options, including call usage if -h or unrecognized flag
while getopts ":o:sg:h" arg; do
    case $arg in
        o) # output prefix
            PREFIX_OPT="${OPTARG}"
            PREFIX_DEFAULT=false
            ;;
        s) # split flag
            SPLIT=true
            ;;
        g) # location of git repo
            GITLOC=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

if $PREFIX_DEFAULT; then
    if $SPLIT; then
        PREFIX=$(basename ${CTG_ASM%.*})".split.yahs"
    else
        PREFIX=$(basename ${CTG_ASM%.*})".yahs"
    fi
else
    if $SPLIT; then
        PREFIX=$PREFIX_OPT".split"
    else
        PREFIX=$PREFIX_OPT
    fi
fi

# Submit with differet template depending on if splitting
if $SPLIT; then
    echo -e "submitting:\n sbatch $GITLOC/FlyComparativeGenomics/VPGRU-split_YaHS_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_R1 $HIC_R2"
    sbatch $GITLOC/FlyComparativeGenomics/VPGRU-split_YaHS_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_R1 $HIC_R2
else
    echo -e "submitting:\n sbatch $GITLOC/FlyComparativeGenomics/VPGRU-YaHS_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_R1 $HIC_R2"
    sbatch $GITLOC/FlyComparativeGenomics/VPGRU-YaHS_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_R1 $HIC_R2
fi
