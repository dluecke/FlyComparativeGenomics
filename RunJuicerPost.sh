#!/bin/bash

# RunJuicerPost.sh a wrapper to incorporate Juicebox Assembly Tools (JBAT) curation
# Submit juicer post job to cluster

usage() {
    echo "USAGE: $0 [-o|-g|-h] CTG_ASM.fa[sta] ASM_HIC.bam YAHS_OUT.agp REVIEW.assembly"
    echo "REQUIRED:"
    echo "  CTG_ASM.fa - contig assembly (input to yahs)"
    echo "  ASM_HIC.bam - BAM alignment of HiC to assembly (used by yahs)"
    echo "  YAHS_OUT.agp - AGP output from yahs (eg *scaffolds_final.agp)"
    echo "  REVIEW.assembly - ASSEMBLY file output from JBAT curation"
    echo "OPTIONAL:"
    echo "  -o STRING prefix for curated scaffold output, default CTG_ASM.yahs"
    echo "      output file STRING.JBAT.final.fa"
    echo "  -g PATH to directory with git repo, default ~"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h"
[ $# -lt 4 ] && usage
[[ "$*" == *"-h "* ]] && usage

# positional variables
CTG_ASM="${@: -4:1}"
HIC_BAM="${@: -3:1}"
AGP="${@: -2:1}"
REVIEW="${@: -1}"

# call usage if inputs don't look right
[[ ! -f $CTG_ASM || "$CTG_ASM" != *".fa"* ]] && { echo "no FASTA file"; usage; }
[[ ! -f $HIC_BAM || "$HIC_BAM" != *".bam" ]] && { echo "no BAM file"; usage; }
[[ ! -f $AGP || "$AGP" != *".agp" ]] && { echo "no AGP file"; usage; }
[[ ! -f $REVIEW || "$REVIEW" != *".assembly" ]] && { echo "no ASSEMBLY file"; usage; }

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

echo -e "submitting:\n sbatch $GITLOC/FlyComparativeGenomics/VPGRU-JuicerPost_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_BAM $AGP $REVIEW"
sbatch $GITLOC/FlyComparativeGenomics/VPGRU-JuicerPost_TEMPLATE.slurm $PREFIX $CTG_ASM $HIC_BAM $AGP $REVIEW