#!/bin/bash

# RunGenomeScope-VPGRU.sh a wrapper to run Kmer counting and GenomeScope2 model
# intended as a QC step for HiFi data, but should work with any .bam
#    (Can extend to .fq, or .fa in future)
# Designed for USDA Ceres VPGRU project space:
#    prefigured conda environment for R packages 
#    GenomeScope2 install in VPGRU project software directory
#    Ceres modules for jellyfish2 and samtools

module unload r # need to use the conda environment R
module load miniconda

# location of Git Repo for GenomeScope
SOFTWARE=/project/vpgru/software/

usage() { 
    echo "USAGE: $0 [-k|-o|-h] SeqFile.bam"
    echo "  -k INT kmer length, default 21"
    echo "  -o STRING name for run label/output, default SeqFile name (no ext)"
    echo "  -h FLAG print usage statement"
    exit 0
}

# function to enable access to genomescope2 conda environment
conda_env_dir() {
    conda config
    echo -e "\nenvs_dirs:\n  /project/vpgru/.conda" >> ~/.condarc
    conda activate -p /project/vpgru/.conda/genomescope2 || \
        { echo "Problem activating conda env genomescope2" ; exit; }
}

# last arg is the seq file
SEQFILE="${@: -1}"
echo $SEQFILE

[ $# -eq 0 ] && usage
echo "cleared arg count check" 

[[ "$SEQFILE" =~ "*.bam" ]] || usage

# default parameter values
K_LEN=21
RUN_ID="${SEQFILE%%.bam}"

while getopts ":hk:o:" arg; do
    case $arg in
        k) # kmer size, default 21
            K_LEN=${OPTARG}
            ;;
        o) # name for RunID and output directory, default filename
            RUN_ID="${OPTARG}"
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

echo "Args Set:"
echo "K_LEN: $K_LEN"
echo "RUN_ID: $RUN_ID"
echo "SEQFILE: $SEQFILE"

# activate environment genomescope2 
[ -z "$(conda info --envs | grep genomescope2)" ] && \
    conda_env_dir || \
    conda activate genomescope 

$SOFTWARE/genomescope2.0/genomescope.R 