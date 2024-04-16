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
SOFTWARE=/project/vpgru/software
# location of Git Repo for this script (and SLURM template)
GITLOC=$(dirname $0)

usage() { 
    echo "USAGE: $0 [-k|-o|-h] SeqFile.bam"
    echo "  -k INT kmer length, default 21"
    echo "  -o STRING name for run label/output, default SeqFile name (no ext)"
    echo "  -h FLAG print usage statement"
    exit 0
}

# function to enable access to genomescope2 conda environment
conda_env_dir() {
    conda config # writes ~/.condarc file, which will be modified
    echo -e "\nenvs_dirs:\n  /project/vpgru/.conda" >> ~/.condarc # add path to genomescope2 env
    source activate -p /project/vpgru/.conda/genomescope2 || \
        { echo "Problem activating conda env genomescope2" ; exit; }
}

# last arg is the seq file
SEQFILE="${@: -1}"
[ $# -eq 0 ] && usage

[[ "$SEQFILE" == *".bam" ]] || { echo "need bam file"; usage; }

# default parameter values
K_LEN=21
SEQNAME=$(basename $SEQFILE)
RUN_ID="${SEQNAME%%.bam}"
FQ_FILE="${SEQFILE%%.bam}.fastq"

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

# activate environment genomescope2, using conda_env_dir function to initiate if necessary 
[ -z "$(conda info --envs | grep genomescope2)" ] && \
    conda_env_dir || \
    source activate genomescope2 

# Print some run info to screen
echo -e "\nSubmitting SLURM job KmerAnalysis-${RUN_ID} to perform Kmer analysis on BAM file:\n $SEQFILE"
echo -e "\nThis will first generate FASTQ file:\n $FQ_FILE"
echo -e "then perform jellyfish count and histo using Kmer length ${K_LEN}, with output files:"
echo -e " ${PWD}/${RUN_ID}.jf\n ${PWD}/${RUN_ID}.histo"
echo -e "\nGenomeScope2 analysis will be performed on ${RUN_ID}.histo"
echo -e "using Kmer length ${K_LEN} and ploidy 2, with output files in directory:\n ${PWD}/${RUN_ID}/\n"

# launch slurm template with proper variables
sbatch --job-name="KmerAnalysis-${RUN_ID}" \
    --mail-user="${USER}@usda.gov" \
    -o "KmerAnalysis-${RUN_ID}.stdout.%j.%N" \
    -e "KmerAnalysis-${RUN_ID}.stderr.%j.%N" \
    --export=ALL,SOFTWARE=${SOFTWARE},BAMFILE=${SEQFILE},KLEN=${K_LEN},RUNID=${RUN_ID} \
    ${GITLOC}/VPGRU-KmerAnalysis_TEMPLATE.slurm

