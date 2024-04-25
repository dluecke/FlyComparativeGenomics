#!/bin/bash

# RunGenomeScope-VPGRU.sh a wrapper to run Kmer counting and GenomeScope2 model
# intended as a QC step for HiFi data, but should work with any .bam or .fq/.fastq
#    (Can extend to .fa in future)
# Designed for USDA Ceres VPGRU project space:
#    prefigured conda environment for R packages 
#    GenomeScope2 install in VPGRU project software directory
#    Ceres modules for jellyfish2 and samtools

usage() { 
    echo "USAGE: $0 [-k|-o|-t|-h] SeqFile.bam"
    echo "  -k INT kmer length, default 21"
    echo "  -o STRING name for run label/output, default SeqFile name (no ext)"
    echo "  -t INT threads, default 32"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args
[ $# -eq 0 ] && usage

# default parameter values
# last arg is the seq file
SEQFILE="${@: -1}"
K_LEN=21
N_CORES=32
SEQNAME=$(basename $SEQFILE)
[[ "SEQNAME" == *".bam" ]] && RUN_ID="${SEQNAME%%.bam}"
[[ "SEQNAME" == *".fq" ]] && RUN_ID="${SEQNAME%%.fq}"
[[ "SEQNAME" == *".fastq" ]] && RUN_ID="${SEQNAME%%.fastq}"

# get options, including call usage if -h flag
while getopts ":hk:o:t:" arg; do
    case $arg in
        k) # kmer size, default 21
            K_LEN=${OPTARG}
            ;;
        o) # name for RunID and output directory, default filename
            RUN_ID="${OPTARG}"
            ;;
        t) # number of threads for SLURM submission
            N_CORES=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done

# call usage if not bam or fastq file
[[ "$SEQFILE" == *".bam" || \
    "$SEQFILE" == *".fastq" || \
    "$SEQFILE" == *".fq"  ]] || { echo "need bam or fastq file"; usage; }

# STARTING SCRIPT ACTIONS

# modules for conda environment (make sure this is good before sbatch)
module unload r # need to use the conda environment R
module load miniconda

# location of Git Repo for GenomeScope
SOFTWARE=/project/vpgru/software
# location of Git Repo for this script (and SLURM template)
GITLOC=$(dirname $0)

# function to enable access to genomescope2 conda environment
conda_env_dir() {
    conda config # writes ~/.condarc file, which will be modified
    echo -e "\nenvs_dirs:\n  /project/vpgru/.conda" >> ~/.condarc # add path to genomescope2 env
    source activate -p /project/vpgru/.conda/genomescope2 || \
        { echo "Problem activating conda env genomescope2" ; exit; }
}

# activate environment genomescope2, using conda_env_dir function to initiate if necessary 
[ -z "$(conda info --envs | grep genomescope2)" ] && \
    conda_env_dir || \
    source activate genomescope2 

# Print some run info to screen
echo -e "\nSubmitting SLURM job KmerAnalysis-${RUN_ID} to perform Kmer analysis on sequence file:\n $SEQFILE"
if [[ "$SEQFILE" == *".bam" ]]; then
    FQ_FILE="${SEQFILE%%.bam}.fastq"
    echo -e "\nThis will first generate FASTQ file:\n $FQ_FILE"
else
    echo -e "\nUsing provided fastq file"
fi
echo -e "Will run jellyfish count and histo using Kmer length ${K_LEN}, with output files:"
echo -e " ${PWD}/${RUN_ID}.jf\n ${PWD}/${RUN_ID}.histo"
echo -e "\nGenomeScope2 analysis will be performed on ${RUN_ID}.histo"
echo -e "using Kmer length ${K_LEN} and ploidy 2, with output files in directory:\n ${PWD}/${RUN_ID}/\n"

# launch slurm template with proper variables
sbatch --job-name="${RUN_ID}-KmerAnalysis" \
    --mail-user="${USER}@usda.gov" \
    -c ${N_CORES} \
    -o "KmerAnalysis-${RUN_ID}.stdout.%j.%N" \
    -e "KmerAnalysis-${RUN_ID}.stderr.%j.%N" \
    --export=ALL,SOFTWARE=${SOFTWARE},SEQFILE=${SEQFILE},KLEN=${K_LEN},RUNID=${RUN_ID},THREADS=${N_CORES} \
    ${GITLOC}/VPGRU-KmerAnalysis_TEMPLATE.slurm

