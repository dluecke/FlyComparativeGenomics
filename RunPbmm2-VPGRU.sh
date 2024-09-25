#!/bin/bash

# RunPbmm2-VPGRU.sh runs PacBio tool for HiFi read minimap2 alignment

usage() { 
    echo "USAGE: $0 [-h|-o|-s|-p|-d|-q|-t] ASSEMBLY.FASTA HIFI.FOFN"
    echo "  HIFI.fofn file of file names, bam|fa|fasta|fa.gz|fq|fastq|fq.gz"
    echo "  -o STRING output prefix for BAM file, default ASSEMBLY.HIFI.align"
    echo "  -s FLAG output SAM file"
    echo "  -p FLAG output PAF file (implies -s)"
    echo "  -d FLAG output depth table via samtools depth"
    echo "  -q STRING SLURM submission queue, default mem but check mem768"
    echo "  -t INT threads, default 32"
    echo "  -h FLAG print usage statement"
    exit 0
}

# call usage if no args or "-h" 
[ $# -eq 0 ] && usage
[[ "$*" == *"-h "* ]] && usage

# call usage if not fasta file in assembly postion
[[ "${@: -2:1}" == *".fasta" || \
    "${@: -2:1}" == *".fa"  ]] || { echo "need fasta assembly file"; usage; }
# call usage if not fofn file in reads postion
[[ "${@: -1}" == *".fofn" ]] || { echo "need fofn reads file"; usage; }

# Run Arguments
# last 2 args the reference and query seq files
ASM_FASTA="${@: -2:1}"
HIFI_FOFN="${@: -1}"
ASM_FN="$(basename $ASM_FASTA)"
HIFI_FN="$(basename $HIFI_FOFN)"
# default run parameters
OUT_PREFIX="${ASM_FN%%.f*a}.${HIFI_FN%.fofn}.align"
SLURM_QUEUE="mem"
SAM_OUT="" # empty string won't trigger SAM out section
PAF_OUT="" # empty string won't trigger PAF out section
GET_DEPTH=""
N_THREAD=32

# get options, including call usage if -h flag
while getopts ":hspdo:q:t:" arg; do
    case $arg in
        s) # flag for SAM out
            SAM_OUT="sam"
            ;;
        p) # flag for PAF out
            PAF_OUT="paf"
            ;;
        d) # flag to calculate depths
            GET_DEPTH="depth"
            ;;
        o) # name for RunID and output directory, default filename
            OUT_PREFIX="${OPTARG}"
            ;;
        q) # submission queue
            SLURM_QUEUE="${OPTARG}"
            ;;
        t) # number of threads for SLURM submission
            N_THREAD=${OPTARG}
            ;;
        h | *) # print help
            usage
            ;;
    esac
done


# Ensure conda env pbmm2 activated before calling sbatch
# function to configure access to pbmm2 conda environment (in vpgru projects)
conda_env_dir() {
    #conda config # writes ~/.condarc file, which will be modified
    echo -e "\nenvs_dirs:\n  - /project/vpgru/.conda" >> ~/.condarc # add path to genomescope2 env
    source activate -p /project/vpgru/.conda/pbmm2 || \
        { echo "Problem activating conda env pbmm2" ; exit; }
}

# check for and activate conda environment pbmm2, 
#  using conda_env_dir function to configure and activate if necessary 
module load miniconda
# check conda working (won't if still on head node)
conda info --envs 2> /dev/null || { echo "Conda not working, need to switch to compute node?" ; usage; }

[ -z "$(conda info --envs | grep pbmm2)" ] && \
    conda_env_dir || \
    source activate pbmm2

# Print info to screen pre-submission
echo -e "Starting pbmm2 to run minimap2, aligning reads in:\n ${HIFI_FOFN}"
echo -e "onto assembly:\n ${ASM_FASTA}\nRun in :\n $PWD"
echo -e "with output BAM file:\n ${OUT_PREFIX}.bam"
echo "using pbmm2 conda environment in /project/vpgru/.conda/pbmm2"

if [[ -n $SAM_OUT || -n $PAF_OUT ]]; then
    echo -e "\nUsing samtools to produce SAM file ${OUT_PREFIX}.sam"
fi

if [[ -n $PAF_OUT ]]; then
    echo -e "\nUsing minimap2 paftools.js to produce ${OUT_PREFIX}.paf.gz"
fi

if [[ -n $GET_DEPTH ]]; then
    echo -e "\nUsing samtools depth to write ${OUT_PREFIX}.depth.tsv"
fi

echo -e "\nSubmitting to the ${SLURM_QUEUE} partition with $N_THREAD cores"
echo -e "Job name:\n pbmm2-${OUT_PREFIX}\n"

# Launch SLURM script via template in Git Repo
# location of Git Repo for this script (and SLURM template)
GITLOC=$(dirname $0)
# launch slurm template with proper variables
sbatch --job-name="pbmm2-${OUT_PREFIX}" \
    --mail-user="${USER}@usda.gov" \
    -p ${SLURM_QUEUE} \
    -n ${N_THREAD} \
    -o "pbmm2.stdout-${OUT_PREFIX}.%j.%N" \
    -e "pbmm2.stderr-${OUT_PREFIX}.%j.%N" \
    --export=ALL,ASM_FASTA=${ASM_FASTA},HIFI_FOFN=${HIFI_FOFN},N_THREAD=${N_THREAD},\
OUT_PREFIX=${OUT_PREFIX},SAM_OUT=${SAM_OUT},PAF_OUT=${PAF_OUT},GET_DEPTH=${GET_DEPTH} \
    ${GITLOC}/VPGRU-pbmm2_TEMPLATE.slurm

