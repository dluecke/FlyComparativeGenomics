#!/bin/bash
#SBATCH -A vpgru

#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=48   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=ceres    # standard node(s)
#SBATCH --job-name="busco5_inse"
#SBATCH --mail-user=david.luecke@usda.gov   # email address
#SBATCH --mail-type=FAIL
#SBATCH --output="busco5_insecta-%j-%N.stdout" # job standard output file (%j replaced by job id)
#SBATCH --error="busco5_insecta-%j-%N.stderr" # job standard error file (%j replaced by job id)

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load busco5/5.7.1

assembly=$1

busco -m genome -c 48 \
    -i ${assembly} \
    -o ${assembly%.*}_insecta \
    -l /project/vpgru/databases/busco/busco_downloads/lineages/insecta_odb10
