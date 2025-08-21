#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="st_view-reg"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 1:00:00
#SBATCH --mem-per-cpu=16G   # replacing mem queue
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 8                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "str.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "str.%j.%N.stderr"     #optional, prints our standard error

# takes BAM file and runs samtools view to extract all reads from ALN $1 in region $2

# USAGE: sbatch samtools_view-region.slurm IN.BAM|IN.SAM "seq:c1-c2"

module load samtools/1.17

IN_ALN=$1
REGION=$2
# default map quality 50
[[ -n $REGION ]] || \
  { echo "USAGE: sbatch samtools_view-region.slurm IN.SAM|IN.BAM seq:c1-c2"; exit; }

# convert region syntax seq_N:0-100 to filename syntax seqN_0to100
REG_TAG=$(echo $REGION | tr -d '_' | tr ':' '_' | sed 's/-/to/')

# output SAM filename
OUT_SAM=${IN_ALN%.*}-${REG_TAG}.sam

echo "Extracting reads mapped to region"
echo "Input alignment: $IN_ALN"
echo "Region: $REGION"
echo "Output SAM: $OUT_SAM"
date

# need to index alignment file, use -c to avoid max seq length issue
if [ ! -f ${IN_ALN}.csi ]; then
    echo "Indexing $IN_ALN"
    echo -e "\nCMD: samtools index -c $IN_ALN\n"
    samtools index -c $IN_ALN
fi

# samtools view options:
# -@ threads used in compression
echo "Extracting SAM region"
echo -e "\nCMD: samtools view -@ 7 -o $OUT_SAM $IN_ALN $REGION\n"
samtools view -@ 7 -o $OUT_SAM $IN_ALN $REGION

# extract region sequence (full length, since coordinates aren't modified in SAM) 
# if unambiguous input fasta can be found and sequence file not already present
# array with all possible fastas (any file ending in 'a' best captures .fa and .fasta)
IN_FA_arr=($(echo $IN_ALN | cut -d'.' -f1)*a)
# if single candidate ref fasta is found check for extracted sequence fasta
if [[ ${#IN_FA_arr[@]} == 1 && -f ${IN_FA_arr[0]} ]]; then
    IN_FA=${IN_FA_arr[0]}

    # sequence name without coordinate breaks
    REG_SEQ=$(echo $REGION | cut -d':' -f1)
    REG_FA=${IN_FA%.*}-${REG_SEQ}.fa

    if [ ! -f $REG_FA ]; then
        echo "extracting sequence $REG_SEQ from $IN_FA into $REG_FA"
        echo -e "\nCMD: samtools faidx -o $REG_FA $IN_FA $REGION\n"
        samtools faidx -o $REG_FA $IN_FA $REGION
    fi
fi

samtools --version
date

