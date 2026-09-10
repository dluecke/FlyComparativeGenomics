#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="gfastats_by_chr"
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p ceres
#SBATCH -t 1:00:00
#SBATCH -o "gfs.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "gfs.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# USAGE: sbatch gfastats-by_chr.slurm <path/to/ASSEMBLY.fasta> <path/to/chr_list.txt>
# chr_list.txt should be a text file with one chromosome name per line
# writes separate gfastats output files for each chromosome in the list
# all outputs into subdirectory gfastats_by_chr-ASSEMBLY/

date

module load apptainer

ASM=$1
CHR_LIST=$2

ASM_FN=$(basename $ASM)

OUTPUT_DIR="gfastats_by_chr-${ASM_FN%.*}"
mkdir -p $OUTPUT_DIR

while read CHR; do
    singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img gfastats $ASM $CHR> $OUTPUT_DIR/${ASM_FN%.*}-${CHR}.gfastats
done < $CHR_LIST

singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img gfastats -v
