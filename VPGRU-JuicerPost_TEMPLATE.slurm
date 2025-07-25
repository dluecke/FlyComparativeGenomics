#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J juicer_post
#SBATCH -o "jp-stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "jp-stderr.%j.%N"
#SBATCH -c 32
#SBATCH -p ceres
#SBATCH -t 24:00:00
#SBATCH -N 1

date

module load apptainer

# string prefix for output
JBAT_OUT=$1

# input files
# contig assembly
CTG_ASM=$2 # need a .fasta.fai for this as well
# HiC to assembly BAM
ALN_BAM=$3
# AGP (*_scaffolds_final.agp)
AGP=$4
# .assembly file post JBAT (*.JBAT.review.assembly)
REVIEW=$5

# need branched version of yahs for juicer pre and post tools to avoid mistakes in coordinate scaling with JBAT review
# see YaHS slurm template for more details
echo -e "command:\n /project/vpgru/software/yahs-juicer/juicer pre -a -o $JBAT_OUT $ALN_BAM $AGP $CTG_ASM.fai >$JBAT_OUT.log 2>&1"
/project/vpgru/software/yahs-juicer/juicer pre -a -o $JBAT_OUT $ALN_BAM $AGP $CTG_ASM.fai >$JBAT_OUT.log 2>&1

echo -e "command:\n /project/vpgru/software/yahs-juicer/juicer post -o $JBAT_OUT.review $REVIEW $JBAT_OUT.liftover.agp $CTG_ASM"
/project/vpgru/software/yahs-juicer/juicer post -o $JBAT_OUT.review $REVIEW $JBAT_OUT.liftover.agp $CTG_ASM

#singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-juicer-tools.img java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre $JBAT_OUT.txt $JBAT_OUT.hic.part <(cat $JBAT_OUT.log | grep PRE_C_SIZE | awk '{print $2" "$3}') && mv $JBAT_OUT.hic.part $JBAT_OUT.hic

echo -e "command:\n singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img gfastats $JBAT_OUT.review.FINAL.fa > $JBAT_OUT.review.FINAL.fa.stats"
singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img gfastats $JBAT_OUT.review.FINAL.fa > $JBAT_OUT.review.FINAL.fa.stats

date
