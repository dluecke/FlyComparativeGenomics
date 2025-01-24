#!/bin/bash
#SBATCH -J juicer_post
#SBATCH -o "jp-stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "jp-stderr.%j.%N"
#SBATCH -c 32
#SBATCH -p short
#SBATCH -N 1

date

module load apptainer

# string prefix for output
JBAT_OUT=$1

# input files
# contig assembly
PCTG_ASM=$2 # need a .fasta.fai for this as well
# HiC to assembly BAM
ALN_BAM=$3
# AGP (*_scaffolds_final.agp)
AGP=$4
# .assembly file post JBAT (*.JBAT.review.assembly)
REVIEW=$5

# need newer version of yahs for juicer post tool
/project/vpgru/software/yahs/juicer pre -a -o $JBAT_OUT $ALN_BAM $AGP $PCTG_ASM.fai >$JBAT_OUT.log 2>&1

/project/vpgru/software/yahs/juicer post -o $JBAT_OUT.review $REVIEW $JBAT_OUT.liftover.agp $PCTG_ASM

#singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-juicer-tools.img java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre $JBAT_OUT.txt $JBAT_OUT.hic.part <(cat $JBAT_OUT.log | grep PRE_C_SIZE | awk '{print $2" "$3}') && mv $JBAT_OUT.hic.part $JBAT_OUT.hic

singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img gfastats $JBAT_OUT.review.FINAL.fa > $JBAT_OUT.review.FINAL.fa.stats

date