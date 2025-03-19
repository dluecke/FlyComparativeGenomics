#!/bin/bash
#SBATCH -J yahs
#SBATCH -o "yahs.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "yahs.%j.%N.stderr"
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -N 1

date

module load apptainer

YAHS_OUT=$1

PCTG_ASM=$2
LEFT_HIC=$3
RIGHT_HIC=$4

# no extension filename of contig assembly
PCTG_NAME=$(basename ${PCTG_ASM%f*a})
# names for HiC alignment sam/bam files
ALN_SAM=$PCTG_NAME"HiC_aln.sam"
ALN_BAM=$PCTG_NAME"HiC_aln.bam"

echo -e "\nStarting bwa alignment of HiC reads:\n $LEFT_HIC\n $RIGHT_HIC"
echo -e "onto contig assembly:\n $PCTG_ASM"
echo -e "commands:\n singularity exec /project/vpgru/software/NextflowSingularityLibrary/staphb-bwa.img bwa index $PCTG_ASM"
echo " singularity exec /project/vpgru/software/NextflowSingularityLibrary/staphb-bwa.img bwa mem -5SP -t 32 $PCTG_ASM $LEFT_HIC $RIGHT_HIC > $ALN_SAM"
singularity exec /project/vpgru/software/NextflowSingularityLibrary/staphb-bwa.img bwa index $PCTG_ASM
singularity exec /project/vpgru/software/NextflowSingularityLibrary/staphb-bwa.img bwa mem -5SP -t 32 $PCTG_ASM $LEFT_HIC $RIGHT_HIC > $ALN_SAM

echo -e "\nConverting alignment from SAM to BAM format"
echo "commands:"
echo " singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools view -S -h -F 2316 $ALN_SAM |\
  singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samblaster |\
  singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools sort -n -@ 32 -o $ALN_BAM"
echo " singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools faidx -o $PCTG_ASM.fai $PCTG_ASM"

singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools view -S -h -F 2316 $ALN_SAM |\
  singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samblaster |\
  singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools sort -n -@ 32 -o $ALN_BAM

singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools faidx -o $PCTG_ASM.fai $PCTG_ASM

echo "BWA steps done"
date

echo -e "\nGenerating stats file for $ALN_BAM"
echo -e "Command:\n singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools stats $ALN_BAM > $ALN_BAM.stats"
singularity exec /project/vpgru/software/NextflowSingularityLibrary/mgibio-alignment_helper-cwl-2.2.1.img samtools stats $ALN_BAM > $ALN_BAM.stats
date

echo -e "\nStarting YaHS scaffolding"
echo -e "command:\n singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-yahs.img yahs -o $YAHS_OUT --no-contig-ec $PCTG_ASM $ALN_BAM"

singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-yahs.img yahs -o $YAHS_OUT --no-contig-ec $PCTG_ASM $ALN_BAM

echo "yahs step finished"
date

echo -e "\nStarting juicer steps to prepare for JBAT curation"
# the in-built yahs juicer (for pre, post) has issues with larger genomes (doesn't calculate scale in same way as juicer_tools) 
# and was dropping ~1/2 of assembly when incorporating review.
# a branch of yahs (https://github.com/zengxiaofei/yahs.git) has a fix for the juicer script, 
# see this thread on yahs issues page: https://github.com/c-zhou/yahs/issues/94
# I added this yahs branch to projects software as yahs-juicer and will use this version for appropriate steps
# also need to use BAM not BIN for this approach, runs fine just slower

echo "commands:"
echo " /project/vpgru/software/yahs-juicer/juicer pre \
  -a -o $YAHS_OUT.JBAT $ALN_BAM ${YAHS_OUT}_scaffolds_final.agp $PCTG_ASM.fai 2> $YAHS_OUT.JBAT/tmp_juicer_pre_JBAT.log"
echo " singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-juicer-tools.img java \
  -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre $YAHS_OUT.JBAT.txt $YAHS_OUT.JBAT.hic.part \
  <(cat $YAHS_OUT.JBAT/tmp_juicer_pre_JBAT.log |\
    grep PRE_C_SIZE | awk '{print $2" "$3}') && \
  mv $YAHS_OUT.JBAT.hic.part $YAHS_OUT.JBAT.hic"

mkdir $YAHS_OUT.JBAT
/project/vpgru/software/yahs-juicer/juicer pre \
  -a -o $YAHS_OUT.JBAT $ALN_BAM ${YAHS_OUT}_scaffolds_final.agp $PCTG_ASM.fai 2> $YAHS_OUT.JBAT/tmp_juicer_pre_JBAT.log

singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-juicer-tools.img java \
  -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre $YAHS_OUT.JBAT.txt $YAHS_OUT.JBAT.hic.part \
  <(cat $YAHS_OUT.JBAT/tmp_juicer_pre_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}') && \
  mv $YAHS_OUT.JBAT.hic.part $YAHS_OUT.JBAT.hic

echo "juicer steps finished"

echo "calculating assembly statistics"
singularity exec /project/vpgru/software/NextflowSingularityLibrary/dmolik-gfastats.img \
   gfastats ${YAHS_OUT}_scaffolds_final.fa > ${YAHS_OUT}_scaffolds_final.fa.stats

echo; date
