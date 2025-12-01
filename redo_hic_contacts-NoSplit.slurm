#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J redo_hic
#SBATCH -o "hic.%j.%N.stdout"  # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "hic.%j.%N.stderr"
#SBATCH -c 16
#SBATCH --mem-per-cpu=16G   # replacing mem queue
#SBATCH -p ceres
#SBATCH -t 36:00:00
#SBATCH -N 1
#SBATCH -C INTEL

# USAGE: sbatch path/to/redo_hic_contacts.slurm ASSEMBLY.fa HIC_R1.fq HIC_R2.fq

# script for regenerating clean HiC contact map on while scaffold assembly (no contig boundaries)
# used to distinguish JBAT artifacts from misassemblies, fine tune curation
# Steps:
# index scaffold assembly
#   using samtools faidx
# map HiC reads to assembly, filter based on HiC methods (same as OTB)
#   using ceres modules bwa_mem2 and samtools, 
#   and samblaster in $SOFTWARE_DIR
# generate .agp and .assembly from fasta
#   using juicebox_scripts (https://github.com/phasegenomics/juicebox_scripts)
#   installed in $SOFTWARE_DIR
# generate paired contacts JBAT.txt and output for chrom.sizes from BAM and fasta index
#   using yahs juicer pre from branch with scaling fix (https://github.com/zengxiaofei/yahs.git)
# generate .hic file from outputs of yahs juicer pre step
#   using singularity image of juicer_tools in $SOFTWARE_DIR/NextflowSingularityLibrary (same as OTB)

# directory with any non-module installs and singularity images
SOFTWARE_DIR=/project/vpgru/software

module load bwa_mem2/2.2.1
module load samtools/1.17
module load apptainer

SCAF_ASM=$1
HIC_R1=$2
HIC_R2=$3

NCORE=16

SCAF_ASM_FN=$(basename $SCAF_ASM)
HIC_BAM=${SCAF_ASM_FN%.*}.HiC.bam

# link for target scaffold assembly
[[ -f $SCAF_ASM_FN ]] !! ln -s $SCAF_ASM .

date

echo -e "\nIndexing contig assembly with samtools:"
which samtools
echo -e "\ncommand:\n samtools faidx $SCAF_ASM_FN"
samtools faidx $SCAF_ASM_FN

echo -e "\nMapping HiC reads $HIC_R1 and $HIC_R2 onto $SCAF_ASM_FN using bwa-mem2:"
which bwa-mem2
echo -e "and filtering using samtools and samblaster:"
which samtools
$SOFTWARE_DIR/samblaster/samblaster --version
echo -e "\nOutput file: $HIC_BAM"
echo -e "\ncommands:\n bwa-mem2 index -p ${SCAF_ASM_FN%.*} $SCAF_ASM_FN"
echo -e " bwa-mem2 mem -5SP -t $NCORE ${SCAF_ASM_FN%.*} $HIC_R1 $HIC_R2 |\
 samtools view -h -F 2316 |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $NCORE -o $HIC_BAM"
bwa-mem2 index -p ${SCAF_ASM_FN%.*} $SCAF_ASM_FN
bwa-mem2 mem -5SP -t $NCORE ${SCAF_ASM_FN%.*} $HIC_R1 $HIC_R2 |\
 samtools view -h -F 2316 |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $NCORE -o $HIC_BAM

echo -e "Generating BAM stats file. Command:\n samtools stats $HIC_BAM > $HIC_BAM.stats"
samtools stats $HIC_BAM > $HIC_BAM.stats

echo -e "\nWriting files ${SCAF_ASM_FN%.*}.agp and ${SCAF_ASM_FN%.*}.assembly with juicebox_scripts:"
grep -H version $SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/setup.py
echo -e "\ncommands:\n $SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/makeAgpFromFasta.py\
 $SCAF_ASM_FN ${SCAF_ASM_FN%.*}.agp"
echo -e " $SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/agp2assembly.py\
 ${SCAF_ASM_FN%.*}.agp ${SCAF_ASM_FN%.*}.assembly"
$SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/makeAgpFromFasta.py\
 $SCAF_ASM_FN ${SCAF_ASM_FN%.*}.agp
$SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/agp2assembly.py\
 ${SCAF_ASM_FN%.*}.agp ${SCAF_ASM_FN%.*}.assembly

echo -e "\nPreparing mapped pair data for juicer_tools using yahs juicer pre:"
echo $SOFTWARE_DIR/yahs-juicer
$SOFTWARE_DIR/yahs-juicer/juicer --help | grep Version
echo -e "\ncommand:"
echo " $SOFTWARE_DIR/yahs-juicer/juicer pre\
 -a -o ${SCAF_ASM_FN%.*}.JBAT $HIC_BAM ${SCAF_ASM_FN%.*}.agp\
 ${SCAF_ASM_FN}.fai 2> ${SCAF_ASM_FN%.*}.JBAT/tmp_juicer_pre_JBAT.log"

mkdir ${SCAF_ASM_FN%.*}.JBAT
$SOFTWARE_DIR/yahs-juicer/juicer pre \
 -a -o ${SCAF_ASM_FN%.*}.JBAT $HIC_BAM ${SCAF_ASM_FN%.*}.agp \
 ${SCAF_ASM_FN}.fai 2> ${SCAF_ASM_FN%.*}.JBAT/tmp_juicer_pre_JBAT.log

echo -e "\nMaking contacts file ${SCAF_ASM_FN%.*}.JBAT.hic using juicer_tools singularity image"
echo -e "command:\n singularity exec $SOFTWARE_DIR/NextflowSingularityLibrary/dmolik-juicer-tools.img\
 java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar\
 pre ${SCAF_ASM_FN%.*}.JBAT.txt ${SCAF_ASM_FN%.*}.JBAT.hic.part\
 <(cat ${SCAF_ASM_FN%.*}.JBAT/tmp_juicer_pre_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}') &&\
 mv ${SCAF_ASM_FN%.*}.JBAT.hic.part ${SCAF_ASM_FN%.*}.JBAT.hic"

singularity exec $SOFTWARE_DIR/NextflowSingularityLibrary/dmolik-juicer-tools.img java \
  -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre ${SCAF_ASM_FN%.*}.JBAT.txt ${SCAF_ASM_FN%.*}.JBAT.hic.part \
  <(cat ${SCAF_ASM_FN%.*}.JBAT/tmp_juicer_pre_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}') && \
  mv ${SCAF_ASM_FN%.*}.JBAT.hic.part ${SCAF_ASM_FN%.*}.JBAT.hic

echo -e "\nFinished\n"
date
