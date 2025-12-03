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

# script for regenerating clean HiC contact map on scaffolded assembly
# used to distinguish JBAT artifacts from misassemblies, fine tune curation
# Steps:
# split scaffolds into contigs (retaining order)
#   using split_fa from purge_dups
# index split contigs assembly
#   using samtools faidx
# map HiC reads to assembly contigs, filter based on HiC methods (same as OTB)
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
module load purge_dups/1.2.6
module load apptainer

SCAF_ASM=$1
HIC_R1=$2
HIC_R2=$3

NCORE=16

SCAF_ASM_FN=$(basename $SCAF_ASM)
CTG_ASM=${SCAF_ASM_FN%.*}.split.fa
HIC_BAM=${CTG_ASM%.*}.HiC.bam

date

echo -e "\nSplitting scaffold assembly $SCAF_ASM into contig assembly $CTG_ASM using split_fa:"
which split_fa
echo -e "\ncommand:\n split_fa $SCAF_ASM > $CTG_ASM"
split_fa $SCAF_ASM > $CTG_ASM

module unload purge_dups # wasn't playing nice with samtools

echo -e "\nIndexing contig assembly with samtools:"
which samtools
echo -e "\ncommand:\n samtools faidx $CTG_ASM"
samtools faidx $CTG_ASM

echo -e "\nMapping HiC reads $HIC_R1 and $HIC_R2 onto $CTG_ASM using bwa-mem2:"
which bwa-mem2
echo -e "and filtering using samtools and samblaster:"
which samtools
$SOFTWARE_DIR/samblaster/samblaster --version
echo -e "\nOutput file: $HIC_BAM"
echo -e "\ncommands:\n bwa-mem2 index -p ${CTG_ASM%.*} $CTG_ASM"
echo -e " bwa-mem2 mem -5SP -t $NCORE ${CTG_ASM%.*} $HIC_R1 $HIC_R2 |\
 samtools view -h -F 2316 |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $NCORE -o $HIC_BAM"
bwa-mem2 index -p ${CTG_ASM%.*} $CTG_ASM
bwa-mem2 mem -5SP -t $NCORE ${CTG_ASM%.*} $HIC_R1 $HIC_R2 |\
 samtools view -h -F 2316 |\
 $SOFTWARE_DIR/samblaster/samblaster |\
 samtools sort -n -@ $NCORE -o $HIC_BAM

echo -e "Generating BAM stats file. Command:\n samtools stats $HIC_BAM > $HIC_BAM.stats"
samtools stats $HIC_BAM > $HIC_BAM.stats

echo -e "\nWriting files ${CTG_ASM%.*}.agp and ${CTG_ASM%.*}.assembly with juicebox_scripts:"
grep -H version $SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/setup.py
echo -e "\ncommands:\n $SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/makeAgpFromFasta.py\
 $CTG_ASM ${CTG_ASM%.*}.agp"
echo -e " $SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/agp2assembly.py\
 ${CTG_ASM%.*}.agp ${CTG_ASM%.*}.assembly"
$SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/makeAgpFromFasta.py\
 $CTG_ASM ${CTG_ASM%.*}.agp
$SOFTWARE_DIR/juicebox_scripts/juicebox_scripts/agp2assembly.py\
 ${CTG_ASM%.*}.agp ${CTG_ASM%.*}.assembly

echo -e "\nPreparing mapped pair data for juicer_tools using yahs juicer pre:"
echo $SOFTWARE_DIR/yahs-juicer
$SOFTWARE_DIR/yahs-juicer/juicer --help | grep Version
echo -e "\ncommand:"
echo " $SOFTWARE_DIR/yahs-juicer/juicer pre\
 -a -o ${CTG_ASM%.*}.JBAT $HIC_BAM ${CTG_ASM%.*}.agp\
 ${CTG_ASM}.fai 2> ${CTG_ASM%.*}.JBAT/tmp_juicer_pre_JBAT.log"

mkdir ${CTG_ASM%.*}.JBAT
$SOFTWARE_DIR/yahs-juicer/juicer pre \
 -a -o ${CTG_ASM%.*}.JBAT $HIC_BAM ${CTG_ASM%.*}.agp \
 ${CTG_ASM}.fai 2> ${CTG_ASM%.*}.JBAT/tmp_juicer_pre_JBAT.log

echo -e "\nMaking contacts file ${CTG_ASM%.*}.JBAT.hic using juicer_tools singularity image"
echo -e "command:\n singularity exec $SOFTWARE_DIR/NextflowSingularityLibrary/dmolik-juicer-tools.img\
 java -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar\
 pre ${CTG_ASM%.*}.JBAT.txt ${CTG_ASM%.*}.JBAT.hic.part\
 <(cat ${CTG_ASM%.*}.JBAT/tmp_juicer_pre_JBAT.log | grep PRE_C_SIZE | awk '{print \$2\" \"\$3}') &&\
 mv ${CTG_ASM%.*}.JBAT.hic.part ${CTG_ASM%.*}.JBAT.hic"

singularity exec $SOFTWARE_DIR/NextflowSingularityLibrary/dmolik-juicer-tools.img java \
  -Xms16384m -Xmx32768m -jar /home/genomics/juicer_tools_1.22.01.jar pre ${CTG_ASM%.*}.JBAT.txt ${CTG_ASM%.*}.JBAT.hic.part \
  <(cat ${CTG_ASM%.*}.JBAT/tmp_juicer_pre_JBAT.log | grep PRE_C_SIZE | awk '{print $2" "$3}') && \
  mv ${CTG_ASM%.*}.JBAT.hic.part ${CTG_ASM%.*}.JBAT.hic

echo -e "\nFinished\n"
date
