#!/bin/bash

# Rephase_XY_contigs.sh collects scaffolds and contigs for X/Y phases, starts rescaffolding
# Takes:
#  Full male diploid assembly with non-autosomes split to contigs,
#  X reference scaffold (eg from female primary)
#  lists for autosome scaffolds to pair with X and Y phases (better hap with X)
#  ctg lists for X and Y from rephase_RTout-V2.R.
#  HiC fastq files R1/R2
# Compiles new phases with autosome scaffolds and contig set
# Starts RagTag scaffolding for X phase with X reference scaffold
# Starts YaHS for Y contigs with yahs_subset.slurm

# USAGE Rephase_XY_contigs.sh INPUT.json ASM_NAME_PREFIX
# INPUT.json is a json file with the following keys:
#   MALE_DIPLOID: Absolute path to the full male diploid assembly
#   X_REF: Absolute path to the X reference scaffold 
#   PHASE_X_AUTOSOMES: Absolute path to the list of autosome scaffolds for X phase
#   PHASE_Y_AUTOSOMES: Absolute path to the list of autosome scaffolds for Y phase
#   X_CTGS: Absolute path to the list of X contigs
#   Y_CTGS: Absolute path to the list of Y contigs
#   HIC_FQ: [Absolute path to HiCR1.fq, Absolute path to HiCR2.fq]

# REQUIREMENTS: samtools
GIT_REPO=~/FlyComparativeGenomics
module load samtools/1.17

if [[ $# -ne 2 ]]; then
  echo "USAGE: Rephase_XY_contigs.sh INPUT.json ASM_NAME_PREFIX"
  exit 1
fi

INPUT_JSON=$1
ASM_ID=$2

# phased scaffold+contig assembly filenames (pre-scaffolding)
X_PHASE_ASM_V0="${ASM_ID}-X_v0.fa"
Y_PHASE_ASM_V0="${ASM_ID}-Y_v0.fa"

INPUT_JSON_FN=$(basename $INPUT_JSON)
LOGFILE="Rephase_XY_contigs-${INPUT_JSON_FN%.json}.log"

if [[ -f $LOGFILE ]]; then
  cat $LOGFILE >> "${LOGFILE}.bak"
  echo -e "\n\n\n" >> "${LOGFILE}.bak"
fi

date > $LOGFILE
echo "Starting Rephase_XY_contigs.sh" >> $LOGFILE
echo "INPUT_JSON: $INPUT_JSON" >> $LOGFILE
echo "ASM_ID: $ASM_ID" >> $LOGFILE
echo

samtools --version >> $LOGFILE

DIPLOID_ASM=$(jq -r '.MALE_DIPLOID' $INPUT_JSON)
X_REF_FA=$(jq -r '.X_REF' $INPUT_JSON) 
X_AUTOSOMES_LIST=$(jq -r '.PHASE_X_AUTOSOMES' $INPUT_JSON)
Y_AUTOSOMES_LIST=$(jq -r '.PHASE_Y_AUTOSOMES' $INPUT_JSON)
X_CTGS_LIST=$(jq -r '.X_CTGS' $INPUT_JSON)
Y_CTGS_LIST=$(jq -r '.Y_CTGS' $INPUT_JSON)
HIC_FQ1=$(jq -r '.HIC_FQ[0]' $INPUT_JSON)
HIC_FQ2=$(jq -r '.HIC_FQ[1]' $INPUT_JSON)

dos2unix $X_CTGS_LIST
dos2unix $Y_CTGS_LIST

echo -e "\nCompiling X phase assembly with autosomes and X contigs" >> $LOGFILE
echo "CMD: samtools faidx $DIPLOID_ASM -r $X_AUTOSOMES_LIST > $X_PHASE_ASM_V0" >> $LOGFILE
samtools faidx $DIPLOID_ASM -r $X_AUTOSOMES_LIST > $X_PHASE_ASM_V0
echo "CMD: samtools faidx $DIPLOID_ASM -r $X_CTGS_LIST >> $X_PHASE_ASM_V0" >> $LOGFILE
samtools faidx $DIPLOID_ASM -r $X_CTGS_LIST >> $X_PHASE_ASM_V0

echo -e "\nCompiling Y phase assembly with autosomes and Y contigs" >> $LOGFILE
echo "CMD: samtools faidx $DIPLOID_ASM -r $Y_AUTOSOMES_LIST > $Y_PHASE_ASM_V0" >> $LOGFILE
samtools faidx $DIPLOID_ASM -r $Y_AUTOSOMES_LIST > $Y_PHASE_ASM_V0
echo "CMD: samtools faidx $DIPLOID_ASM -r $Y_CTGS_LIST >> $Y_PHASE_ASM_V0" >> $LOGFILE
samtools faidx $DIPLOID_ASM -r $Y_CTGS_LIST >> $Y_PHASE_ASM_V0

echo -e "\nStatistics for V0 phase assemblies" >> $LOGFILE
echo "CMD: $GIT_REPO/gfastats.slurm $X_PHASE_ASM_V0" >> $LOGFILE
echo "CMD: $GIT_REPO/gfastats.slurm $Y_PHASE_ASM_V0" >> $LOGFILE
sbatch $GIT_REPO/gfastats.slurm $X_PHASE_ASM_V0
sbatch $GIT_REPO/gfastats.slurm $Y_PHASE_ASM_V0

mkdir -p ragtag_Xctgs
cd ragtag_Xctgs
echo -e "\nStarting RagTag scaffolding for X phase assembly in $PWD" >> ../$LOGFILE
echo "Compiling X contig sequences: ${ASM_ID}-X_ctgs.fa" >> ../$LOGFILE
echo "CMD: samtools faidx $DIPLOID_ASM -r $X_CTGS_LIST > ${ASM_ID}-X_ctgs.fa" >> ../$LOGFILE
samtools faidx $DIPLOID_ASM -r $X_CTGS_LIST > ${ASM_ID}-X_ctgs.fa
echo -e "\nSubmitting RagTag scaffolding for X phase assembly" >> ../$LOGFILE
echo "CMD: sbatch $GIT_REPO/ragtag_scaffold.slurm $X_REF_FA ${ASM_ID}-X_ctgs.fa" >> ../$LOGFILE
sbatch $GIT_REPO/ragtag_scaffold.slurm $X_REF_FA ${ASM_ID}-X_ctgs.fa
cd ..

mkdir yahs_Yctgs
cd yahs_Yctgs
echo -e "\nStarting YaHS scaffolding for Y phase assembly" >> ../$LOGFILE
echo "CMD: sbatch $GIT_REPO/yahs_subset.slurm ../$Y_PHASE_ASM_V0 $HIC_FQ1 $HIC_FQ2 $Y_CTGS_LIST $ASM_ID-Y" >> ../$LOGFILE
sbatch $GIT_REPO/yahs_subset.slurm ../$Y_PHASE_ASM_V0 $HIC_FQ1 $HIC_FQ2 $Y_CTGS_LIST $ASM_ID"-Y"
cd ..

echo
date >> $LOGFILE
