#!/bin/bash

# Rephase_XY_scaffolds.sh finished rescaffolding for X/Y phases, starts split HiC mapping
# Takes same JSON inputs as Rephase_XY_contigs.sh, of which uses:
#  Full male diploid assembly with non-autosomes split to contigs
#  X reference scaffold (eg from female primary)
#  lists for autosome scaffolds to pair with X and Y phases (better hap with X)
#  HiC fastq files R1/R2
# Also uses outputs from Rephase_XY_contigs.sh:
#  yahs_Yctgs/[ASM_TAG]-Y_scaffolds_final.fa
#  ragtag_Xctgs/ragtag_output/ragtag.scaffold.fasta
# Runs ragtag_scaffold.sh to scaffold Y yahs output onto reference X scaffold
# Compiles X and Y phase scaffolds from autosomes and ragtag output
# Starts final HiC mapping with redo_hic_contacts.slurm

# Better to run from compute node for ragtag_scaffold.sh step
# USAGE Rephase_XY_scaffolds.sh INPUT.json ASM_NAME_PREFIX
# Both inputs should be same same as for previous Rephase_XY_contigs.sh run
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

if [[ $# -ne 2 || ! -d ragtag_Xctgs || ! -d yahs_Yctgs ]]; then
  echo "USAGE: Rephase_XY_scaffolds.sh INPUT.json ASM_NAME_PREFIX"
  echo "Need to have run Rephase_XY_contigs.sh first in same directory"
  echo "Both inputs should be same same as previous Rephase_XY_contigs.sh"
  exit 1
fi

INPUT_JSON=$1
ASM_ID=$2

# phased scaffold+contig assembly filenames (pre-scaffolding)
X_PHASE_ASM_V1="${ASM_ID}-X_v1.fa"
Y_PHASE_ASM_V1="${ASM_ID}-Y_v1.fa"

INPUT_JSON_FN=$(basename $INPUT_JSON)
LOGFILE="Rephase_XY_scaffolds-${INPUT_JSON_FN%.json}.log"

if [[ -f $LOGFILE ]]; then
  cat $LOGFILE >> "${LOGFILE}.bak"
  echo -e "\n\n\n" >> "${LOGFILE}.bak"
fi

date > $LOGFILE
echo "Starting Rephase_XY_scaffolds.sh" >> $LOGFILE
echo "INPUT_JSON: $INPUT_JSON" >> $LOGFILE
echo "ASM_ID: $ASM_ID" >> $LOGFILE
echo

samtools --version >> $LOGFILE

DIPLOID_ASM=$(jq -r '.MALE_DIPLOID' $INPUT_JSON)
X_REF_FA=$(jq -r '.X_REF' $INPUT_JSON) 
X_AUTOSOMES_LIST=$(jq -r '.PHASE_X_AUTOSOMES' $INPUT_JSON)
Y_AUTOSOMES_LIST=$(jq -r '.PHASE_Y_AUTOSOMES' $INPUT_JSON)
HIC_FQ1=$(jq -r '.HIC_FQ[0]' $INPUT_JSON)
HIC_FQ2=$(jq -r '.HIC_FQ[1]' $INPUT_JSON)

X_PHASE_RAGTAG_OUT="$PWD/ragtag_Xctgs/ragtag_output/ragtag.scaffold.fasta"
Y_PHASE_YAHS_OUT="$PWD/yahs_Yctgs/${ASM_ID}-Y_scaffolds_final.fa"

if [[ ! -f $X_PHASE_RAGTAG_OUT || ! -f $Y_PHASE_YAHS_OUT ]]; then
  echo "Missing ragtag or yahs output files, check previous steps" >> $LOGFILE
  exit 1
fi

echo "X phase ragtag output: $X_PHASE_RAGTAG_OUT" >> $LOGFILE
echo "Y phase yahs output: $Y_PHASE_YAHS_OUT" >> $LOGFILE

mkdir -p ragtag_Yctgs
cd ragtag_Yctgs

echo -e "\nScaffolding Y yahs output onto X reference scaffold with ragtag_scaffold.sh in $PWD" >> $LOGFILE
echo "CMD: $GIT_REPO/ragtag_scaffold.sh $X_REF_FA $Y_PHASE_YAHS_OUT 2" >> $LOGFILE
$GIT_REPO/ragtag_scaffold.sh $X_REF_FA $Y_PHASE_YAHS_OUT 2 >> $LOGFILE 2>&1
date >> $LOGFILE
Y_PHASE_RAGTAG_OUT="$PWD/ragtag_output/ragtag.scaffold.fasta"
if [[ ! -f $Y_PHASE_RAGTAG_OUT ]]; then
  echo "Missing ragtag output file, check ragtag_scaffold.sh step" >> $LOGFILE
  exit 1
fi
echo "RagTag Y scaffolding finished" >> $LOGFILE

echo -e "\nCompiling X phase assembly with autosomes and ragtag output" >> $LOGFILE
echo "CMD: samtools faidx $DIPLOID_ASM -r $X_AUTOSOMES_LIST > $X_PHASE_ASM_V1" >> $LOGFILE
samtools faidx $DIPLOID_ASM -r $X_AUTOSOMES_LIST > $X_PHASE_ASM_V1
echo "CMD: cat $X_PHASE_RAGTAG_OUT >> $X_PHASE_ASM_V1" >> $LOGFILE
cat $X_PHASE_RAGTAG_OUT >> $X_PHASE_ASM_V1
echo "CMD: sbatch $GIT_REPO/gfastats.slurm $X_PHASE_ASM_V1" >> $LOGFILE
sbatch $GIT_REPO/gfastats.slurm $X_PHASE_ASM_V1

echo -e "\nCompiling Y phase assembly with autosomes and ragtag output" >> $LOGFILE
echo "CMD: samtools faidx $DIPLOID_ASM -r $Y_AUTOSOMES_LIST > $Y_PHASE_ASM_V1" >> $LOGFILE
samtools faidx $DIPLOID_ASM -r $Y_AUTOSOMES_LIST > $Y_PHASE_ASM_V1
echo "CMD: cat $Y_PHASE_RAGTAG_OUT >> $Y_PHASE_ASM_V1" >> $LOGFILE
cat $Y_PHASE_RAGTAG_OUT >> $Y_PHASE_ASM_V1
echo "CMD: sbatch $GIT_REPO/gfastats.slurm $Y_PHASE_ASM_V1" >> $LOGFILE
sbatch $GIT_REPO/gfastats.slurm $Y_PHASE_ASM_V1

echo -e "\nRescaffolded X and Y phase assemblies finished" >> $LOGFILE
echo -e "\nStarting final HiC mapping with redo_hic_contacts.slurm" >> $LOGFILE

echo -e "CMD: sbatch $GIT_REPO/redo_hic_contacts.slurm $X_PHASE_ASM_V1 $HIC_FQ1 $HIC_FQ2" >> $LOGFILE
sbatch $GIT_REPO/redo_hic_contacts.slurm $X_PHASE_ASM_V1 $HIC_FQ1 $HIC_FQ2 >> $LOGFILE 2>&1
echo -e "CMD: sbatch $GIT_REPO/redo_hic_contacts.slurm $Y_PHASE_ASM_V1 $HIC_FQ1 $HIC_FQ2" >> $LOGFILE
sbatch $GIT_REPO/redo_hic_contacts.slurm $Y_PHASE_ASM_V1 $HIC_FQ1 $HIC_FQ2 >> $LOGFILE 2>&1

echo -e "\nRephase_XY_scaffolds.sh finished" >> $LOGFILE
date >> $LOGFILE

