#!/bin/bash

# Reciprocal_RagTag.sh runs one leg of reciprocal ragtag
# Takes reference assembly (to use as curration guide), 
#  query assembly (to currate), and 
#  reads for query assembly
# Runs 2 rounds of RagTag correct/scaffold with GFF shielding "better in qry" regions
#  Round 1 between qry and ref
#  Round 2 beteen round1 output and original qry assembly
# GFFs produced by RagTagShield-WriteGFF.R, 
#  which uses nucmer COORDS output to ID alignment regions where qry assembly has greater contiguity
#  writes GFF formated file with these regions to give RagTag correct, preventing breaks in these regions

# USAGE: sbatch Reciprocal_RagTag.slurm REF.fa QRY.fa QRY_READS.fq N_pass (number for which pass in pipeline 1-3)
REF_ASM=$1
QRY_ASM=$2
QRY_READS=$3
# added in .sh version for running in SnakeMake pipeline
PASS_TAG=$4

# check reads file (3rd argument) exists 
[[ -f $QRY_READS ]] || \
  { echo "USAGE: sbatch Reciprocal_RagTag.slurm REF.fa QRY.fa QRY_READS.fq"; exit; }

date
echo "Starting Reciprocal_RagTag"
echo "reference: $REF_ASM"
echo "subject/query: $QRY_ASM"
echo "reads: $QRY_READS"

module load mummer/4.0.0rc1
module load r/4.4.1
module load ragtag/2.1.0
module load samtools/1.17

# location of FlyComparativeGenomics/ for RagTagShield-writeGFF.R, gfastats.slurm, 
GIT_REPOS=~

# Print message for start of jobs, smk pipeline appends to a continuous log so want to distinguish runs
echo -e "\n\n======================================" | tee /dev/stderr
date | tee /dev/stderr
echo -e "Starting new Reciprocal_RagTag job.\nREF_ASM:\
 ${REF_ASM}\nQRY_ASM: ${QRY_ASM}\nQRY_READS: ${QRY_READS}\nPASS_TAG: ${PASS_TAG}\n\n" | tee /dev/stderr

# ROUND 1 - rescaffold QRY_ASM based on REF_ASM

# Check for prior successful round1 run using ragtag_round1-pct_chrs.txt, if round1 successful will return integer
R1_CHR_N=$(grep --no-messages ChromNumber ragtag_round1-pct_chrs.txt | cut -f2)

# skip round 1 if integer round 1 chromosome count and round 1 output exists
if [[ $R1_CHR_N =~ $re_int && -f ragtag_output-round1/ragtag.scaffold.fasta ]]; then

  echo -n "\nFound output from successful round 1."
  echo "Skipping to round 2 using ${QRY_ASM} and ragtag_output-round1/ragtag.scaffold.fasta"

else

  # Alignment between ref and qry for RagTagShield
  REFASM_FILENAME=$(basename $REF_ASM)
  QRYASM_FILENAME=$(basename $QRY_ASM)
  ALIGNMENT1_NAME=${REFASM_FILENAME%.*}-vs-${QRYASM_FILENAME%.*}-c5000

  echo -e "\nRunning genome alignment to identify regions to preserve in $QRY_ASM"
  echo -e "Commands:\n nucmer -p ${ALIGNMENT1_NAME} --threads=32 -c 5000 ${REF_ASM} ${QRY_ASM}"
  echo -e " show-coords ${ALIGNMENT1_NAME}.delta > ${ALIGNMENT1_NAME}.coords\n"

  nucmer -p ${ALIGNMENT1_NAME} --threads=32 -c 5000 ${REF_ASM} ${QRY_ASM}
  show-coords ${ALIGNMENT1_NAME}.delta > ${ALIGNMENT1_NAME}.coords

  # RagTagShield for first round of ragtag
  echo -e "\nWriting GFF for regions to preserve in $QRY_ASM"
  echo -e "Command:\n Rscript $GIT_REPOS/FlyComparativeGenomics/RagTagShield-writeGFF.R ${ALIGNMENT1_NAME}.coords"

  Rscript $GIT_REPOS/FlyComparativeGenomics/RagTagShield-writeGFF.R ${ALIGNMENT1_NAME}.coords

  # RagTag round 1
  # run correct to break contigs if alignment and reads support
  # --debug to provide more info on breaks
  # --intra only breaks when query maps to single ref sequence, 
  #    prevents splitting out unplaced reference scaffolds
  # -f minimum alignment length to consider
  # -t threads
  # --gff defines regions to preserve in qry assembly
  # -T corr for high quality long reads
  # -o output directory
  echo -e "\nStarting RagTag Round 1"
  echo -e "Commands:\n ragtag.py correct --debug --intra -f 5000 -t 32 --gff ${ALIGNMENT1_NAME}.coords-BetterInQry.gff -R $QRY_READS -T corr -o ./ragtag_output-round1 $REF_ASM $QRY_ASM"
  echo -e " ragtag.py scaffold -r -t 32 -o ./ragtag_output-round1 $REF_ASM ragtag_output-round1/ragtag.correct.fasta\n"

  ragtag.py correct --debug --intra -f 5000 -t 32 \
    --gff ${ALIGNMENT1_NAME}.coords-BetterInQry.gff \
    -R $QRY_READS -T corr \
    -o ./ragtag_output-round1 \
    $REF_ASM $QRY_ASM

  # run scaffold to reorder corrected QRY scaffolds
  # -r to estimate gap size
  ragtag.py scaffold -r -t 32 -o ./ragtag_output-round1 $REF_ASM ragtag_output-round1/ragtag.correct.fasta

  # A known issue in RagTag can introduce duplicate sequence header names: https://github.com/malonge/RagTag/issues/94#issuecomment-1006880604
  # This case is slightly different.  Not using -u, but the "_RagTag" label on modified scaffolds can still cause duplicates in round 2.
  # Similar issue can also crop up between passes in the larger pipeline (3 total passes, 2 rounds each)
  # Easiest fix is editing this tag in the round1 output. Also has the benefit of tracking which pass and round introduced changes
  echo -e "\nEditing sequence headers in pass $PASS_TAG round 1 output to avoid duplicates in round 2 or future passes"
  echo -e "Command:\n sed -i 's/RagTag/RTp${PASS_TAG}r1/' ragtag_output-round1/ragtag.scaffold.fasta"

  sed -i "s/RagTag/RTp${PASS_TAG}r1/" ragtag_output-round1/ragtag.scaffold.fasta

  # stats for new assembly
  echo -e "\nRunning gfastats and get_pct_chrs.sh on ragtag_output-round1/ragtag.scaffold.fasta"
  sbatch $GIT_REPOS/FlyComparativeGenomics/gfastats.slurm ragtag_output-round1/ragtag.scaffold.fasta
  samtools faidx ragtag_output-round1/ragtag.scaffold.fasta
  $GIT_REPOS/FlyComparativeGenomics/get_pct_chrs.sh ragtag_output-round1/ragtag.scaffold.fasta.fai \
    > ragtag_round1-pct_chrs.txt

fi

# ROUND 2 - rescaffold round 1 output based on QRY_ASM (original input)
# RagTagShield will preserve regions of improvement

# Alignment between ragtag scaffolds and initial assembly for RagTagShield
ALIGNMENT2_NAME=${QRYASM_FILENAME%.*}-vs-ragtag_round1-c5000

echo -e "\nRunning genome alignment to find improved regions to preserve in ragtag_output-round1/ragtag.scaffold.fasta"
echo -e "Commands:\n nucmer -p ${ALIGNMENT2_NAME} --threads=32 -c 5000 ${QRY_ASM} ragtag_output-round1/ragtag.scaffold.fasta"
echo -e " show-coords ${ALIGNMENT2_NAME}.delta > ${ALIGNMENT2_NAME}.coords\n"

nucmer -p ${ALIGNMENT2_NAME} --threads=32 -c 5000 ${QRY_ASM} ragtag_output-round1/ragtag.scaffold.fasta
show-coords ${ALIGNMENT2_NAME}.delta > ${ALIGNMENT2_NAME}.coords

# RagTagShield for 2nd round
echo -e "\nWriting GFF for regions to preserve in ragtag_output-round1/ragtag.scaffold.fasta"
echo -e "Command:\n Rscript $GIT_REPOS/FlyComparativeGenomics/RagTagShield-writeGFF.R ${ALIGNMENT2_NAME}.coords"
Rscript $GIT_REPOS/FlyComparativeGenomics/RagTagShield-writeGFF.R ${ALIGNMENT2_NAME}.coords

# RagTag round 2 - QRY_ASM is now reference for ragtag_round1
echo -e "\nStarting RagTag Round 2"
echo -e "Commands:\n ragtag.py correct --debug --intra -f 5000 -t 32 --gff ${ALIGNMENT2_NAME}.coords-BetterInQry.gff -R $QRY_READS -T corr -o ./ragtag_output-round2 $QRY_ASM ragtag_output-round1/ragtag.scaffold.fasta"
echo -e " ragtag.py scaffold -r -t 32 -o ./ragtag_output-round2 $QRY_ASM ragtag_output-round2/ragtag.correct.fasta\n"

ragtag.py correct --debug --intra -f 5000 -t 32 \
  --gff ${ALIGNMENT2_NAME}.coords-BetterInQry.gff \
  -R $QRY_READS -T corr \
  -o ./ragtag_output-round2 \
  $QRY_ASM ragtag_output-round1/ragtag.scaffold.fasta

ragtag.py scaffold -r -t 32 -o ./ragtag_output-round2 $QRY_ASM ragtag_output-round2/ragtag.correct.fasta

echo -e "\nEditing sequence headers in output of pass $PASS_TAG round 2 to avoid potential duplicates in subsequent passes"
echo -e "Command:\n sed -i 's/RagTag/RTp${PASS_TAG}r2/' ragtag_output-round2/ragtag.scaffold.fasta"
sed -i "s/RagTag/RTp${PASS_TAG}r2/" ragtag_output-round2/ragtag.scaffold.fasta

# stats for new assembly
sbatch $GIT_REPOS/FlyComparativeGenomics/gfastats.slurm ragtag_output-round2/ragtag.scaffold.fasta
samtools faidx ragtag_output-round2/ragtag.scaffold.fasta
$GIT_REPOS/FlyComparativeGenomics/get_pct_chrs.sh ragtag_output-round2/ragtag.scaffold.fasta.fai \
  > ragtag_round2-pct_chrs.txt

