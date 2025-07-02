#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="recip_ragtag"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 5-00:00:00
#SBATCH --mem-per-cpu=16G   # replacing mem queue
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "rrt.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "rrt.%j.%N.stderr"     #optional, prints our standard error

# SLURM submission to run one leg of reciprocal ragtag
# Takes reference assembly (to use as curration guide), 
#  query assembly (to currate), and 
#  reads for query assembly
# Runs 2 rounds of RagTag correct/scaffold with GFF shielding "better in qry" regions
#  Round 1 between qry and ref
#  Round 2 beteen round1 output and original qry assembly
# GFFs produced by RagTagShield-WriteGFF.R, 
#  which uses nucmer COORDS output to ID alignment regions where qry assembly has greater contiguity
#  writes GFF formated file with these regions to give RagTag correct, preventing breaks in these regions

# USAGE: sbatch Reciprocal_RagTag.slurm REF.fa QRY.fa QRY_READS.fq
REF_ASM=$1
QRY_ASM=$2
QRY_READS=$3

# check reads file (3rd argument) exists 
[[ -f $QRY_READS ]] || (echo "USAGE: sbatch Reciprocal_RagTag.slurm REF.fa QRY.fa QRY_READS.fq"; exit;)

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

# ROUND 1 - rescaffold QRY_ASM based on REF_ASM

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

# stats for new assembly
sbatch $GIT_REPOS/FlyComparativeGenomics/gfastats.slurm ragtag_output-round1/ragtag.scaffold.fasta
samtools faidx ragtag_output-round1/ragtag.scaffold.fasta
$GIT_REPOS/FlyComparativeGenomics/get_pct_chrs.sh ragtag_output-round1/ragtag.scaffold.fasta.fai \
  > ragtag_round1-pct_chrs.txt


# ROUND 2 - rescaffold round 1 output based on QRY_ASM (original input)
# RagTagShield will preserve regions of improvement

# Alignment between ragtag scaffolds and initial assembly for RagTagShield
ALIGNMENT2_NAME=${QRYASM_FILENAME%.*}-vs-ragtag_round1-c5000

nucmer -p ${ALIGNMENT2_NAME} --threads=32 -c 5000 ${REF_ASM} ${QRY_ASM}

show-coords ${ALIGNMENT2_NAME}.delta > ${ALIGNMENT2_NAME}.coords

# RagTagShield for 2nd round

Rscript $GIT_REPOS/FlyComparativeGenomics/RagTagShield-writeGFF.R ${ALIGNMENT2_NAME}.coords

# RagTag round 2 - QRY_ASM is now reference for ragtag_round1
ragtag.py correct --debug --intra -f 5000 -t 32 \
  --gff ${ALIGNMENT2_NAME}.coords-BetterInQry.gff  \
  -R $QRY_READS -T corr \
  -o ./ragtag_output-round2 \
  $QRY_ASM ragtag_output-round1/ragtag.scaffold.fasta

ragtag.py scaffold -r -t 32 -o ./ragtag_output-round2 $QRY_ASM ragtag_output-round2/ragtag.correct.fasta

# stats for new assembly
sbatch $GIT_REPOS/FlyComparativeGenomics/gfastats.slurm ragtag_output-round2/ragtag.scaffold.fasta
samtools faidx ragtag_output-round2/ragtag.scaffold.fasta
$GIT_REPOS/FlyComparativeGenomics/get_pct_chrs.sh ragtag_output-round2/ragtag.scaffold.fasta.fai \
  > ragtag_round2-pct_chrs.txt

