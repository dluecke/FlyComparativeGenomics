#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J rscript_plotasm
#SBATCH -o "pac.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "pac.%j.%N.stderr"
#SBATCH -c 2
#SBATCH -p ceres
#SBATCH -t 1:00:00
#SBATCH -N 1

# PlotsForAssemblyCuration.slurm takes Mummer output and makes standard set of plots
# to visualize assembly agreement, ID unplaced scaffolds, and check for duplication

# Intended to run via command line submission on Ceres

# USAGE: sbatch PlotsForAssemblyCuration.slurm INFILES_ALN.csv INFILES_SELFALN.csv [OUTFILENAME default Alignments-INFILES]

# INFILES_ALN.csv has alignment id, refASM, qryASM, COORDS, FAI, ref,qry axis labels, no header:
#     pctgs,asmM,asmF,asmM-vs-asmF.coords,asmM.fa.fai,asmF.fa.fai,male primary,female primary
#     hapF1,asmF,hapF1,asmF-vs-Fh1.coords,asmF.fa.fai,Fh1.fa.fai,female primary,female hap1
# Primary vs Primary needs to be first in this list

# INFILES_SELFALN.csv is for self alignments (unplaced scaffolds onto chrs) 
#  can run from single assembly via self_align-ChrsVsUnplaced.sh
# alignment id (all must match refASM/qryASM and haps must match alignment ids in INFILES_ALN.csv), 
#  COORDS, and FAI for self alignments, no headers:
#     Fasm,FChrs-vs-FUnplaced.coords,FChrs.fa.fai,FUnplaced.fa.fai
#     hapF1,Fh1Chrs-vs-Fh1Unplaced.coords,Fh1Chrs.fa.fai,Fh1Unplaced.fa.fai

module load r/4.4.1

IN_ALN=$1
IN_SELF=$2

OUTNAME="Alignments-"$(echo $IN_ALN | cut -d'_' -f1)

# repo location
FCG=~/FlyComparativeGenomics

Rscript $FCG/PlotsForAssemblyCuration.R $IN_ALN $IN_SELF $OUTNAME
