#!/bin/bash

# PrepareCurationPlots.sh takes paths to female and male primary and hap assemblies
#  and prepares for PlotsForAssemblyCuration.R
# Preps and starts all necessary alignments, including self-align of small scaffolds onto chr-scale
# Writes the CSV submission files used by PlotsForAssemblyCuration.R

# Input is CSV with columns assembly ID, absolute paths, and labels (no commas)
# NO HEADER, ROWS NEED TO BE IN THIS ORDER (or flip male/female):
# Fpri,/path/to/Fpri.fa,Spp female primary scaffolds
# hapF1,/path/to/hapF1.fa,Spp female hap1 scaffolds
# hapF2,/path/to/hapF2.fa,Spp female hap2 scaffolds
# Mpri,/path/to/Mpri.fa,Spp male primary scaffolds
# hapM1,/path/to/hapM1.fa,Spp male hap1 scaffolds
# hapM2,/path/to/hapM2.fa,Spp male hap2 scaffolds

# also need N_CHR, number of chromosomes to use for self-alignment

IN_CSV=$1
N_CHR=$2
[[ -f $IN_CSV && $N_CHR -gt 0 ]] || \
 { echo "Usage: PrepareCurationPlots.sh asmIdsPathsLabels.csv N_CHR"; exit; }

# options for nucmer
C_VAL=5000
# location of git repo
FCG=~/FlyComparativeGenomics

# filenames for PlotsForAssemblyCuration.R submission CSVs
ALN_CSV=${IN_CSV%.*}_ALN.csv
SELF_CSV=${IN_CSV%.*}_SELF.csv


# arrays for IDs, files (will make links in working dir), and labels
arr_ID=()
arr_LINK=()
arr_LAB=()

while read line; do 

    ID=$(echo $line | cut -d',' -f1)
    arr_ID+=("$ID")

    PTH=$(echo $line | cut -d',' -f2)
    ln -s $PTH
    LINK=$(basename $PTH)
    arr_LINK+=("$LINK")

    LAB=$(echo $line | cut -d',' -f3)
    arr_LAB+=("$LAB")

done < $IN_CSV

# Between-Assembly Alignments, based on row number in CSV_IN

# primaries
# submit nucmer alignment job
$FCG/RunMummerAlignment-VPGRU.sh \
  -o "${arr_ID[0]}-vs-${arr_ID[3]}" \
  -c $C_VAL \
  ${arr_LINK[0]} \
  ${arr_LINK[3]}
# write line for Rscript submission CSV
echo "p_asm,${arr_ID[0]},${arr_ID[3]},${arr_ID[0]}-vs-${arr_ID[3]}-c${C_VAL}.coords,${arr_LINK[0]}.fai,${arr_LINK[3]}.fai,${arr_LAB[0]},${arr_LAB[3]}" >> $ALN_CSV

# haps against primaries
$FCG/RunMummerAlignment-VPGRU.sh \
  -o "${arr_ID[0]}-vs-${arr_ID[1]}" \
  -c $C_VAL \
  ${arr_LINK[0]} \
  ${arr_LINK[1]}
echo "${arr_ID[1]},${arr_ID[0]},${arr_ID[1]},${arr_ID[0]}-vs-${arr_ID[1]}-c${C_VAL}.coords,${arr_LINK[0]}.fai,${arr_LINK[1]}.fai,${arr_LAB[0]},${arr_LAB[1]}" >> $ALN_CSV

$FCG/RunMummerAlignment-VPGRU.sh \
  -o "${arr_ID[0]}-vs-${arr_ID[2]}" \
  -c $C_VAL \
  ${arr_LINK[0]} \
  ${arr_LINK[2]}
echo "${arr_ID[2]},${arr_ID[0]},${arr_ID[2]},${arr_ID[0]}-vs-${arr_ID[2]}-c${C_VAL}.coords,${arr_LINK[0]}.fai,${arr_LINK[2]}.fai,${arr_LAB[0]},${arr_LAB[2]}" >> $ALN_CSV

$FCG/RunMummerAlignment-VPGRU.sh \
  -o "${arr_ID[3]}-vs-${arr_ID[4]}" \
  -c $C_VAL \
  ${arr_LINK[3]} \
  ${arr_LINK[4]}
echo "${arr_ID[4]},${arr_ID[3]},${arr_ID[4]},${arr_ID[3]}-vs-${arr_ID[4]}-c${C_VAL}.coords,${arr_LINK[3]}.fai,${arr_LINK[4]}.fai,${arr_LAB[3]},${arr_LAB[4]}" >> $ALN_CSV

$FCG/RunMummerAlignment-VPGRU.sh \
  -o "${arr_ID[3]}-vs-${arr_ID[5]}" \
  -c $C_VAL \
  ${arr_LINK[3]} \
  ${arr_LINK[5]}
echo "${arr_ID[5]},${arr_ID[3]},${arr_ID[5]},${arr_ID[3]}-vs-${arr_ID[5]}-c${C_VAL}.coords,${arr_LINK[3]}.fai,${arr_LINK[5]}.fai,${arr_LAB[3]},${arr_LAB[5]}" >> $ALN_CSV


# Self-Alignments
for i in "${!arr_LINK[@]}"; do

    # job to split scaffolds (Chr vs unplaced) and run nucmer
    $FCG/self_align-ChrsVsUnplaced.sh ${arr_LINK[$i]} $N_CHR

    # write line for Rscript submission CSV
    echo "${arr_ID[$i]},${arr_LINK[$i]%.*}-Chrs-vs-unplaced-c${C_VAL}.coords,${arr_LINK[$i]%.*}-Chrs.fa.fai,${arr_LINK[$i]%.*}-unplaced.fa.fai" >> $SELF_CSV

done
