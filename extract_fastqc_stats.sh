#!/bin/bash

# extract_fastqc_stats.sh writes a csv file with all summary statistics and test results
# for all fastqc results in working directory (unzipped *_fastqc/ output directories)
# output csv CURRENTDIRECTORY-fastqc_stats.csv
# written for output from FastQC v0.12.1

# output file name
CSV_OUT="$(basename $PWD)-fastqc_stats.csv"

# header row for output
echo "Filename,FileType,Encoding,TotalSequences,TotalBases,\
FlaggedSequences,SequenceLength,PercentGC,\
BasicStatistics,PerBaseQual,PerTileQual,PerSequenceQualScores,\
PerBaseSequenceContent,PerSequenceGcContent,PerBaseNContent,\
SequenceLengthDistribution,SequenceDuplicationLevels,\
OverrepresentedSequences,AdapterContent" > $CSV_OUT

# all fastqc results directories
FASTQC_DIRS=($(ls -d *fastqc))

for FQC_D in ${FASTQC_DIRS[@]}; do

    # array for data from fastqc_data.txt, includes filename, need to remove spaces
    CSV_ROW=($(head -n11 ${FQC_D}/fastqc_data.txt | tail -n8 | cut -f2 | tr -d ' '))

    # append data from summary.txt, PASS/WARN/FAIL test results
    CSV_ROW+=($(cut -f1 ${FQC_D}/summary.txt))

    # print row array as comma separated, replace final comma with line break
    printf '%s,' "${CSV_ROW[@]}" | sed 's/,$/\n/'

done >> $CSV_OUT
