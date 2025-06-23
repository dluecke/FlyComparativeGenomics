#!/bin/bash

# get_pct_chrs.sh reports percent of total assmebly length in chromosome-scale scaffolds
# takes fasta index FAI of assembly and optional MIN_LEN_CHR value

# check first argument is fasta index *.fai
if [[ "$1" != *".fai" ]]; then
    echo "USAGE: get_pct_chrs.sh ASM.fa.fai [MIN_LEN_CHR(in bp, default 1e7)]"
    exit
fi

FAI=$1

# default min chromosome lenth to 10Mbp if not provided
[[ -n $2 ]] && MIN_LEN_CHR=$2 || MIN_LEN_CHR=10000000

# array of total chromosome length, total assembly length, and pct in chromosomes
ARR_CHRS=($(awk -v chrlen=$MIN_LEN_CHR '
        {totlen += $2} 
        $2 > chrlen {totchr += $2; nchr += 1} 
    END {
        print totchr, totlen, totchr/totlen, nchr
    }
    ' $FAI))

# print output to screen
echo "Calculating percent of total assembly length in chromosome-scale sequences"
echo -e "Assembly index:\n $FAI"
echo -e "Chromosome length threshold:\n $MIN_LEN_CHR\n"
echo "All sequences in $FAI longer than ${MIN_LEN_CHR} bp:"
awk -v chrlen=$MIN_LEN_CHR '$2>chrlen {print " "$1"\t"$2}' $FAI

echo -e "\nAssembly statistics (can grep and cut -f2 to extract values):"
echo -e "AssemblyFAI:\t$FAI"
echo -e "MinLenChrom:\t$MIN_LEN_CHR"
echo -e "TotalLength:\t${ARR_CHRS[1]}"
echo -e "ChromLength:\t${ARR_CHRS[0]}"
echo -e "ChromNumber:\t${ARR_CHRS[3]}"
echo -e "PctInChroms:\t${ARR_CHRS[2]}"
