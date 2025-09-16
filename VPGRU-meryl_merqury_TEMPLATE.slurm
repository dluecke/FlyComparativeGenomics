#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="merqury"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 1-00:00:00
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                 #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "mqy.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "mqy.%j.%N.stderr"     #optional, prints our standard error

date

module load merqury
which meryl
which merqury.sh
which R

# meryl directory names are seq file names with extension replaced by .meryl
# matches default meryl behavior 
READS_FN=$(basename $READS_FASTQ)
READS_MERYL="${READS_FN%.f*}.meryl"
ASM_FN=$(basename $ASM_FA)
# meryl db names for assemblies, used in -l (reads only) step:
#   meryl difference ${READS_MERYL} ${ASM_MERYL}${ASM2_MERYL} output ${READSONLY}
ASM_MERYL="${ASM_FN%.f*a}.meryl"
if [[ -n $ASM2_FA ]]; then
    ASM2_FN=$(basename $ASM2_FA)
    ASM2_MERYL=" ${ASM2_FN%.f*a}.meryl" # only reads only -l step: meryl difference ${READS_MERYL} ${ASM_MERYL}${ASM2_MERYL}
else
    ASM2_MERYL=""
fi

echo -e "Meryl and Merqury jobs running on -p ceres with -n ${THREADS}\n"

# Meryl kmer counting of HiFi reads
echo -e "Running meryl count on reads:\n ${READS_FASTQ}"
echo -e "with kmer length $K_LEN\n"
echo "CMD: meryl k=${K_LEN} count ${READS_FASTQ} output ${READS_MERYL}"
meryl k=${K_LEN} count ${READS_FASTQ} output ${READS_MERYL}

# Meryl kmer countings on maternal/paternal reads (if provided)
if [[ -n $TRIO ]]; then
    READS_MAT="${ASM_FN%.f*}-mat_reads.meryl"
    READS_PAT="${ASM_FN%.f*}-pat_reads.meryl"
    echo -e "\nTrio mode. Running meryl count on parental reads with kmer length ${K_LEN}."
    echo -e "Maternal:"
    ls $MATERNAL
    echo -e "output prefix ${READS_MAT}\n"
    echo "CMD: meryl k=${K_LEN} count ${MATERNAL} output ${READS_MAT}"
    meryl k=${K_LEN} count ${MATERNAL} output "${READS_MAT}" 
    echo -e "\nPaternal:"
    ls $PATERNAL
    echo -e "output prefix ${READS_PAT}\n"
    echo "CMD: meryl k=${K_LEN} count ${PATERNAL} output ${READS_PAT}"
    meryl k=${K_LEN} count ${PATERNAL} output "${READS_PAT}"

    # Merqury using parental hap-mers 
    # Primary only (no second assembly provided)
    if [[ -z $ASM2_FA ]]; then
        echo -e "\nTrio plus Primary assembly mode (single assembly)."
        echo -e "NOTE: You should have both haplotypes from trio binning - did you mean to include -d ASM2.fasta???"
        echo -e "Running merqury comparing reads meryl output:\n ${READS_MERYL}"
        echo -e "and parental reads:\n ${READS_MAT}\n ${READS_PAT}"
        echo -e "to assembly file:\n ${ASM_FA}"
        echo -e "with output prefix ${MERQURY_OUT}\n"
        echo "CMD: merqury.sh ${READS_MERYL} ${READS_MAT} ${READS_PAT} ${ASM_FA} ${MERQURY_OUT}"
        merqury.sh ${READS_MERYL} ${READS_MAT} ${READS_PAT} ${ASM_FA} ${MERQURY_OUT}
    # Diploid mode, second assembly provided
    else
        echo -e "\nTrio plus Diploid assembly mode (2 assembly files)."
        echo -e "Running merqury comparing reads meryl output:\n ${READS_MERYL}"
        echo -e "and parental reads:\n ${READS_MAT}\n ${READS_PAT}"
        echo -e "to assembly files:\n ${ASM_FA}\n ${ASM2_FA}"
        echo -e "with output prefix ${MERQURY_OUT}\n"
        echo "CMD: merqury.sh ${READS_MERYL} ${READS_MAT} ${READS_PAT} ${ASM_FA} ${ASM2_FA} ${MERQURY_OUT}"
        merqury.sh ${READS_MERYL} ${READS_MAT} ${READS_PAT} ${ASM_FA} ${ASM2_FA} ${MERQURY_OUT}
    fi
else
    # Merqury without parental hap-mers
    # Primary only (no second assembly provided)
    if [[ -z $ASM2_FA ]]; then
        echo -e "\nPrimary assembly mode (single assembly)"
        echo -e "Running merqury comparing reads meryl output:\n ${READS_MERYL}"
        echo -e "to assembly file:\n ${ASM_FA}"
        echo -e "with output prefix ${MERQURY_OUT}\n"
        echo "CMD: merqury.sh ${READS_MERYL} ${ASM_FA} ${MERQURY_OUT}"
        merqury.sh ${READS_MERYL} ${ASM_FA} ${MERQURY_OUT}
    # Diploid mode, second assembly provided
    else
        echo -e "\nDiploid assembly mode (2 assembly files)"
        echo -e "Running merqury comparing reads meryl output:\n ${READS_MERYL}"
        echo -e "to assembly files:\n ${ASM_FA}\n ${ASM2_FA}"
        echo -e "with output prefix ${MERQURY_OUT}\n"
        echo "CMD: merqury.sh ${READS_MERYL} ${ASM_FA} ${ASM2_FA} ${MERQURY_OUT}"
        merqury.sh ${READS_MERYL} ${ASM_FA} ${ASM2_FA} ${MERQURY_OUT}
    fi
fi

# Replaced by merqury.sh command with parental reads
# if [[ -n $TRIO ]]; then
#     echo -e "Running merqury hapmers.sh on maternal, paternal, and F1 meryl dbs\n"
#     echo "CMD: $MERQURY/trio/hapmers.sh ${READS_MAT} ${READS_PAT} ${READS_MERYL}"
#     $MERQURY/trio/hapmers.sh ${READS_MAT} ${READS_PAT} ${READS_MERYL}

# fi

if [[ -n $READS_ONLY ]]; then
    # ReadsOnly outputs names for fastq file name
    READSONLY_MERYL="${READS_FN%.f*}.ReadsOnly.meryl"
    READSONLY_FASTQ="${READS_FN%.f*}.ReadsOnly.fastq.gz"
    READSONLY_NAMES="${READS_FN%.f*}.ReadsOnly.names.txt"

    echo -e "Running meryl difference to find kmers observed only in hifi reads:\n ${READSONLY_MERYL}\n"
    echo "CMD: meryl difference ${READS_MERYL} ${ASM_MERYL}${ASM2_MERYL} output ${READSONLY_MERYL}"
    meryl difference ${READS_MERYL} ${ASM_MERYL}${ASM2_MERYL} output ${READSONLY_MERYL}

    echo -e "\nRunning meryl-lookup to extract fastq.gz file of reads with kmers not observed in assembly:\n ${READSONLY_FASTQ}\n"
    echo "CMD: meryl-lookup -include -sequence ${READS_FASTQ} -mers ${READSONLY_MERYL} | pigz -c > ${READSONLY_FASTQ}"
    meryl-lookup -include -sequence ${READS_FASTQ} -mers ${READSONLY_MERYL} | \
        pigz -c > ${READSONLY_FASTQ}

    echo -e "\nExtracting read names from ${READSONLY_FASTQ}:\n ${READSONLY_NAMES}\n"
    echo "CMD: zgrep -e "^@" ${READSONLY_FASTQ} > ${READSONLY_NAMES}"
    zgrep -e "^@" ${READSONLY_FASTQ} > ${READSONLY_NAMES}
fi


date
