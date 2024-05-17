#!/bin/bash
#SBATCH -J KmerAnalysis
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -p short
#SBATCH -o "KmerAnalysis.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "KmerAnalysis.stderr.%j.%N" 		#optional but it prints our standard error
#SBATCH --mail-user=David.Luecke@usda.gov   # going to replace with $USER in the submission shell
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date

module load samtools
module load jellyfish2
jellyfish --version
$SOFTWARE/genomescope2.0/genomescope.R --version

if [[ "$SEQFILE" == *".bam" ]]; then
    FASTQ="${SEQFILE%%.bam}.fastq"
    if [ ! -f $FASTQ ]; then
        echo -e "\nConverting BAM file:\n ${SEQFILE}\nto FASTQ file:\n ${FASTQ}\n"
        samtools fastq --threads $((THREADS-1)) $SEQFILE > $FASTQ
    else
        echo -e "\nUsing pre-existing FASTQ file:\n ${FASTQ}\n"
    fi
elif [[ "$SEQFILE" == *".fq" || "$SEQFILE" == *".fastq" ]]; then
    FASTQ=$SEQFILE
    echo -e "\nUsing provided FASTQ file:\n ${FASTQ}\n"
else
    echo "file extension not recognized: $SEQFILE"
    exit;
fi

# estimating jellyfish count --size parameter, based on default 32 cores - jellyfish2 automatically adjusts
SIZE=$(jellyfish mem -m $KLEN --mem 96G | cut -f1 -d' ')

echo -e "Starting jellyfish count with -m $KLEN and -s $SIZE, output into:\n ${PWD}/${RUNID}.jf\n"
jellyfish count -C -m $KLEN -s $SIZE -t ${THREADS} $FASTQ -o ${RUNID}.jf
echo -e "Starting jellyfish histo to write:\n ${PWD}/${RUNID}.histo\n"
jellyfish histo -t ${THREADS} ${RUNID}.jf > ${RUNID}.histo

echo -e "Starting genomescope analysis on ${RUNID}.histo with -k ${KLEN} and -p 2."
echo -e "Output in:\n ${PWD}/${RUNID}\n"
$SOFTWARE/genomescope2.0/genomescope.R -i ${RUNID}.histo -o $RUNID -k $KLEN -p 2 --fitted_hist

# Write line for csv containing sequence file identifier and estimated error rate
sleep 10
TrayID=$(echo $FASTQ | cut -d'.' -f1)
EstErr=$(grep err: KmerAnalysis-${RUNID}.stdout | cut -d' ' -f5 | cut -d':' -f2)
echo "${TrayID},${EstErr}" > ${RUNID}_Line-HiFiTray_EstError.csv
echo -e "\nWrote ${FASTQ} line for HiFiTray_EstError.csv:\n ${RUNID}_Line-HiFiTray_EstError.csv\n"

date

