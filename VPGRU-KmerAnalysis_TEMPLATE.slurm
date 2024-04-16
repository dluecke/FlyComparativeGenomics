#!/bin/bash
#SBATCH -J KmerAnalysis
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -o "KmerAnalysis.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "KmerAnalysis.stderr.%j.%N" 		#optional but it prints our standard error
#SBATCH --mail-user=David.Luecke@usda.gov   # going to replace with $USER in the submission shell
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date

module load samtools
module load jellyfish2
jellyfish --version
$SOFTWARE/genomescope2.0/genomescope.R --version

FASTQ="${BAMFILE%%.bam}.fastq"

echo -e "\nConverting BAM file:\n ${BAMFILE}\nto FASTQ file:\n ${FASTQ}"
samtools fastq $BAMFILE > $FASTQ

SIZE=$(jellyfish mem -m $KLEN --mem 96G | cut -f1 -d' ')
echo -e "\nStarting jellyfish count with -m $KLEN and -s $SIZE, output into:\n ${PWD}/${RUNID}.jf"
jellyfish count -C -m $KLEN -s $SIZE -t 32 $FASTQ -o ${RUNID}.jf
echo -e "\nStarting jellyfish histo to write:\n ${PWD}/${RUNID}.histo"
jellyfish histo -t 32 ${RUNID}.jf > ${RUNID}.histo

echo -e "\nStarting genomescope analysis on ${RUNID}.histo with -k ${KLEN} and -p 2."
echo -e "Output in:\n ${PWD}/${RUNID}\n"
$SOFTWARE/genomescope2.0/genomescope.R -i ${RUNID}.histo -o $RUNID -k $KLEN -p 2 --fitted_hist

date



