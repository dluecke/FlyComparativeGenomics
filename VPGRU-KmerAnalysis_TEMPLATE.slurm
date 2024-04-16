#!/bin/bash
#SBATCH -J KmerAnalysis-$RUNID
#SBATCH -N 1
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -o "otb.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "otb.stderr.%j.%N" 		#optional but it prints our standard error
#SBATCH --mail-user=${USER}@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date

module load samtools
module load jellyfish2
jellyfish --version
$SOFTWARE/genomescope2.0/genomescope.R --version

FASTQ="${BAMFILE%%.bam}.fastq"

echo "samtools fastq $BAMFILE > $FASTQ"

SIZE=$(jellyfish mem -m $KLEN --mem 96G | cut -f1 -d' ')
echo "jellyfish count -C -m $KLEN -s $SIZE -t 32 $FASTQ -o ${RUNID}.jf"
echo "jellyfish histo -t 32 ${RUNID}.jf > ${RUNID}.histo"

echo "SIZE is $SIZE"
jellyfish mem -m 21 -s 32G

echo "$SOFTWARE/genomescope2.0/genomescope.R -i ${RUNID}.histo -o $RUNID -k $KLEN -p 2 --fitted_hist"



date



