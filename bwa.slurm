#!/bin/bash
#SBATCH -J bwa
#SBATCH -o "bwa.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "bwa.%j.%N.stderr"
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -N 1
#SBATCH -C INTEL

# USAGE: sbatch path/to/bwa.slurm ASSEMBLY INDEX READS.tsv
# INDEX the string for bwa index
# READS.tsv in this format:
# SampleID  path/to/R1.fastq    path/to/R2.fastq

date

module load bwa_mem2/2.2.1
module load samtools/1.17

ASSEMBLY=$1 # assembly file
INDEX=$2 # prefix string for index 
THREADS=32

echo "Mapping reads:"
cat $3
echo "onto assembly $ASSEMBLY with bwa_mem2"
which bwa_mem2
echo -e "\nCalculating depth with samtools"
which samtools

echo "indexing assembly, command:"
echo " bwa-mem2 index -p $INDEX $ASSEMBLY"
bwa-mem2 index -p $INDEX $ASSEMBLY
echo "indexing complete"
date

# function for mapping steps, call with SAMPLE, R1, R2
Map_GetDepths() {
  SAMPLE=$1
  R1=$2
  R2=$3
  OUT="${ASSEMBLY%.f*a}.$SAMPLE"
  echo -e "\nAligning sample $SAMPLE to $INDEX, command:"
  echo " bwa-mem2 mem -t $THREADS $INDEX $R1 $R2 | samtools sort -o $OUT.bam"
  bwa-mem2 mem -t $THREADS $INDEX $R1 $R2 | samtools sort -o $OUT.bam
  echo "Finished writing $OUT.bam"
  echo "Submitting bcfcall.slurm job"
  sbatch ~/FlyComparativeGenomics/bcfcall.slurm $OUT.bam $ASSEMBLY
  echo -e "\nGetting statistics for $OUT.bam, command:\n samtools stats $OUT.bam > $OUT.bam.stats"
  samtools stats $OUT.bam > $OUT.bam.stats
  echo -e "\nCalculating depth in $OUT.bam, command:"
  echo " samtools depth -a $OUT.bam > $OUT.depth.tsv"
  samtools depth -a $OUT.bam > $OUT.depth.tsv
  echo "Writing $OUT.depth_by_scaffold.tsv"
  ~/FlyComparativeGenomics/get_depth_by_scaffold.sh $OUT.depth.tsv
  echo -e "Done with $SAMPLE.\n"
  date
}

while read SAMPLE_LINE; do
    S=$(echo $SAMPLE_LINE | cut -f1)
    R1=$(echo $SAMPLE_LINE | cut -f2)
    R2=$(echo $SAMPLE_LINE | cut -f3)
    Map_GetDepths $S $R1 $R2
done < $3
