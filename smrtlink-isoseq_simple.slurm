#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="isoseq"
#SBATCH -N 1
#SBATCH -n 17
#SBATCH -p ceres
#SBATCH -t 16:00:00
#SBATCH -o "isoseq.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "isoseq.%j.%N.stderr"               #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# takes demuxed isoseq BAM and extracts transcripts fasta
# USAGE: sbatch smrtlink-isoseq_simple.slurm SAMPLE.BAM [FILE_PREFIX, default bam name without extension]

# written for older isoseq data (SequelII, 2023), newer Kinnex data takes a different approach

module load smrtlink/25.3.0.273777

IN_BAM=$1
if [ -z "$2" ]; then
    FILE_PREFIX=$(basename "$IN_BAM" .bam)
else
    FILE_PREFIX=$2
fi

date
echo "Extracting transcripts from $IN_BAM, output prefix: $FILE_PREFIX"

lima --version
isoseq --version

# get isoseq primers (unclear if these are precisely correct, but seem to work for required primers in isoseq refine step)
if [ ! -f IsoSeq_v2_primers_12.fasta ]; then
    wget https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta
fi

# write primer fasta, empirically seen in SequelIIe_r64312e_20230809_195524_CmacIsoSeq data
if [ ! -f IsoSeqX_3p.fa ]; then
    echo -e ">IsoSeqX_3p\AAGCAGTGGTATCAACGCAGAGTAC" > IsoSeqX_3p.fa
fi

# lima to remove primers
echo "running lima to remove primers"
echo "CMD: lima --num-threads 16 $IN_BAM IsoSeqX_3p.fa ${FILE_PREFIX}_trimmed.bam"
lima --num-threads 16 $IN_BAM IsoSeqX_3p.fa ${FILE_PREFIX}_trimmed.bam

# isoseq refine to get full-length non-chimeric reads
echo "running isoseq refine to get full-length non-chimeric reads"
echo "CMD: isoseq refine --num-threads 16 ${FILE_PREFIX}_trimmed.bam IsoSeq_v2_primers_12.fasta ${FILE_PREFIX}_refined.bam --require-polya"
isoseq refine --num-threads 16 ${FILE_PREFIX}_trimmed.bam IsoSeq_v2_primers_12.fasta ${FILE_PREFIX}_refined.bam --require-polya

# isoseq cluster to get high-quality isoforms, including fasta output
echo "running isoseq cluster2 to get high-quality isoforms"
echo "CMD: isoseq cluster2 --singletons ${FILE_PREFIX}_refined.bam ${FILE_PREFIX}.bam ${FILE_PREFIX}.fasta"
isoseq cluster2 --singletons ${FILE_PREFIX}_refined.bam ${FILE_PREFIX}.bam ${FILE_PREFIX}.fasta

date