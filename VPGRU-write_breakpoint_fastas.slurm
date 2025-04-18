#!/bin/bash
#SBATCH --job-name="breakpoint_fastas"
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p short
#SBATCH -o "bpfa.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "bpfa.%j.%N.stderr"
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# Script to take assembly, blast hits, read alignments and fastqs and write fastas for breakpoint region alignments
# by David Luecke 4/18/25
# Intended for nuMt verification, blast hits from mitochondrial sequence and aligned HiFi reads, both vs same reference
# outputs fastas:
#   assembly region centered on cluster of blast hits
#   reads which span either breakpoint between hit cluster and rest of assembly
# uses scripts from FlyComparativeGenomics/ git repo: get_reads_spanning_BLAST_hits.sh (wraps get_spanning_reads.sh)
# requires samtools and seqkit

module load samtools/1.17
module load seqkit/2.4.0

# Usage: 
# sbatch [--export=OPTIONAL=$VAL,X="Y",Z=123] VPGRU-write_breakpoint_fastas.slurm ASM.fa BLAST.out6 ALN.bam READS.fq

# SET VARIABLES
# location of other scripts same as for this script in FlyComparativeGenomics repo
FCG_PATH=$(dirname $0)

# Required variables, given as positional 
[[ -n $1 ]] && REF_ASM=$1 || (echo "Nothing in position 1"; exit; )
[[ -n $2 ]] && BLAST_HITS=$2 || (echo "Nothing in position 2"; exit; )
[[ -n $3 ]] && ALN_BAM=$3 || (echo "Nothing in position 3"; exit; )
[[ -n $4 ]] && READS_FQ=$4 || (echo "Nothing in position 4"; exit; )

# Optional variables, can be set with sbatch --export on submission
# minimum space between blast hit clusters (bp)
[[ -n $MIN_CLUSTER_GAP ]] || MIN_CLUSTER_GAP=50000
# threshold score for a hit in a cluster (bit score)
[[ -n $SEED_BIT_THRESH ]] || SEED_BIT_THRESH=5000
# distance on either end of cluster boundary aligned reads must extend (bp)
[[ -n $READS_BOUNDARY_SPAN ]] || READS_BOUNDARY_SPAN=5000
# total length of region centered on hits to extract from assembly (bp)
#  if cluster is closer to end of sequence than half this length will truncate (instead of extending to full length in other direction)
[[ -n $TOTAL_REGION_LENGTH ]] || TOTAL_REGION_LENGTH=200000

# HIT CLUSTERING
$FCG_PATH/get_reads_spanning_BLAST_hits.sh $BLAST_HITS $ALN_BAM $MIN_CLUSTER_GAP $SEED_BIT_THRESH $READS_BOUNDARY_SPAN

# FASTAS FOR EACH CLUSTER
for CLUSTERFILE in HitClusterInfo-*txt; do

    CLUSTER_ID=$(echo $CLUSTERFILE | cut -d'-' -f2-)
    CLUSTER_ID=${CLUSTER_ID%.*}
    mkdir $CLUSTER_ID-fastas

    # Assembly Sequence Region
    # sequence name
    CLUSTER_SEQ=$(grep "sequence name:" $CLUSTERFILE | cut -f2)
    # cluster coordinates
    CLUSTER_BEG=$(grep "region coordinates:" $CLUSTERFILE | cut -f2 | cut -d' ' -f1)
    CLUSTER_END=$(grep "region coordinates:" $CLUSTERFILE | cut -f2 | cut -d' ' -f3)
    CLUSTER_MID=$(( (CLUSTER_BEG + CLUSTER_END) / 2 ))
    # region coordinates
    REGION_BEG=$(( CLUSTER_MID - TOTAL_REGION_LENGTH/2 ))
    REGION_BEG=$(( REGION_BEG > 0 ? REGION_BEG : 0 )) # need 0 floor for start coord, samtools faidx won't take negative coordinate
    REGION_END=$(( CLUSTER_MID + TOTAL_REGION_LENGTH/2 )) # ok if upper end higher than seq length, samtools faidx will truncate
    # get region fasta with samtools faidx
    samtools faidx $REF_ASM ${CLUSTER_SEQ}:${REGION_BEG}-${REGION_END} \
        > $CLUSTER_ID-fastas/${REF_ASM%.*}-${CLUSTER_SEQ}_${REGION_BEG}to${REGION_END}.fa

    # Cluster Boundary Spanning Reads
    # filename with read id list 
    READS_LIST=$(grep "list of read IDs spanning either boundary:" $CLUSTERFILE | cut -f2)
    # write fasta of reads with seqkit grep and fq2fa
    seqkit grep -f $READS_LIST $READS_FQ | seqkit fq2fa \
        > $CLUSTER_ID-fastas/${READS_FQ%.*}-${CLUSTER_ID}_span${READS_BOUNDARY_SPAN}bp.fa

done
