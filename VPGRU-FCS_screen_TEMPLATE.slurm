#!/bin/bash
#SBATCH --job-name="FCS-VPGRU"   #name of this job
#SBATCH -p mem              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "stderr.%j.%N"     #optional, prints our standard error

date

module load fcs-genome

if [[ "$IN_SEQFILE" == *".bam" ]]; then
    module load samtools
    SEQFILE="${IN_SEQFILE%%.bam}.fasta"
    echo -e "\nConverting BAM file:\n ${IN_SEQFILE}\nto FASTA file:\n ${SEQFILE}\n"
    which samtools
    samtools fasta --threads $((THREADS-1)) $IN_SEQFILE > $SEQFILE
elif [[ "$IN_SEQFILE" == *".fq" || "$IN_SEQFILE" == *".fastq" ]]; then
    module load seqkit
    SEQFILE="${IN_SEQFILE%%.f*q}.fasta"
    echo -e "\nConverting FASTQ file:\n ${IN_SEQFILE}\nto FASTA file:\n ${SEQFILE}\n"
    which seqkit
    seqkit --threads $THREADS fq2fa $IN_SEQFILE -o $SEQFILE
elif [[ "$IN_SEQFILE" == *".fa" || "$IN_SEQFILE" == *".fasta" ]]; then
    SEQFILE=$IN_SEQFILE
else
    echo "file extension not recognized: $IN_SEQFILE"
    exit;
fi

echo -e "\nRunning FCS-genome screen on file:\n ${SEQFILE}"
echo -e "using taxonomic ID:\n $TAXID"
echo -e "with output directed to:\n $OUTDIR"
echo -e "\nRun on the mem partition with $THREADS tasks\n"

which fcs.py
echo 
fcs.py screen genome --fasta $SEQFILE --out-dir $OUTDIR --tax-id $TAXID

echo "FCS-genome screen complete.  Collecting taxa summaries."
${FCG_REPO}/collect_FCS_taxa_summaries.sh $OUTDIR 

date