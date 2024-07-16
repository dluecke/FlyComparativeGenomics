#!/bin/bash
#SBATCH --job-name="pbmm2"   #name of this job
#SBATCH -p mem              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "stderr.%j.%N"     #optional, prints our standard error

date

# Print info to stdout for reference
echo -e "Starting pbmm2 to run minimap2, aligning reads in:\n ${HIFI_FOFN}"
echo -e "onto assembly:\n ${ASM_FASTA}\nRun in :\n $PWD"
echo -e "with output BAM file ${OUT_PREFIX}.bam"
echo -e "using pbmm2 conda environment in /project/vpgru/.conda/pbmm2\n"

echo -e "\nRan on the ${SLURM_QUEUE} partition with $N_THREAD cores"
echo -e "Job name:\n pbmm2-${OUT_PREFIX}\n"

pbmm2 --version
echo -e "\nCommand:\n pbmm2 align --preset HIFI --log-level INFO --sort ${ASM_FASTA} ${HIFI_FOFN} ${OUT_PREFIX}.bam\n"
pbmm2 align --preset HIFI --log-level INFO --sort ${ASM_FASTA} ${HIFI_FOFN} ${OUT_PREFIX}.bam

if [[ -n $SAM_OUT || -n $PAF_OUT ]]; then
    echo -e "\nUsing samtools to produce SAM file ${OUT_PREFIX}.sam"
    module load samtools
    samtools --version
    echo -e "\nCommand:\n samtools view -h ${OUT_PREFIX}.bam > ${OUT_PREFIX}.sam\n"
    samtools view -h ${OUT_PREFIX}.bam > ${OUT_PREFIX}.sam
fi

if [[ -n $PAF_OUT ]]; then
    echo -e "\nUsing minimap2 paftools.js to produce ${OUT_PREFIX}.paf.gz"
    module load minimap2
    echo -e "minimap2 version:"
    minimap2 --version
    echo -e "\nCommand:\n paftools.js sam2paf ${OUT_PREFIX}.sam | gzip -c - > ${OUT_PREFIX}.paf.gz"
    paftools.js sam2paf ${OUT_PREFIX}.sam | gzip -c - > ${OUT_PREFIX}.paf.gz
fi

date
