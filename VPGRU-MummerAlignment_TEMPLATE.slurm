#!/bin/bash
#SBATCH -A vpgru
#SBATCH -J mummer
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -p ceres
#SBATCH -t 6:00:00 
#SBATCH -o "mummer.stdout.%j.%N"          # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "mummer.stderr.%j.%N"          #optional but it prints our standard error
#SBATCH --mail-user=David.Luecke@usda.gov   # going to replace with $USER in the submission shell
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date

module load mummer/4.0.0rc1
module load samtools/1.17 # just need for faidx

echo "mummer aligner loaded, version:"
mummer --version
echo "samtools loaded for faidx to get sequence lengths, module:"
which samtools

echo -e "\nPerforming nucmer alignment ${ALIGNMENT_NAME}:"
echo "  Reference: ${REFERENCE}"
echo "  Query: ${QUERY}"
echo "  Seed Match Length: ${C_VAL}"
echo "  Job Name: ${RUN_ID}-mummer"
echo "  Partition:  ${PARTITION}"
echo "  Cores/Threads:  ${THREADS}"
echo -e "Command:\n nucmer -p ${ALIGNMENT_NAME} --threads=${THREADS} ${MAXMATCH} -c ${C_VAL} ${REFERENCE} ${QUERY}"

nucmer -p ${ALIGNMENT_NAME} --threads=${THREADS} ${MAXMATCH} -c ${C_VAL} ${REFERENCE} ${QUERY}

echo "Finished alignment, using show-coords and mummerplot to generate alignment coordinates"
echo -e "Command:\n show-coords ${ALIGNMENT_NAME}.delta > ${ALIGNMENT_NAME}.coords"
show-coords ${ALIGNMENT_NAME}.delta > ${ALIGNMENT_NAME}.coords
echo -e "Command:\n mummerplot -p ${ALIGNMENT_NAME} -R ${REFERENCE} -Q ${QUERY} --postscript ${ALIGNMENT_NAME}.delta"
mummerplot -p ${ALIGNMENT_NAME} -R ${REFERENCE} -Q ${QUERY} --postscript ${ALIGNMENT_NAME}.delta
echo "Finished mummerplot. Rough draft dotplot file: ${ALIGNMENT_NAME}.ps"

samtools faidx $REFERENCE
samtools faidx $QUERY

echo "Use ${ALIGNMENT_NAME}.coords, $REFERENCE.fai, $QUERY.fai files\
 with make_l.coords() and PlotDFCoords() functions in alignment_region_dotplot-FUNCTIONS.R"

date
