#!/bin/bash
#SBATCH -J mummer
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -p short 
#SBATCH -o "mummer.stdout.%j.%N"          # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "mummer.stderr.%j.%N"          #optional but it prints our standard error
#SBATCH --mail-user=David.Luecke@usda.gov   # going to replace with $USER in the submission shell
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

date

module load mummer/4.0.0rc1
echo "mummer aligner loaded, version:"
mummer --version

echo "Performing nucmer alignment ${ALIGNMENT_NAME}:"
echo "  Reference: ${REFERENCE}"
echo "  Query: ${QUERY}"
echo "  Seed Match Length: ${C_VAL}"
echo "  Job Name: ${RUN_ID}-mummer"
echo "  Partition:  ${PARTITION}"
echo "  Cores/Threads:  ${N_CORES}

nucmer -p ${ALIGNMENT_NAME} --threads=${THREADS} -c ${C_VAL} ${REFERENCE} ${QUERY}

echo "Finished alignment, using mummerplot to generate alignment coordinates"
mummerplot -p ${ALIGNMENT_NAME} -R ${REFERENCE} -Q ${QUERY} --postscript ${ALIGNMENT_NAME}.delta
echo "Finished mummerplot. Rough draft dotplot file: ${ALIGNMENT_NAME}.ps"

echo "Extracting alignment coordinates into tsv format"
${TOOLS_PATH}/alignment_and_visualization/convert_gnuplot_to_tsv.sh ${ALIGNMENT_NAME}
echo "Finished extracting tsv files for use by alignment_region_dotplot.R:"
echo "  ${ALIGNMENT_NAME}.plot.tsv"
echo "  ${ALIGNMENT_NAME}.breaks.tsv"
echo "------------------------"

date
