#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="tblastx"
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p short
#SBATCH -o "tbx.stdout.%j.%N"           # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "tbx.stderr.%j.%N"           #optional but it prints our standard error
#SBATCH --mail-user=david.luecke@usda.gov
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

module load blast+/2.15.0

TARGET=$1
QUERY=$2

if [ ! -f $TARGET.ndb ]; then
    makeblastdb -dbtype nucl -in $TARGET
fi

tblastx -num_threads 32 \
    -db $TARGET \
    -query $QUERY \
    -outfmt 6 \
    -out ${QUERY%.f*a}-${TARGET%.f*a}-tblastx.out6
