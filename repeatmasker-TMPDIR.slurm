#!/bin/bash
#SBATCH -J repeatMask-TMPDIR
#SBATCH -o "rm-%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "rm-%j.%N.stderr"
#SBATCH -c 32
#SBATCH -p medium
#SBATCH -N 1

date
module load repeatmodeler/2.0.4
module load repeatmasker/4.1.5

ASM=$1
DB_NAME=$2
WORKDIR=$PWD

# copy of script with only repeatmasker steps, repeatmodeler commented out
#echo "Running repeatmodeler on assembly $1 with output $2"
#which RepeatModeler

echo "Running in TMPDIR: ${TMPDIR} to avoid I/O limit"
cp $ASM ${TMPDIR}
cp RM_*/consensi.fa.classified ${TMPDIR}
cd ${TMPDIR}

#echo "command:"
#echo "BuildDatabase -name $DB_NAME $ASM"
#BuildDatabase -name $DB_NAME $ASM

#echo "command:"
#echo "RepeatModeler -threads 32 -database $DB_NAME"
#RepeatModeler -threads 32 -database $DB_NAME
#date

echo "Running RepeatMasker on $ASM with library RM_*/consensi.fa.classified"
which RepeatMasker

#cp RM_*/consensi.fa.classified .
echo "command:"
echo "RepeatMasker -pa 32 -xsmall -dir repmask_RepModLib_${DB_NAME} -e rmblast -lib consensi.fa.classified $ASM"
RepeatMasker -pa 32 -xsmall -dir repmask_RepModLib_${DB_NAME} -e rmblast -lib consensi.fa.classified $ASM

date

echo "Moving results from $TMPDIR to $WORKDIR"
echo "commands:"
echo " rsync -rv ${TMPDIR}/repmask_RepModLib_${DB_NAME} $WORKDIR"
#echo " rsync -rv ${TMPDIR}/RM* $WORKDIR"
echo " rsync -rv ${TMPDIR}/${DB_NAME}* $WORKDIR"
rsync -rv ${TMPDIR}/repmask_RepModLib_${DB_NAME} $WORKDIR
#rsync -rv ${TMPDIR}/RM* $WORKDIR
rsync -rv ${TMPDIR}/${DB_NAME}* $WORKDIR

date
