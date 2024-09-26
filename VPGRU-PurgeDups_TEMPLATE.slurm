#!/bin/bash
#SBATCH --job-name="purge_dups"   #name of this job
#SBATCH -p medium              #name of the partition (queue) you are submitting to
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "pd.stdout.%j.%N"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "pd.stderr.%j.%N"     #optional, prints our standard error

date

module load purge_dups
module load minimap2

# Print info to stdout
echo -e "Starting purge_dups on assembly:\n ${PRI_ASM}"
echo -e "using HiFi alignment:\n ${PAF_IN}\nRun in :\n $PWD"
echo -e "with output file prefix:\n ${PD_OUT}"
echo -e "using purge_dups and minimap2 modules\n"
echo -e "Submitting to $QUEUE partition with $N_THREAD threads"

mkdir $PD_OUT

echo "Running initial purge_dups steps"
echo -e "Commands:\n pbcstat -O $PD_OUT $PAF_IN"
echo -e " calcuts $PD_OUT/PB.stat > cutoffs.$PD_OUT 2>calcuts.$PD_OUT.log"
pbcstat -O $PD_OUT $PAF_IN
calcuts $PD_OUT/PB.stat > cutoffs.$PD_OUT 2>calcuts.$PD_OUT.log

echo "Starting assembly self-alignment"
echo -e "Commands:\n split_fa $PRI_ASM > $PRI_ASM.split"
echo " minimap2 -xasm5 -DP $PRI_ASM.split $PRI_ASM.split | gzip -c - > $PRI_ASM.split.self.paf.gz"
split_fa $PRI_ASM > $PRI_ASM.split
minimap2 -xasm5 -DP $PRI_ASM.split $PRI_ASM.split | \
    gzip -c - > $PRI_ASM.split.self.paf.gz

echo "Running main purge_dups step"
echo -e "Command:\n purge_dups -2 -T cutoffs.$PD_OUT -c $PD_OUT/PB.base.cov $PRI_ASM.split.self.paf.gz > dups.$PD_OUT.bed 2> purge_dups.$PD_OUT.log"
purge_dups -2 -T cutoffs.$PD_OUT -c $PD_OUT/PB.base.cov \
    $PRI_ASM.split.self.paf.gz > \
    dups.$PD_OUT.bed 2> purge_dups.$PD_OUT.log

echo "Purging assembly"
echo -e "Command:\n get_seqs -p $PD_OUT -e dups.$PD_OUT.bed $PRI_ASM"
get_seqs -p $PD_OUT -e dups.$PD_OUT.bed $PRI_ASM

echo 
date
