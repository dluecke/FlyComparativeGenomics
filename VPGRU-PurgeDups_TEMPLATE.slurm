#!/bin/bash
#SBATCH --job-name="purge_dups"   #name of this job
#SBATCH -p ceres              #name of the partition (queue) you are submitting to
#SBATCH -t 6:00:00
#SBATCH -N 1                  #number of nodes in this job
#SBATCH -n 32                #number of cores/tasks in this job, you get all 20 physical cores with 2 threads per core with hyper-threading
#SBATCH --mail-user=David.Luecke@usda.gov  #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH -o "pd.%j.%N.stdout"     # standard output, %j adds job number to output file name and %N adds the node name
#SBATCH -e "pd.%j.%N.stderr"     #optional, prints our standard error

date

module load purge_dups/1.2.6
module load minimap2

# Print info to stdout
echo -e "Starting purge_dups on assembly:\n ${PRI_ASM}"
echo -e "using HiFi alignment:\n ${PAF_IN}\nRun in :\n $PWD"
echo -e "with output file prefix:\n ${PD_OUT}"
echo -e "using purge_dups and minimap2 modules\n"
echo -e "Submitting to $QUEUE partition with $N_THREAD threads"

purge_dups -h 2>&1 | head -n5
echo -e "minimap2 version: $(minimap2 --version)\n"

mkdir $PD_OUT

# name of split fasta for self-alignment (avoids writing to same location as PRI_ASM if given as path)
ASM_SPLIT="$(basename $PRI_ASM).split"

echo "Running initial purge_dups steps"
echo -e "Commands:\n pbcstat -O $PD_OUT $PAF_IN"
echo -e " calcuts $PD_OUT/PB.stat > cutoffs.$PD_OUT 2>calcuts.$PD_OUT.log"
pbcstat -O $PD_OUT $PAF_IN
calcuts $PD_OUT/PB.stat > cutoffs.$PD_OUT 2>calcuts.$PD_OUT.log

if [[ -n $M_CUT ]]; then
    # get -l and -u values from auto-cutoff file
    CUT_VALS=($(awk '{print $1, $NF}' cutoffs.$PD_OUT))
    echo "Finding manual cutoff thresholds"
    echo -e "command:\n calcuts -l ${CUT_VALS[0]} -m $M_CUT -u ${CUT_VALS[1]} $PD_OUT/PB.stat \
        > cutoffs_manual.$PD_OUT 2>calcuts_manual.$PD_OUT.log"
    calcuts -l ${CUT_VALS[0]} -m $M_CUT -u ${CUT_VALS[1]} $PD_OUT/PB.stat \
        > cutoffs_manual.$PD_OUT 2>calcuts_manual.$PD_OUT.log
    CUTOFFS=cutoffs_manual.$PD_OUT
else
    CUTOFFS=cutoffs.$PD_OUT
fi

if [[ -f $ASM_SPLIT.self.paf.gz ]]; then
    echo "Using previously generated self-alignment file $ASM_SPLIT.self.paf.gz"
else
    echo "Starting assembly self-alignment"
    echo -e "Commands:\n split_fa $PRI_ASM > $ASM_SPLIT"
    echo " minimap2 -xasm5 -DP $ASM_SPLIT $ASM_SPLIT | gzip -c - > $ASM_SPLIT.self.paf.gz"
    split_fa $PRI_ASM > $ASM_SPLIT
    minimap2 -xasm5 -DP $ASM_SPLIT $ASM_SPLIT | \
        gzip -c - > $ASM_SPLIT.self.paf.gz
fi

echo "Running main purge_dups step"
echo -e "Command:\n purge_dups -2 -T $CUTOFFS -c $PD_OUT/PB.base.cov $ASM_SPLIT.self.paf.gz \
    > dups.$PD_OUT.bed 2> purge_dups.$PD_OUT.log"
purge_dups -2 -T $CUTOFFS -c $PD_OUT/PB.base.cov \
    $ASM_SPLIT.self.paf.gz > \
    dups.$PD_OUT.bed 2> purge_dups.$PD_OUT.log

echo "Purging assembly"
echo -e "Command:\n get_seqs -p $PD_OUT -e dups.$PD_OUT.bed $PRI_ASM"
get_seqs -p $PD_OUT -e dups.$PD_OUT.bed $PRI_ASM

echo 
date
