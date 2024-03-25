#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p long
#SBATCH -o "otb.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "otb.stderr.%j.%N" 		#optional but it prints our standard error
#SBATCH --mail-user=<username>
#SBATCH --mail-type=START,END,FAIL               #will receive an email when job ends or fails

# Modified by DL to incorporate troubleshoot suggestions from DM while OTB still in beta (March 2024)
# more info at https://github.com/molikd/otb/wiki

module load nextflow/22.04.3

# need to use an over verion of apptainer (singularity replacement)
module load apptainer/1.1.8

# We need an assembly name, generally this is just the name of the organism"
Assembly_Name="SpeciesSexOptions"

# Hi-C Foward and reverse reads
Forward="RawData/*_R1.fastq.gz"
Reverse="RawData/*_R2.fastq.gz"

# HiFi Reads, ours are stored as BAM files
# If globbing give as much of the name as possible, the combined read file is named based on this
CCS='RawData/AsMuchInfoAsPossible*.bam'

# Current OTB needs a single HiFi input file, so if there's multiple we need to concatenate
if [[ $(ls $CCS | wc -w) -gt 1 ]]; then
  HiFi_fastq="${CCS%%\**}.combined.reads.fastq"
  # takes a long time so make sure this isn't already done
  if [[ ! -f "$HiFi_fastq" ]]; then
    echo "combining:"
    ls $CCS
    echo "into $HiFi_fastq"
    echo "this might take a minute"
    module load samtools
    for bam_file in $CCS; do
      samtools fastq $bam_file
    done > $HiFi_fastq
    module unload samtools
    echo "finished! combined reads are in $HiFi_fastq"
  fi
  echo "Using $HiFi_fastq as HiFi reads for OTB"
  CCS=${HiFi_fastq}
fi

# Singularity path assignment isn't working, pre-installed but need to link them where expected
ln -s /project/vpgru/software/NextflowSingularityLibrary/*.img work/singularity/

# HiFiASM can purge duplicate haplotigs to save memory. If trying trio-binning reduce to 1 or 0. See -l parameter in https://hifiasm.readthedocs.io/en/latest/parameter-reference.html#parameter-reference
PurgeDups=2

#Comment/Uncommment for busco
#Busco="--busco" #Busco will be run
Busco="" #Busco will not be run

#Comment/Uncommment for Yahs
Yahs="-y" #Yahs will be run
#Yahs="" #Yahs will not be run

#Comment/Uncomment for Polishing (only select one of)
#Polish_Type="" #No polishing
Polish_Type="simple" #Simple Polishing
#Polish_Type="dv" #Deep Variant Polishing
#Polish_Type="merfin" #merfin Polishing

#Comment/Uncomment for Type (only select one of)
#HiFi_Type="phasing"
HiFi_Type="default"
#HiFi_Type="trio"

#Comment/Uncomment for Runner (only select one of)
#Runner="slurm_usda"
Runner="slurm_usda_mem"

Threads="32"

Busco_Location="-l /project/ag100pest/software/OTB_test/busco"
Busco_DB="-p insecta_odb10"

if [[ -z "$BUSCO" ]]; then
  ./otb.sh -n ${Assembly_Name} -f "$( echo ${Forward})" -r "$(echo ${Reverse})" -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --purge-dups ${PurgeDups} ${Yahs} ${Busco} --polish_type ${Polish_Type} --runner ${Runner} -c -s
else
  ./otb.sh -n ${Assembly_Name} -f "$( echo ${Forward})" -r "$(echo ${Reverse})" -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --purge-dups ${PurgeDups} ${Yahs} ${Busco} ${Busco_Location} ${Busco_DB} --polish_type ${Polish_Type} --runner ${Runner} -c -s
fi
