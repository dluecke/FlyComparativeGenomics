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
# To use this template first get OTB with `git clone git@github.com:molikd/otb.git` and copy template into new otb/ directory

module load nextflow/22.04.3

# need to use an older verion of apptainer (singularity replacement)
module load apptainer/1.1.8

# We need an assembly name, generally this is just the name of the organism
Assembly_Name="Species_Sex-Assembly_Type-N"

# READS
# HiFi Reads, ours are stored as BAM files, can use fastq
# If globbing give as much of the name as possible, the combined read file is named based on this
#CCS='/project/vpgru/PathToHiFi/AsMuchInfoAsPossible*.bam'
CCS='/project/vpgru/ComparativeGenomics/Musca_autumnalis/m64312e_22022.combined.reads.fastq'

# Hi-C Foward and reverse reads (won't use in Lite mode)
# 3/26/24 haven't yet checked if globbed Hi-C files are supported
#Forward="/project/vpgru/PathToHiC/*_R1.fastq.gz"
#Reverse="/project/vpgru/PathToHiC/*_R2.fastq.gz"
Forward='/project/vpgru/ComparativeGenomics/Musca_autumnalis/HiC/sa22218_S3HiC_R1.fastq.gz'
Reverse='/project/vpgru/ComparativeGenomics/Musca_autumnalis/HiC/sa22218_S3HiC_R2.fastq.gz'

# PIPELINE OPTIONS
# hifiasm phasing of diploid genome
#Comment/Uncomment for Type (only select one of)
HiFi_Type="default" # only use HiFi reads to generate primary/alternate contigs
#HiFi_Type="phasing" # attempt to generate haplotigs using Hi-C if available
#HiFi_Type="trio"    # use maternal/paternal reads to generate maternal/paternal chromosomes

# HiFiASM can purge duplicate haplotigs to save memory. If using trio-binning reduce to 1 or 0. See -l parameter in https://hifiasm.readthedocs.io/en/latest/parameter-reference.html#parameter-reference
PurgeDups=3

# Lite mode skips all Hi-C steps (phasing, polishing, scaffolding)
# if set to lite will automatically skip Hi-C read paths and next 2 options
Lite_Run=""
#Lite_Run="lite"

#Comment/Uncomment for Polishing (only select one of)
#Polish_Type="" #No polishing
Polish_Type="--polish_type simple" #Simple Polishing
#Polish_Type="--polish_type dv" #Deep Variant Polishing
#Polish_Type="--polish_type merfin" #merfin Polishing

#Comment/Uncommment for Yahs scaffolding with Hi-C reads
#Yahs="" #Yahs will not be run
Yahs="-y" #Yahs will be run

#Comment/Uncomment for Runner (only select one of)
#Runner="slurm_usda"
Runner="slurm_usda_mem"

Threads="32"

#Comment/Uncommment for busco quality check for benchmark single copy genes
#Busco="--busco" #Busco will be run
Busco="" #Busco will not be run

Busco_Location="-l /project/ag100pest/software/OTB_test/busco"
Busco_DB="-p insecta_odb10"

# RUNNING OTB (NO MORE TEMPLATE EDITS NEEDED BELOW HERE)
# Modifications by DL 
# Singularity path assignment isn't working, pre-installed but need to link them where expected
mkdir work
mkdir work/singularity
ln -s /project/vpgru/software/NextflowSingularityLibrary/*.img work/singularity/

# HiFi read processing pre-OTB
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
      # use samtools fastq to convert BAMs, then cat loop output to new file
      # --threads is number of additional, so 2 total
      samtools fastq --threads 1 $bam_file
    done > $HiFi_fastq
    module unload samtools
    echo "finished! combined reads are in $HiFi_fastq"
  fi
  echo "Using $HiFi_fastq as HiFi reads for OTB"
  CCS=${HiFi_fastq}
fi

# Call OTB with specified options
if [[ -z "$Lite_Run" ]]; then
# normal run (-z checks for empty string)
  if [[ -z "$Busco" ]]; then
  # no BUSCO
    ./otb.sh -n ${Assembly_Name} -f "$( echo ${Forward})" -r "$(echo ${Reverse})" -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --purge-dups ${PurgeDups} ${Yahs} ${Busco} ${Polish_Type} --runner ${Runner} -c -s
  else
  # include BUSCO
    ./otb.sh -n ${Assembly_Name} -f "$( echo ${Forward})" -r "$(echo ${Reverse})" -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --purge-dups ${PurgeDups} ${Yahs} ${Polish_Type} ${Busco} ${Busco_Location} ${Busco_DB} --runner ${Runner} -c -s
  fi
else
# lite run
  if [[ -z "$Busco" ]]; then
  # no BUSCO
    ./otb.sh --lite -n ${Assembly_Name} -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --purge-dups ${PurgeDups} --runner ${Runner} -c -s
  else
  # include BUSCO
    ./otb.sh --lite -n ${Assembly_Name} -in "$(echo ${CCS})" -m ${HiFi_Type} -t ${Threads} --purge-dups ${PurgeDups} ${Busco} ${Busco_Location} ${Busco_DB} --runner ${Runner} -c -s
  fi
fi