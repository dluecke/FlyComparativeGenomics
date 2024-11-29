#!/bin/bash
#SBATCH --job-name="split_fa"
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p short
#SBATCH -o "split_fa.stdout.%j.%N" 		# standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "split_fa.stderr.%j.%N" 		#optional but it prints our standard error

# using purge_dups split_fa tool to de-scaffold assembly - useful for blobtools BLAST step

module load purge_dups

split_fa $1 > $1.split.fa
