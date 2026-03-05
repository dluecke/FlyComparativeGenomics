#!/bin/bash
#SBATCH -A vpgru
#SBATCH --job-name="wdt"
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -p ceres
#SBATCH -t 12:00:00
#SBATCH -o "wdt.%j.%N.stdout"               # standard out %j adds job number to outputfile name and %N adds the node
#SBATCH -e "wdt.%j.%N.stderr"               #optional but it prints our standard error

# touch all files in working directory or subdirectories

find . -type f -exec touch {} +
