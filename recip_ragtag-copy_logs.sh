#!/bin/bash

# recip_ragtag-copy_logs.sh finds and copies log files from reciprocal ragtag steps
# makes mirror directory structure and puts logs in same location to avoid same-name confusion
# written to be called at end of Snakemake pipeline, but can be used standalone
# search strings for each file set should be mutually exclusive, otherwise may get recursive copying

# mirror directory stucture in log_copies/
# printf format removes leading "./"
for d in $(find . -type d -printf '%P\n'); do
    mkdir -p log_copies/$d
done

# read mapping and assembly statistics
for s in $(find . -name "*.stats"); do
    cp $s log_copies/$s
done

# output from slurm jobs
for r in $(find . -name "rrt.*"); do
    cp $r log_copies/$r
done

# AGPs from ragtag correct and scaffold
for a in $(find . -name "*.agp"); do 
    cp $a log_copies/$a
done

# RagTagShield output: dfJumps and GFFs used in ragtag correct
for c in $(find . -name "*.coords-*"); do
    cp $c log_copies/$c
done

# internal ragtag .log files
for l in $(find . -name "*.log"); do
    cp $l log_copies/$l
done

# info .txt files, including pct_chrs.txt
for t in $(find . -name "*.txt"); do
    cp $t log_copies/$t
done
