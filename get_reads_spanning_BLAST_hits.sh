#!/bin/bash

# get_reads_spanning_BLAST_hits.sh identifies BLAST hit clusters
# and gets reads mapping to that region (plus boundary spanning reads separately)
# Intended to find HiFi reads mapping to NuMt boundaries.
# takes positional arguments:
#   BLAST output in outfmt 6
#   BAM file for reads mapped to sequences matching BLAST target
#  Optional (all or none, maybe use case getopts later):
#   Minimum bp distance between clusters, DEFAULT 50000
#   bit score threshold to include cluster, DEFAULT 5000
#   boundary span bp distance for get_spanning_reads.sh, DEFAULT 5000
# outputs for each cluster:
#   SAM file with all hits 
#   SAM files for all boundary spanning hits via get_spanning_reads.sh
#   union list of reads from both boundaries
# requires:
#  samtools, written with version 1.17

# path for get_spanning_reads.sh
GIT_PATH=~/FlyComparativeGenomics

# required positional arguments
BLAST_HITS=$1
BAM_IN=$2
# optional positions, set defaults
[[ -n $3 ]] && MIN_GAP=$3 || MIN_GAP=50000
[[ -n $4 ]] && BIT_THRESH=$4 || BIT_THRESH=5000
[[ -n $5 ]] && BOUNDARY_SPAN=$5 || BOUNDARY_SPAN=5000

# need to keep each target sequence separate, loop through each SEQ with hit above bit score threshold
while read SEQ; do

    # array with start and stop of each possible cluster regardless of hit quality
    ALL_CLUSTER_EDGES=($(\
        awk -v seq=$SEQ '
            $2 == seq
        ' $BLAST_HITS | \
        cut -f9,10 | tr '\t' '\n' | sort -n | \
        awk -v OFS='\n' -v gap=$MIN_GAP '
            NR == 1 {
                print $1; last=$1
            };
            NR > 1 {
                if($1-last > gap){
                    print last, $1
                }; 
                last = $1}
            END {print $1} \
        ' \
    ))

    # array with starts of seed hits (ok if in same cluster)
    SEED_STARTS=($(\
        awk -v seq=$SEQ -v b=$BIT_THRESH '
            $2 == seq && $12 > b { 
                print $9 
            }
        ' $BLAST_HITS | sort -n \
    ))
    
    # arrays for boundaries of clusters in SEQ
    CLUSTER_BEGS=()
    CLUSTER_ENDS=()
    # loop through seed hits
    for s in ${SEED_STARTS[@]}; do
        # check every potential cluster edge
        for i in $(seq 0 $((${#ALL_CLUSTER_EDGES[@]}-1))); do
            # loop until seed coordinate greater than or equal to edge position
            if [[ ! ${ALL_CLUSTER_EDGES[$i]} -lt $s ]]; then
                # break to next seed if the same cluster as previous so no duplicate clusters
                [[ $i == $LASTi ]] && break || LASTi=$i
                # i will usually be even unless seed equals cluster end boundary 
                if ((i%2 == 0)); then # even index, typical case
                    CLUSTER_BEGS+=(${ALL_CLUSTER_EDGES[$i]})
                    CLUSTER_ENDS+=(${ALL_CLUSTER_EDGES[$((i+1))]})
                else # odd index is boundary end, need to shift values
                    CLUSTER_BEGS+=(${ALL_CLUSTER_EDGES[$((i-1))]})
                    CLUSTER_ENDS+=(${ALL_CLUSTER_EDGES[$i]})
                fi
                break
            fi
        done
    done

    # loop through all clusters/regions in SEQ
    for j in $(seq 0 $((${#CLUSTER_BEGS[@]}-1))); do

        # cluster edges
        REGION_BEG=${CLUSTER_BEGS[$j]}
        REGION_END=${CLUSTER_ENDS[$j]}

        # outfile names
        # all reads mapped to cluster region
        REGION_SAM=${BLAST_HITS%.*}-bit${BIT_THRESH}-${MIN_GAP}gap-${BAM_IN%.*}-${SEQ}_${REGION_BEG}to${REGION_END}.sam
        # reads spanning boundaries, default names from get_spanning_reads.sh
        BOUNDARYBEG_SAM=${BAM_IN%.*}-${SEQ}_${REGION_BEG}span${BOUNDARY_SPAN}bp.sam
        BOUNDARYEND_SAM=${BAM_IN%.*}-${SEQ}_${REGION_END}span${BOUNDARY_SPAN}bp.sam
        # list of reads spanning either boundary (input to seqtk subseq)
        BOUNDARY_READS=${BLAST_HITS%.*}_bit${BIT_THRESH}-${MIN_GAP}gap-${BAM_IN%.*}-${SEQ}_${REGION_BEG}or${REGION_END}span${BOUNDARY_SPAN}bp.lst

        # samtools view to extract reads mapped to cluster region
        samtools view $BAM_IN ${SEQ}:${REGION_BEG}-${REGION_END} > $REGION_SAM

        # boundary spanning reads
        $GIT_PATH/get_spanning_reads.sh $BAM_IN $SEQ $REGION_BEG $BOUNDARY_SPAN
        $GIT_PATH/get_spanning_reads.sh $BAM_IN $SEQ $REGION_END $BOUNDARY_SPAN

        # list reads from both boundaries, non-redundant
        cat $BOUNDARYBEG_SAM $BOUNDARYEND_SAM | cut -f1 | sort -u > $BOUNDARY_READS

    done

done < <(awk -v b=$BIT_THRESH '$12 > b {print $2}' $BLAST_HITS | sort -u)
