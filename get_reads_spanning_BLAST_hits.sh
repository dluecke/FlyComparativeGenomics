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
#   report file for cluster with edge coords, hits, boundary spanning reads
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

    # want an easy-to-read tag for cluster id
    # alphanumeric-only version of SEQ name
    SEQ_ALNUM=$(echo $SEQ | tr -d -c '[:alnum]')
    # using letters since sequences and coordinates already numeric
    if [[ ! ${#CLUSTER_BEGS[@]} -gt 26 ]]; then
        LETTER_ID=( {A..Z} )
    else # if more than 26 clusters use two-letter (AA AB .. AZ BA BB .. etc) which allows 676
        LETTER_ID=( {A..Z}{A..Z} )
    fi

    # loop through all clusters/regions with high scoring hits in SEQ
    for j in $(seq 0 $((${#CLUSTER_BEGS[@]}-1))); do

        # cluster id tag
        CLUSTERTAG=${SEQ_ALNUM}_${LETTER_ID[$j]}

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
        # report filename
        REGION_REPORT=HitClusterInfo-${CLUSTERTAG}.txt

        # samtools view to extract reads mapped to cluster region
        samtools view $BAM_IN ${SEQ}:${REGION_BEG}-${REGION_END} > $REGION_SAM

        # boundary spanning reads
        $GIT_PATH/get_spanning_reads.sh $BAM_IN $SEQ $REGION_BEG $BOUNDARY_SPAN
        $GIT_PATH/get_spanning_reads.sh $BAM_IN $SEQ $REGION_END $BOUNDARY_SPAN

        # list reads from both boundaries, non-redundant
        cat $BOUNDARYBEG_SAM $BOUNDARYEND_SAM | cut -f1 | sort -u > $BOUNDARY_READS

        # cluster reporting
        N_CLUSTER_HITS=$(awk -v s=$SEQ -v b=$REGION_BEG -v e=$REGION_END '$2 == s && $9 >= b && $9 <= e' $BLAST_HITS | wc -l)
        {
        echo -e "Details for cluster ${LETTER_ID[$j]} of BLAST hits on $SEQ\n"
        echo -e "\nClustering Parameters:"
        echo -e "cluster best hit score threshold:\t$BIT_THRESH"
        echo -e "minimum distance between clusters:\t$MIN_GAP"
        echo -e "span distance for boundary reads:\t$BOUNDARY_SPAN"
        echo -e "\nBLAST Results"
        echo -e "input hits file:\t$BLAST_HITS"
        echo -e "sequence name:\t$SEQ"
        echo -e "number of hits:\t$N_CLUSTER_HITS"
        echo -e "region length:\t$((REGION_END - REGION_BEG + 1))"
        echo -e "region coordinates:\t$REGION_BEG - $REGION_END"
        echo -e "\nBLAST hits in cluster:"
        awk -v s=$SEQ -v b=$REGION_BEG -v e=$REGION_END '
            $2 == s && $9 >= b && $9 <= e
        ' $BLAST_HITS | sort -n -k9,10
        echo -e "\nAligned Reads"
        echo -e "input alignment file:\t$BAM_IN"
        echo -e "number of reads in region:\t$(wc -l $REGION_SAM | awk '{print $1}')"
        echo -e "all alignments in region:\t$REGION_SAM"
        echo -e "number of reads spanning region start:\t$(wc -l $BOUNDARYBEG_SAM | awk '{print $1}')"
        echo -e "number of reads spanning region end:\t$(wc -l $BOUNDARYEND_SAM | awk '{print $1}')"
        echo -e "alignments spanning region start:\t$BOUNDARYBEG_SAM"
        echo -e "alignments spanning region end:\t$BOUNDARYEND_SAM"
        echo -e "list of read IDs spanning either boundary:\t$BOUNDARY_READS"
        echo -e "IDs of reads spanning region start:"
        printf "\t%s\n" $(cut -f1 $BOUNDARYBEG_SAM)
        echo -e "IDs of reads spanning region end:"
        printf "\t%s\n" $(cut -f1 $BOUNDARYEND_SAM)
        } > $REGION_REPORT

    done

done < <(awk -v b=$BIT_THRESH '$12 > b {print $2}' $BLAST_HITS | sort -u)
