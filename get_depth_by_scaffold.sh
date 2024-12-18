#!/bin/bash

# get depth by scaffold
# takes samtools depth output TSV with 3 columns: seqid position depth
# make sure to use -a option for samtools depth to report sites with depth=0
# outputs depth_by_scaffold.tsv with 3 columns: seqid n_positions sum(depth)/n_positions
DEPTHFILE=$1

if [[ $DEPTHFILE =~ "depth" ]]; then
    OUTFILE=$(echo $DEPTHFILE | sed "s/depth/depth_by_scaffold/")
else
    OUTFILE="$DEPTHFILE.depth_by_scaffold.tsv"
fi

awk -v OFS='\t' '{
   sum[$1] += $3; count[$1]++
} END { 
    n=asorti(count, indices); 
    for (i=1; i<=n; i++) {
        print indices[i], count[indices[i]], sum[indices[i]]/count[indices[i]]
    } 
}' $DEPTHFILE > $OUTFILE

