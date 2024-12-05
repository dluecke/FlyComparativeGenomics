#!/bin/bash

# get depth by windows
WINDOW=100000
DEPTHFILE=$1
OUTFILE=$(echo $DEPTHFILE | sed 's/depth/depth_by_windows/')
awk -v win=$WINDOW -v OFS='\t' '{w=int($2/win)} {print $0, $1"-"w}' $DEPTHFILE | \
    awk -v OFS='\t' '{
        sum[$4] += $3; count[$4]++
    } END { 
        n=asorti(count, indices); 
        for (i=1; i<=n; i++) {
            print indices[i], count[indices[i]], sum[indices[i]]/count[indices[i]]
        } 
    }' > temp_depth_windows.tsv

awk -v OFS='\t' '{
        split($1, arr, "-")
    } {
        print length(arr[1]), arr[1], arr[2], $0
    }' temp_depth_windows.tsv | \
    sort -k1,1n -k2,2 -k3,3n | cut -f2-6 > $OUTFILE

