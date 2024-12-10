#!/bin/bash

# get depth by windows
DEPTHFILE=$1

[[ -n $2 ]] && WIN_KB=$2 || WIN_KB=100

WINDOW=$((WIN_KB*1000))

OUTFILE=$(echo $DEPTHFILE | sed "s/depth/depth_by_windows${WIN_KB}kb/")

awk -v win=$WINDOW -v OFS='\t' '{w=int($2/win)} {print $0, $1"-"w}' $DEPTHFILE | \
    awk -v OFS='\t' '{
        sum[$4] += $3; count[$4]++
    } END { 
        n=asorti(count, indices); 
        for (i=1; i<=n; i++) {
            print indices[i], count[indices[i]], sum[indices[i]]/count[indices[i]]
        } 
    }' | \
    awk -v OFS='\t' '{
        split($1, arr, "-")
    } {
        print length(arr[1]), arr[1], arr[2], $0
    }' | \
    sort -k1,1n -k2,2 -k3,3n | cut -f2-6 > $OUTFILE
