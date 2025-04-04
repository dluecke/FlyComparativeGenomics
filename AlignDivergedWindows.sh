#!/bin/bash

# AlignDivergedWindows.sh finds and aligns windows of divergence between neo-sex chromosomes
# Use to find sequence markers specific to heterogametic sex (XY males, ZW female) 
# Resort to this if alignment doesn't provide clear divergence windows

# Positional input:
# 1. name of target scaffold (for now needs to be same in all 3 assemblies)
# 2-4. nucmer alignment COORDS files:
#  hap1 vs hap2 of heterogametic sex
#  homogametic primary vs heterogametic hap1
#  homo primary vs het hap2
# 5-7. corresponding fastas for those assemblies:
#  het hap1
#  het hap2
#  homo primary

# Finds longest windows of elevated divergence between haps 
# and extracts homologous windows in het haps and homo primary 
# Aligns sequences using muscle5
# Outputs multifasta of windows and alignment MSA file

# Written to run on SciNet Ceres cluster 4/3/2025

module load samtools/1.17
module load muscle/5.1.0

# input target scaffold
FOCAL_SCAF=$1

# input coords files
COORDS_Het1vsHet2=$2
COORDS_HomPvsHet1=$3
COORDS_HomPvsHet2=$4

# input assemblies
ASM_Het1=$5
ASM_Het2=$6
ASM_HomP=$7

# constant (for now) variables 
WINDOW_SIZE=20000 # window size in bp
N_HITS=5 # number of regions to align
# can either filter by max percent identity then take longest hits
MAX_ID=95 # percent identity max for divergence windows
# or filter by minimum hit length then take lowest identity
MIN_HIT_LENGTH=7500 # shortest hit length for alignment seed

# function to get corresponding coordinates for hit boundaries in female
# adapted from answer here https://stackoverflow.com/questions/17853037/in-a-column-of-numbers-find-the-closest-value-to-some-target-value
# with tweaks and formatting by Copilot
# Function to find the closest value
find_closest_value() {
  local target=$1      # Target value
  local column_t=$2    # Column to search for closest value in file
  local column_o=$3    # Column with coordinate in aligned assembly
  local coords_file=$4 # Input coords file

  # Run AWK and print the output (capture when calling function)
  awk -v c=$column_t -v co=$column_o -v t=$target '
      NR == 1 {
          # Process the first line
          d = $c - t
          d = (d < 0 ? -d : d)   # Absolute difference
          v = $c                 # Closest value
          vo = $co      # Coordinate value from other assembly
          next                   # Skip to next line
      }
      {
          # Process subsequent lines
          m = $c - t
          m = (m < 0 ? -m : m)   # Absolute difference
          if (m < d) {
              d = m              # Update smallest difference
              v = $c             # Update closest value
              vo = $co  # Update column 1 value
          }
      }
      END {
          # Print the target, closest value and matching coordinate
          print t, v, vo
      }
  ' "$coords_file"
}

# Get rows for diverged windows between heterogametic haps, seeds for alignment region
# using maxid, for now opting for minlength approach
#awk -v maxid=$MAX_ID -v scaf=$FOCAL_SCAF '
#    $10 <= maxid && $12 == scaf && $13 == scaf
#    ' $COORDS_Het1vsHet2 | sort -nr -k7,7 | head -n${N_HITS} > \
#    $COORDS_Het1vsHet2.seedhits
# using minlength
awk -v minlength=$MIN_HIT_LENGTH -v scaf=$FOCAL_SCAF '
    $7 >= minlength && $12 == scaf && $13 == scaf
    ' $COORDS_Het1vsHet2 | sort -n -k10,10 | head -n${N_HITS} > \
    $COORDS_Het1vsHet2.seedhits
# possible that these hits are close enough that resulting alignments overlap

# for each seed hit pull region and run alignment
i=1
while read HIT; do

    # seed hit boundary coords and length for both het haps 
    S_Het1=$(echo $HIT | awk '{print $1}')
    E_Het1=$(echo $HIT | awk '{print $2}')
    L_Het1=$(echo $HIT | awk '{print $7}')
    S_Het2=$(echo $HIT | awk '{print $4}')
    E_Het2=$(echo $HIT | awk '{print $5}')
    L_Het2=$(echo $HIT | awk '{print $8}')

    # Closest corresponding coords in homogametic primary 
    # ideally would do this with the hap corresponding to neo-X(Z), but not typically known
    # find in both hap1 and hap2 alignments and go with the hit that has closest length
    # hap1 alignment
    # find_closest_value output lines
    ClosestOut_h1S1=$(find_closest_value $S_Het1 4 1 <(awk -v scaf=$FOCAL_SCAF '
                                                    $12 == scaf && $13 == scaf
                                                    ' $COORDS_HomPvsHet1))
    ClosestOut_h1E1=$(find_closest_value $E_Het1 5 2 <(awk -v scaf=$FOCAL_SCAF '
                                                    $12 == scaf && $13 == scaf
                                                    ' $COORDS_HomPvsHet1))
    # adjust hit to difference between hap1 target and closest match
    Sh1_HomP=$(( $(echo $ClosestOut_h1S1 | awk {'print $3'}) \
               + $(echo $ClosestOut_h1S1 | awk {'print $1'}) \
               - $(echo $ClosestOut_h1S1 | awk {'print $2'}) ))
    Eh1_HomP=$(( $(echo $ClosestOut_h1E1 | awk {'print $3'}) \
               + $(echo $ClosestOut_h1E1 | awk {'print $1'}) \
               - $(echo $ClosestOut_h1E1 | awk {'print $2'}) ))
    Lh1_HomP=$(( Eh1_HomP - Sh1_HomP ))
    
    # hap2 alignment
    # find_closest_value output lines
    ClosestOut_h2S1=$(find_closest_value $S_Het2 4 1 <(awk -v scaf=$FOCAL_SCAF '
                                                    $12 == scaf && $13 == scaf
                                                    ' $COORDS_HomPvsHet1))
    ClosestOut_h2E1=$(find_closest_value $E_Het2 5 2 <(awk -v scaf=$FOCAL_SCAF '
                                                    $12 == scaf && $13 == scaf
                                                    ' $COORDS_HomPvsHet1))
    # adjust hit to difference between hap1 target and closest match
    Sh2_HomP=$(( $(echo $ClosestOut_h2S1 | awk {'print $3'}) \
               + $(echo $ClosestOut_h2S1 | awk {'print $1'}) \
               - $(echo $ClosestOut_h2S1 | awk {'print $2'}) ))
    Eh2_HomP=$(( $(echo $ClosestOut_h2E1 | awk {'print $3'}) \
               + $(echo $ClosestOut_h2E1 | awk {'print $1'}) \
               - $(echo $ClosestOut_h2E1 | awk {'print $2'}) ))
    Lh2_HomP=$(( Eh2_HomP - Sh2_HomP ))

    # compare hom implied hit lengths from both hap assemblies
    # difference and absolute value (drop any leading '-')
    dLh1=$((Lh1_HomP - L_Het1))
    dLh1=${dLh1#-}
    dLh2=$((Lh2_HomP - L_Het2))
    dLh2=${dLh2#-}
    if [[ $dLh1 -gt $dLh2 ]]; then
        S_HomP=$Sh2_HomP
        E_HomP=$Eh2_HomP
        L_HomP=$Lh2_HomP
    else
        S_HomP=$Sh1_HomP
        E_HomP=$Eh1_HomP
        L_HomP=$Lh1_HomP
   fi

    # lengths to pad coordinate ends to get to window length
    PAD_Het1=$(( (WINDOW_SIZE - L_Het1)/2 ))
    PAD_Het2=$(( (WINDOW_SIZE - L_Het2)/2 ))
    PAD_HomP=$(( (WINDOW_SIZE - L_HomP)/2 ))

    # extract region from het1
    S1=$((S_Het1 - PAD_Het1))
    E1=$((E_Het1 + PAD_Het1))
    samtools faidx -o "${ASM_Het1%.*}-${FOCAL_SCAF}-${S1}to${E1}.fa" \
        -r <(echo "${FOCAL_SCAF}:$S1-$E1") $ASM_Het1

    # extract region from het2
    S2=$((S_Het2 - PAD_Het2))
    E2=$((E_Het2 + PAD_Het2))
    samtools faidx -o "${ASM_Het2%.*}-${FOCAL_SCAF}-${S2}to${E2}.fa" \
        -r <(echo "${FOCAL_SCAF}:$S2-$E2") $ASM_Het2

    # extract region from homP
    SP=$((S_HomP - PAD_HomP))
    EP=$((E_HomP + PAD_HomP))
    samtools faidx -o "${ASM_HomP%.*}-${FOCAL_SCAF}-${SP}to${EP}.fa" \
        -r <(echo "${FOCAL_SCAF}:$SP-$EP") $ASM_HomP
    
    # add assembly info to sequence names and concatenate into single fasta 
    sed "s/scaffold/${ASM_Het1%.*}-scaffold/" ${ASM_Het1%.*}-${FOCAL_SCAF}-${S1}to${E1}.fa > DivergedWindow$i.fa
    sed "s/scaffold/${ASM_Het2%.*}-scaffold/" ${ASM_Het2%.*}-${FOCAL_SCAF}-${S2}to${E2}.fa >> DivergedWindow$i.fa
    sed "s/scaffold/${ASM_HomP%.*}-scaffold/" ${ASM_HomP%.*}-${FOCAL_SCAF}-${SP}to${EP}.fa >> DivergedWindow$i.fa

    # submit alignment job to SLURM
    sbatch ~/FlyComparativeGenomics/muscle.slurm DivergedWindow$i.fa

    ((i+=1))
done < $COORDS_Het1vsHet2.seedhits
