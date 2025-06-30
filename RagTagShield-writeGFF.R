# RagTagShield-writeGFF.R writes GFFs for regions to shield from RagTag correct

# Part of modifications to RagTag process for improving two similar-quality draft assemblies
# of same species (eg female and male) by placing small scaffolds into gaps based on 
# alignment position in other assembly
# RagTag correct may over-correct in cases where "query" assembly has contiguous sequence 
# missing from "reference" (even when using --intra, although that option does help)
# This approach aims to identify and shield these regions by writing GFF to pass 
# to RagTag correct --gff option
# Can then use RagTag in both directions to fill assembly gaps with unplaced scaffolds
# without creating new gaps in areas current "query" is better than "reference"

# Takes .coords file from nucmer whole genome alignment
# Uses tidyverse and slider packages (installed in Ceres VPGRU project space)


# get command line input
args = commandArgs(trailingOnly = TRUE)
# check if arg count correct (before loading packages)
if (length(args) != 1) {
  stop("USAGE: Rscript RagTagShield-writeGFF.R ALIGNMENT.coords")
}


# path to VPGRU packages
.libPaths(new="/project/vpgru/software/R_packages")

library(tidyverse)
library(slider)

OUT_GFF.QRY = paste0(args[1], "-BetterInQry.gff", sep = "")
OUT_GFF.REF = paste0(args[1], "-BetterInRef.gff", sep = "")


# function to read coords file and produce df.jumps
# includes scaled difference between alignment gap sizes in two assemblies
make_df.jumps_fromCoords = function(IN_COORDS_FILE){
  # read coords format, skip first 5 header lines
  df.coords <- read.table(IN_COORDS_FILE, stringsAsFactors = F, skip = 5)
  # drop formatting separating columns
  df.coords <- df.coords[,-c(3,6,9,11,14)]
  # rename to work with downstream subsetting/plot functions
  colnames(df.coords) <- c("R1", "R2", "Q1", "Q2", 
                           "ref_length", "qry_length", 
                           "pctID",
                           "ref_SeqName", "qry_SeqName")
  # assign orientation, in coords format reverse matches shown in query start/stop reversal
  df.coords$orientation <- apply(df.coords, 1, FUN = function(x){
    if( x[3] <= x[4] ){return("forward")}
    if( x[3] >= x[4] ){return("reverse")}
  })
  
  # most common alignment pairs (allow IDing homologs without requiring same name)
  df.seq_pairs = df.coords %>% 
    group_by(qry_SeqName) %>% 
    count(ref_SeqName, qry_SeqName) %>% 
    slice(which.max(n)) %>% 
    group_by(ref_SeqName) %>%
    slice(which.max(n)) %>%
    as.data.frame()
  row.names(df.seq_pairs) = df.seq_pairs$ref_SeqName
  
  df.jumps = df.coords %>%
    filter(qry_SeqName == df.seq_pairs[ref_SeqName,]$qry_SeqName,
           orientation == "forward",
           ref_length > median(ref_length)) %>% 
    group_by(ref_SeqName) %>% filter(n() > 100) %>% 
    mutate(coord_diff = R1-Q1) %>%
    # filter points more than 1 sd off the diagonal (for whole sequence)
    filter( abs(coord_diff) < median(abs(coord_diff)) + 1*sd(abs(coord_diff)) ) %>%
    # sliding windows for more precise off-diagonal filtering
    mutate(
      # weighted mean of coordinate difference by region length
      wtmn.cdiff.win = slide_dbl(
        .x = coord_diff, 
        .w = mean(c(ref_length, qry_length)),
        weighted.mean,
        .before = 50,
        .after = 49,
        .complete = F),
      # SD of coordinate differences
      sd.cdiff.win = slide_dbl(
        .x = coord_diff,
        sd,
        .before = 50,
        .after = 49,
        .complete = F)
    ) %>%
    # keep if within 1 SD of window weighted mean
    filter(coord_diff > wtmn.cdiff.win - sd.cdiff.win,
           coord_diff < wtmn.cdiff.win + sd.cdiff.win) %>%
    # ensure coordinates always increasing to remove inverted alignments
    filter(sapply(row_number(), function(i){
      all(R2[i] < R1[(i+1):n()]) &&
        all(Q2[i] < Q1[(i+1):n()])
    })) %>%
    # values for jumps between aligned regions
    mutate(Rnext = lead(R1),
           Qnext = lead(Q1),
           Rjump = Rnext - R2, 
           Qjump = Qnext - Q2,
           # difference in alignment gaps between two assemblies, scaled by average gap size
           d_jump = (Rjump - Qjump) / ((Rjump + Qjump)/2) ) %>%
    arrange(ref_SeqName)
  
  return(df.jumps)
}

df.JUMPS = make_df.jumps_fromCoords(args[1])

# query-better outliers in gap size difference 
df.JUMPS_QryBetter = df.JUMPS %>%
  filter(d_jump < median(d_jump, na.rm = T) - 2*sd(d_jump, na.rm = T))

# gff formatted table for "QryBetter" regions to be shielded
df.GFF_QryBetter = data.frame(
  SEQ = df.JUMPS_QryBetter$qry_SeqName,
  SOURCE = "RagTagShield",
  FEATURE = "BetterInQry",
  # coordinates, buffer by 100bp in case of alignment ambiguity
  START = df.JUMPS_QryBetter$Q2 - 100,
  END = df.JUMPS_QryBetter$Qnext + 100,
  SCORE = ".",
  STRAND = ".",
  FRAME = ".")

# write output GFF file
write.table(df.GFF_QryBetter, file = OUT_GFF.QRY, sep = '\t',
            col.names = F, row.names = F, quote = F)

# same for better in ref, allows single pre-RagTag alignment to inform both directions
df.JUMPS_RefBetter = df.JUMPS %>%
  filter(d_jump > median(d_jump, na.rm = T) + 2*sd(d_jump, na.rm = T))

# gff formatted table for "QryBetter" regions to be shielded
df.GFF_RefBetter = data.frame(
  SEQ = df.JUMPS_RefBetter$ref_SeqName,
  SOURCE = "RagTagShield",
  FEATURE = "BetterInRef",
  # coordinates, buffer by 100bp in case of alignment ambiguity
  START = df.JUMPS_RefBetter$R2 - 100,
  END = df.JUMPS_RefBetter$Rnext + 100,
  SCORE = ".",
  STRAND = ".",
  FRAME = ".")

# write output GFF file
write.table(df.GFF_RefBetter, file = OUT_GFF.REF, sep = '\t',
            col.names = F, row.names = F, quote = F)
