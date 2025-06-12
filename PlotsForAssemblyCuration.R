# PlotsForAssemblyCuration.R takes Mummer output and makes standard set of plots
# to visulize assembly agreement, ID unplaced scaffolds, and check for duplication

# Intended to run via command line submission on Ceres

# USAGE: Rscript PlotsForAssemblyCuration.R INFILES_ALN.csv INFILES_SELFALN.csv 

# INFILES_ALN.csv has alignment id, COORDS, FAI, and ref and qry axis labels, no header:
#     pctgs,asmM-vs-asmF.coords,asmM.fa.fai,asmF.fa.fai,male primary,female primary
#     hapF1,asmF-vs-Fh1.coords,asmF.fa.fai,Fh1.fa.fai,female primary,female hap1

# INFILES_SELFALN.csv is for self alignments (unplaced scaffolds onto chrs) 
#  can run from single assembly via self_align-ChrsVsUnplaced.sh
# alignment id (haps must match ids in INFILES_ALN.csv), 
#  COORDS, and FAI for self alignments, no headers:
#     Fasm,FChrs-vs-FUnplaced.coords,FChrs.fa.fai,FUnplaced.fa.fai
#     hapF1,Fh1Chrs-vs-Fh1Unplaced.coords,Fh1Chrs.fa.fai,Fh1Unplaced.fa.fai


# path to VPGRU packages
.libPaths(new="/project/vpgru/software/R_packages")

# load custom functions
source("~/FlyComparativeGenomics/alignment_region_dotplot-FUNCTIONS.R")

# get command line input
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("USAGE: Rscript PlotsForAssemblyCuration.R INFILES_ALN.csv INFILES_SELFALN.csv")
}

# read main alignments files
INFILES_ALN = read.csv(args[1], row.names = 1, header = FALSE, stringsAsFactors = FALSE)

# split out labels dataframe
df.AxisLabels = INFILES_ALN[,c(4,5)]
colnames(df.AxisLabels) = c("RefLab", "QryLab")

# convert rest to list structure used by custom functions
l.INFILES_ALN = setNames(split(INFILES_ALN[,c(1:3)], seq(nrows(INFILES_ALN))), 
                         rownames(INFILES_ALN))

# read self alignment files
INFILES_SELF = read.csv(args[2], row.names = 1, header = FALSE, stringsAsFactors = FALSE)

# convert to list structure
l.INFILES_SELF = setNames(split(INFILES_SELF, seq(nrows(INFILES_SELF))), 
                          rownames(INFILES_SELF))


write.table(data.frame(l.INFILES_ALN), file = "test-ALN.tsv", sep = '\t', row.names = F)
write.table(df.AxisLabels, file = "test-labels.tsv", sep = '\t', row.names = F)
write.table(data.frame(l.INFILES_SELF), file = "test-SELF.tsv", sep = '\t', row.names = F)
