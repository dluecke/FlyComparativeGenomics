# PlotsForAssemblyCuration.R takes Mummer output and makes standard set of plots
# to visulize assembly agreement, ID unplaced scaffolds, and check for duplication

# Intended to run via command line submission on Ceres

# USAGE: Rscript PlotsForAssemblyCuration.R INFILES_ALN.csv INFILES_SELFALN.csv OUTFILENAME

# INFILES_ALN.csv has alignment id, COORDS, FAI, and ref and qry axis labels, no header:
#     pctgs,asmM-vs-asmF.coords,asmM.fa.fai,asmF.fa.fai,male primary,female primary
#     hapF1,asmF-vs-Fh1.coords,asmF.fa.fai,Fh1.fa.fai,female primary,female hap1

# INFILES_SELFALN.csv is for self alignments (unplaced scaffolds onto chrs) 
#  can run from single assembly via self_align-ChrsVsUnplaced.sh
# alignment id (haps must match ids in INFILES_ALN.csv), 
#  COORDS, and FAI for self alignments, no headers:
#     Fasm,FChrs-vs-FUnplaced.coords,FChrs.fa.fai,FUnplaced.fa.fai
#     hapF1,Fh1Chrs-vs-Fh1Unplaced.coords,Fh1Chrs.fa.fai,Fh1Unplaced.fa.fai

# get command line input
args = commandArgs(trailingOnly = TRUE)
# check if arg count correct (before loading packages)
if (length(args) != 3) {
  stop("USAGE: Rscript PlotsForAssemblyCuration.R INFILES_ALN.csv INFILES_SELFALN.csv OUTFILENAME")
}


# path to VPGRU packages
.libPaths(new="/project/vpgru/software/R_packages")

# load custom functions
source("~/FlyComparativeGenomics/alignment_region_dotplot-FUNCTIONS.R")

# output filenames
OUT_TSV = paste0(args[3], ".tsv")
OUT_PDF = paste0(args[3], ".pdf")

# read main alignments files
INFILES_ALN = read.csv(args[1], row.names = 1, header = FALSE, stringsAsFactors = FALSE)

# split out labels dataframe
df.AxisLabels = INFILES_ALN[,c(4,5)]
colnames(df.AxisLabels) = c("RefLab", "QryLab")

# convert rest to list structure used by custom functions
l.INFILES_ALN = setNames(split(as.matrix(INFILES_ALN[,c(1:3)]), 
                               seq(nrow(INFILES_ALN))), 
                         rownames(INFILES_ALN))

# read self alignment files
INFILES_SELF = read.csv(args[2], row.names = 1, header = FALSE, stringsAsFactors = FALSE)

# convert to list structure
l.INFILES_SELF = setNames(split(as.matrix(INFILES_SELF), 
                                seq(nrow(INFILES_SELF))), 
                          rownames(INFILES_SELF))

# Process Data
# l.coords and l.homologs
l.ALN = make_l.l_coords.l_homologs(l.INFILES_ALN)
l.SELF = make_l.l_coords.l_homologs(l.INFILES_SELF)

# plots
l.PLOTS = plot_l.l_coords.l_homologs(l.ALN, df.AxisLabels)

# collect unplaced scaffolds by assembly
l.UNPLACED = make_l.unplaced(l.PLOTS)

# plot self align and call duplicates
l.SELF.PLOTS = plot_l.unplaced_self(l.UNPLACED, l.SELF)

# dataframe with unplaced and duplicate by assembly/scaffold
df.UNPLACED_DUPS = make_df.unplaced_dups(l.UNPLACED, l.SELF.PLOTS)

# output unplaced/duplicate df to table
write.table(df.UNPLACED_DUPS, file = OUT_TSV, sep = '\t',
            row.names = F, quote = F)

# plots to PDF
pdf(file = OUT_PDF)
l.PLOTS
l.SELF.PLOTS
dev.off()

