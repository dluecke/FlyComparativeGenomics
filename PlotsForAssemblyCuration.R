# PlotsForAssemblyCuration.R takes Mummer output and makes standard set of plots
# to visulize assembly agreement, ID unplaced scaffolds, and check for duplication

# Intended to run via command line submission on Ceres

# USAGE: Rscript PlotsForAssemblyCuration.R INFILES_ALN.csv INFILES_SELFALN.csv OUTFILENAME

# INFILES_ALN.csv has alignment id, refASM, qryASM, COORDS, FAI, ref,qry axis labels, no header:
#     pctgs,asmM,asmF,asmM-vs-asmF.coords,asmM.fa.fai,asmF.fa.fai,male primary,female primary
#     hapF1,asmF,hapF1,asmF-vs-Fh1.coords,asmF.fa.fai,Fh1.fa.fai,female primary,female hap1

# INFILES_SELFALN.csv is for self alignments (unplaced scaffolds onto chrs) 
#  can run from single assembly via self_align-ChrsVsUnplaced.sh
# alignment id (all must match refASM/qryASM and haps must match alignment ids in INFILES_ALN.csv), 
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

# alignment to assembly map
df.ALNtoASM = INFILES_ALN[,c(1,2)]
colnames(df.ALNtoASM) = c("RefAsm", "QryAsm")

# split out labels dataframe
df.AxisLabels = INFILES_ALN[,c(6,7)]
colnames(df.AxisLabels) = c("RefLab", "QryLab")

# convert rest to list structure used by custom functions
l.INFILES_ALN = setNames(split(as.matrix(INFILES_ALN[,c(3:5)]), 
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

# all chromosomes plot for each alignment
lapply(l.PLOTS, function(ALN) ALN$p.Chr)

# by-assembly by-chromosome alignments, unplaced scaffolds, self align (duplicates)
lapply(names(l.UNPLACED), 
       # per assembly
       function(ASMi){
         # for haps ASM and ALN have same name, can access l.PLOT elements directly
         if(ASMi %in% names(l.PLOTS)){
           l.l.p = lapply(names(l.UNPLACED[[ASMi]]),
                  # per chromosome
                  function(CHR){
                    # chromosome hap onto primary
                    l.p.C = l.PLOTS[[ASMi]]$l.p.Chr[[CHR]]
                    # hap unplaced scaffolds onto primary
                    l.p.CP = l.PLOTS[[ASMi]]$l.p.ChrPlacement[[CHR]]$p.IntoQry
                    # hap self align scaffolds onto chromosome, shows duplicates
                    l.s.p = l.SELF.PLOTS[[ASMi]][[CHR]]$plot
                    R = list(l.p.C, l.p.CP, l.s.p)
                    return(R)
                  })
         # primary ASM and ALN names different, need to map to alignments
         } else {
           # alignments where ASM is reference
           v.ALN_AsRef = row.names(df.ALNtoASM)[df.ALNtoASM$RefAsm == ASMi]
           # alignments where ASM is query
           v.ALN_AsQry = row.names(df.ALNtoASM)[df.ALNtoASM$QryAsm == ASMi]
           # non-hap (primaries) alignment (name not in asm list), 
           #  should only be 1 but take [1] to make sure
           PRI_ALN = c(v.ALN_AsRef, v.ALN_AsQry)[
             !c(v.ALN_AsRef, v.ALN_AsQry) %in% names(l.UNPLACED)
             ][1]
           l.l.p = lapply(names(l.UNPLACED[[ASMi]]),
                  function(CHR){
                    # chromosome primary onto primary
                    l.p.C = l.PLOTS[[PRI_ALN]]$l.p.Chr[[CHR]]
                    # do AsQry first, either empty or has primaries aln
                    l.p.AQ = list()
                    if(length(v.ALN_AsQry) > 0){
                      l.p.AQ = lapply(v.ALN_AsQry, function(ALN_Q){
                        # unplaced primary scaffolds onto other primary 
                        return(l.PLOTS[[ALN_Q]]$l.p.ChrPlacement[[CHR]]$p.IntoQry)
                      })
                    }
                    l.p.AR = list()
                    if(length(v.ALN_AsRef) > 0){
                      l.p.AR = lapply(v.ALN_AsRef, function(ALN_R){
                        # unplaced primary scaffolds onto hap
                        return(l.PLOTS[[ALN_R]]$l.p.ChrPlacement[[CHR]]$p.IntoRef)
                      })
                    }
                    # pri self align scaffolds onto chromosome, shows duplicates
                    l.s.p = l.SELF.PLOTS[[ASMi]][[CHR]]$plot
                    R = list(l.p.C, l.p.AQ, l.p.AR, l.s.p)
                    return(R)
                  })
         }
         return(l.l.p)
       }
)


dev.off()

