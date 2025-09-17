# SexSpecificScaffolds-FUNCTIONS.R
# functions for testing scaffold depth and variant data for sex chromosome signals

library(tidyverse)

# INPUT FILE FORMATS
###
# n_masked files - 1 of each CSV per assembly
#
# ASM.nmasked.csv: 
#   scaffold,length,masked
# row per scaffold with scaf name, total length, and n hardmasked bases (Ns)
# also need for unmasked assemblies to get count of Ns from scaffold gaps etc
# written by hardmask.slurm (same time as converting softmask to hardmask)
# if run on unmasked will just do the N counts
#
# ASM.windowsNkb_nmasked.csv:
#   window,length,masked
# same format as above but row per window
# written by get_nmasked_by_window.sh on hardmasked or unmasked assembly
#
###
# depth files - 1 of each TSV per seq library per assembly mapping replicate
#
# ASM-SAMPLE.depth_by_scaffold.tsv:
#   scaffold  n_sites avg_depth
# row per scaffold with scaf name, # sites with depth reported, and average depth
# written by get_depth_by_scaffold.sh
# 
# ASM-SAMPLE.depth_by_windowsXkb.tsv:
#   scaffold  window_number window_ID n_sites avg_depth
# written from depth.tsv file by get_depth_by_window.sh, called by window-depth.slurm
#
###
# variant files - 1 of each per seq library per assembly
#
# ASM-SAMPLE.bam.nseg_by_scaffold.csv:
#   scaffold,n_seg_sites
# generated from VCF file of BAM read alignments during
# bcfcall.slurm job generating VCF
# number of segregating sites per scaffold (dumb count, no fr(var) considered)
#
# ASM-SAMPLE.bam.vcf.nvar_by_windowsXkb.tsv:
#   scaffold  window_number window_ID n_seg_sites
    

# normalize_depth() function to add columns for by-unmasked-bp and by-library 
#  normalized depths per region (scaffold or window)
# Accounting for Masked/Unmasked BP:
#  depending on samtools view options (-a, -aa) the depth.tsv file may include:
#   - no 0s regardless of site masking
#   - 0s if site is unmasked but without mapped reads (true 0)
#   - 0s if site is masked (artifact 0)
#  the raw average doesn't account for this difference, but the 
#  unmasked depth average should include "true" but not "artifact" 0s
#    (here "masked" applies to all N bases, including eg scaffold gaps in unmasked assembly)
#  Can convert the raw_avg from any of these scenarios into true unmasked average depth
#  with number of sites used in the average and number of unmasked bases in sequence:
#    true_unmask_avg = raw_avg * (n_sites_used / n_unmasked_sites)
# By-Library normalizing is needed to compare/combine values across libraries with 
#  different overall depths. Performed on unmasked depth average from above.
#  For each sequence:
#    avg_norm_SEQ_depth = avg_raw_SEQ_depth / avg_raw_LIBRARY_depth
#  since the "raw" depths are themselves averages across all unmasked sites on sequence 
#  the library average needs to be weighted by unmasked sequence length:
#     LIBRARY_depth = sum( SCAFi_depth * SCAFi_unmasked ) / sum( SCAFi_unmasked )
# Takes: 
#  df.depth dataframe with seqID row names, $n_sites (number sites in average), and $avg_raw_depth
#  df.nmasked dataframe with seqID row names, $length (total bp), and $unmasked (non-N bp)
# Returns:
#  df.depth dataframe with seqID row names, $avg_raw_depth $unmasked $avg_unmask_depth $avg_norm_depth
# RETURN DF ROW NAMES TAKEN FROM DF_NMASK WITH MISSING VALUES SET TO 0 so make sure they match
normalize_depth <- function(DF_RAW_DEPTH, DF_NMASK){
  DF_NORM_DEPTH <- data.frame(
    unmasked = DF_NMASK$unmasked,
    n_sites = DF_RAW_DEPTH[rownames(DF_NMASK),]$n_sites,
    avg_raw_depth = DF_RAW_DEPTH[rownames(DF_NMASK),]$avg_raw_depth,
    row.names = row.names(DF_NMASK)
  )
  # missing scaffolds are treated as "no reads mapped"
  # unmasked, n_sites, and depth all set to 0
  DF_NORM_DEPTH[is.na(DF_NORM_DEPTH)] = 0
  # adjust for 0 artifacts
  DF_NORM_DEPTH$avg_unmask_depth = 
    DF_NORM_DEPTH$avg_raw_depth * (DF_NORM_DEPTH$n_sites / DF_NORM_DEPTH$unmasked)
  AVG_DEPTH = weighted.mean(DF_NORM_DEPTH$avg_unmask_depth, DF_NORM_DEPTH$unmasked)
    #sum( DF_NORM_DEPTH$avg_unmask_depth * DF_NORM_DEPTH$unmasked , na.rm = T) /
     #           sum( DF_NORM_DEPTH$unmasked , na.rm = T)
  DF_NORM_DEPTH$avg_norm_depth = 
    DF_NORM_DEPTH$avg_unmask_depth / AVG_DEPTH
  return(DF_NORM_DEPTH)
}

# make_df.sexbias() function to build the PER SCAFFOLD df.sexbias dataframe 
#   from multiple reps per sex (adapted from make_df.sexbias_NoReps(), below)
# Takes depth and variant calls per scaffold (see above for formatting)
#   read mapping against unmasked or hardmasked assembly,
#   need "nmasked.csv" either way, see above
# Does by-library normalizations, averages across libraries,
#  and calculates per scaffold statistics
# 2 arguments:
#  Filename of ASM.NMASK.CSV 
#  List (of lists of lists) with library depth and variant data
#    Format: LISTNAME$STAT$SEX$SAMPLE
make_df.sexbias <- function(FILE_SCAF_NMASKED, L_INFILES){
  
  # df for scaffold length, masked/unmasked bp per scaffold
  DF_SCAFFOLDS <- read.csv(FILE_SCAF_NMASKED, header = F,
                      col.names = c('scaffold', 'length', 'masked'))
  row.names(DF_SCAFFOLDS) = DF_SCAFFOLDS$scaffold  # useful to be able to access by either name or column value
  DF_SCAFFOLDS$unmasked <- DF_SCAFFOLDS$length - DF_SCAFFOLDS$masked
  # proportion of total assembly
  DF_SCAFFOLDS$pr_length <- DF_SCAFFOLDS$length / sum(DF_SCAFFOLDS$length)
  DF_SCAFFOLDS$pr_unmasked <- DF_SCAFFOLDS$unmasked / sum(DF_SCAFFOLDS$unmasked)
  # order by scaffold length
  DF_SCAFFOLDS <- DF_SCAFFOLDS[order(-DF_SCAFFOLDS$length),]
  
  # long df with sex, sample, scaffold id and lengths, raw and normalized depth
  DF_DEPTH <- lapply(L_INFILES$depth, function(sex){
    lapply(sex, function(FILE){
      # call normalize_depth() function from above while reading file
      sample_depth <- normalize_depth(
        read.table(FILE, header = F, row.names = 1, 
                   col.names = c('scaffold', 'n_sites', 'avg_raw_depth')), 
        DF_SCAFFOLDS) %>% arrange(desc(n_sites))
      sample_depth$scaffold <- row.names(sample_depth)
      rownames(sample_depth) = NULL
      return(sample_depth)
    }) %>% bind_rows(.id ="sample")
  }) %>% bind_rows(.id="sex") %>% 
    select(sex, sample, scaffold, n_sites, unmasked,
           avg_raw_depth, avg_unmask_depth, avg_norm_depth)
  
  # variant count long df with sex, sample, scaffold id and lengths
  # just a count so not as complicated as depth (no unmasked average and normalizing needed)
  DF_NVARI <- lapply(L_INFILES$nvari, function(sex){
    lapply(sex, function(FILE){
      # call normalize_depth() function from above while reading file
      df.vari <- read.csv(FILE, header = F, row.names = 1, 
                          col.names = c('scaffold','n_variants'))
      sample_nvari <- data.frame(
        unmasked = DF_SCAFFOLDS$unmasked,
        n_variants = df.vari[rownames(DF_SCAFFOLDS),],
        row.names = row.names(DF_SCAFFOLDS)
      ) %>% arrange(desc(unmasked))
      # still set unobserved scaffolds to "no variants"
      sample_nvari[is.na(sample_nvari)] = 0
      sample_nvari$scaffold <- row.names(sample_nvari)
      rownames(sample_nvari) = NULL
      return(sample_nvari)
    }) %>% bind_rows(.id ="sample")
  }) %>% bind_rows(.id="sex") %>% 
    select(sex, sample, scaffold, unmasked, n_variants)
  # per kb variant frequence
  DF_NVARI$fr_var = 1000 * DF_NVARI$n_variants / DF_NVARI$unmasked
  
  # wide DF to collect all values by scaffold
  # sample depth data
  DF_DEPTH_WIDE <- DF_DEPTH %>% select(-sex, -unmasked, -n_sites) %>% 
    pivot_wider(names_from = sample, 
                values_from = c(avg_raw_depth, avg_unmask_depth, avg_norm_depth)) %>% 
    as.data.frame()
  # sample variant data
  DF_NVARI_WIDE <- DF_NVARI %>% select(-sex, -unmasked) %>%
    pivot_wider(names_from = sample,
                values_from = c(n_variants, fr_var)) %>%
    as.data.frame()
  # list of depth statistics by sex
  L_DEPTH_STATS <- lapply(names(L_INFILES$depth), function(SEX){ 
    depth_stats <- DF_DEPTH %>% filter(sex == SEX) %>% 
      group_by(scaffold) %>% 
      summarise(mean_depth = mean(avg_norm_depth), 
                sd_depth = sd(avg_norm_depth) ) %>% as.data.frame()
  names(depth_stats) = c(names(depth_stats)[1], paste(names(depth_stats), SEX, sep='.')[2:3])
    rownames(depth_stats) = depth_stats$scaffold
    return(depth_stats)
  })
  # list of variant stats by sex
  L_NVARI_STATS <- lapply(names(L_INFILES$nvari), function(SEX){
    nvari_stats <- DF_NVARI %>% filter(sex == SEX) %>% 
      group_by(scaffold) %>% 
      summarise(mean_frVar = mean(fr_var), 
                sd_frVar = sd(fr_var),
                mean_nVar = mean(n_variants),
                sd_nVar = sd(n_variants)) %>% as.data.frame()
    names(nvari_stats) = c(names(nvari_stats)[1], paste(names(nvari_stats), SEX, sep='.')[2:5] )
    rownames(nvari_stats) = nvari_stats$scaffold
    return(nvari_stats)
  })
  
  # combine by columns, use merge to drop redundant columns
  DF_WIDE <- merge(
    DF_SCAFFOLDS,
    merge(
      merge(DF_DEPTH_WIDE, DF_NVARI_WIDE),
      merge(
        merge(L_DEPTH_STATS[[1]], L_DEPTH_STATS[[2]]),
        merge(L_NVARI_STATS[[1]], L_NVARI_STATS[[2]])
      )
    )
  ) %>% arrange(desc(length)) %>% as.data.frame()
  
  # top level bias stats, using list positions to avoid hard coding sex ids 
  DF_WIDE$Depth_bias <- L_DEPTH_STATS[[1]][DF_WIDE$scaffold, 2] / L_DEPTH_STATS[[2]][DF_WIDE$scaffold, 2]
  DF_WIDE$FrVar_bias <- L_NVARI_STATS[[1]][DF_WIDE$scaffold, 2] / L_NVARI_STATS[[2]][DF_WIDE$scaffold, 2]
  
  rownames(DF_WIDE) = DF_WIDE$scaffold
  
  L_SEXBIAS <- list(
    SexBias = DF_WIDE,
    DataDepth = DF_DEPTH,
    DataVariants = DF_NVARI
  )
  
  return(L_SEXBIAS)
  
}

# make_df.sexbias-NoReps() a function to build the df.sexbias dataframe
# depth and variant calls from read mapping against hardmasked assembly
# Initial version of function written for single HiFi read set per sex (NoReps)
# takes a list of 5 filenames, all generated by Ceres scripts: 
#   nmasked: nmasked.csv file with per-scaffold lengths and masked bp counts
#   mean_depth.female: depth_per_scaffold.tsv from female reads against hardmasked assembly
#   mean_depth.male: depth_per_scaffold.tsv from male reads against hardmasked assembly
#     depth_per_scaffold.tsv w/ 3 columns: scaffold n_sites depth -- via get_depth_scaffold.sh
#     (needed to account for ambiguity reporting in masked or depth=0 sites)
#   nvariF: nseg_per_scaffold.csv from female reads, count of variant sites of hardmasked scaffolds
#   nvariM: nseg_per_scaffold.csv from male reads, count of variant sites of hardmasked scaffolds
# reads in all files, sorts by scaffold length, calculates for each scaffold:
#   - total length, masked bp, unmasked length
#   - proportion of total assembly length
#   - proportion of unmasked assembly length
#   - female and male normalized depths (to unmasked length-weighted average depth)
#   - female and male frequency of variant sites
#   - F/M ratio of normalized depth
#   - F/M ratio of variant frequency 
make_df.sexbias_NoReps <- function(l.IN_FILES){
  # get scaffolds and lengths from asm fasta index
  NMASKED <- read.csv(l.IN_FILES$nmasked, header = F,
                      col.names = c('scaffold', 'length', 'masked'))
  NMASKED$unmasked <- NMASKED$length - NMASKED$masked
  NMASKED <- NMASKED[order(-NMASKED$unmasked),]
  # input dataframes
  IN_DEPTHF <- read.table(l.IN_FILES$mean_depth.female, header = F, row.names = 1, 
                          col.names = c('scaffold', 'n_sites', 'depth'))
  IN_DEPTHM <- read.table(l.IN_FILES$mean_depth.male, header = F, row.names = 1, 
                          col.names = c('scaffold', 'n_sites', 'depth'))
  IN_NVARIF <- read.csv(l.IN_FILES$nvariF, header = F, row.names = 1, 
                        col.names = c('scaffold','n_variants'))
  IN_NVARIM <- read.csv(l.IN_FILES$nvariM, header = F, row.names = 1, 
                        col.names = c('scaffold','n_variants'))
  # scaffolds ordered by length
  SCAFFOLDS <- NMASKED$scaffold
  # start output df ordered by scaffold length
  SEXBIAS <- data.frame(length = NMASKED$length,
                        masked = NMASKED$masked,
                        unmasked = NMASKED$unmasked,
                        nsitesF = IN_DEPTHF[SCAFFOLDS,]$n_sites,
                        mean_depth.female = IN_DEPTHF[SCAFFOLDS,]$depth,
                        nsitesM = IN_DEPTHM[SCAFFOLDS,]$n_sites,
                        mean_depth.male = IN_DEPTHM[SCAFFOLDS,]$depth,
                        mean_nVar.female = IN_NVARIF[SCAFFOLDS,],
                        mean_nVar.male = IN_NVARIM[SCAFFOLDS,],
                        row.names = SCAFFOLDS)
  # treat missing values as 0
  SEXBIAS[is.na(SEXBIAS)] = 0
  # proportion of total length per scaffold
  SEXBIAS$pr_length <- SEXBIAS$length / sum(SEXBIAS$length)
  # proportion of unmasked length per scaffold
  SEXBIAS$pr_unmasked <- SEXBIAS$unmasked / sum(SEXBIAS$unmasked)
  # unmasked length-weighted normalized depths
  SEXBIAS$depthF_n <- SEXBIAS$mean_depth.female / sum(SEXBIAS$mean_depth.female * SEXBIAS$pr_unmasked)
  SEXBIAS$depthM_n <- SEXBIAS$mean_depth.male / sum(SEXBIAS$mean_depth.male * SEXBIAS$pr_unmasked)
  # frequency of variant sites per kb
  SEXBIAS$fr_varF <- 1000 * SEXBIAS$mean_nVar.female / SEXBIAS$unmasked
  SEXBIAS$fr_varM <- 1000 * SEXBIAS$mean_nVar.male / SEXBIAS$unmasked
  # frequency of variants per kb, normalized to library average
  #SEXBIAS$fr_varF_n <- SEXBIAS$fr_varF / (1000 * sum(SEXBIAS$mean_nVar.female)/sum(SEXBIAS$unmasked))
  #SEXBIAS$fr_varM_n <- SEXBIAS$fr_varM / (1000 * sum(SEXBIAS$mean_nVar.male)/sum(SEXBIAS$unmasked))
  # sex bias (female / male)
  SEXBIAS$Depth_bias <- SEXBIAS$depthF_n / SEXBIAS$depthM_n
  SEXBIAS$FrVar_bias <- SEXBIAS$fr_varF / SEXBIAS$fr_varM
  
  return(SEXBIAS)
  
}


# Function to combine stats from df.SexBias across list of scaffolds
# Used in make_df.HomologSexBias() inside sapply call
ClusterStats <- function(SCAF_LIST, DF.SEXBIAS){
  cbind( nscafs = length(SCAF_LIST), 
         longest = SCAF_LIST[order(-DF.SEXBIAS[SCAF_LIST,]$length)][1],
         DF.SEXBIAS[SCAF_LIST,] %>% 
           summarise(total.bp = sum(length), 
                     unmasked.bp = sum(unmasked),
                     pr_length = sum(pr_length),
                     pr_unmasked = sum(pr_unmasked),
                     # length-weighted average across scaffolds
                     depth.F = weighted.mean(mean_depth.female, unmasked),
                     depth.M = weighted.mean(mean_depth.male, unmasked),
                     # frequency across scaffolds
                     var_count.F = sum(mean_nVar.female), 
                     var_count.M = sum(mean_nVar.male),
                     var_fr.F = sum(mean_nVar.female)/sum(unmasked),
                     var_fr.M = sum(mean_nVar.male)/sum(unmasked)
                     #depth_n.F = sum(unmasked*depthF_n/sum(unmasked), na.rm = T), 
                     #depth_n.M = sum(unmasked*depthM_n/sum(unmasked), na.rm = T), 
                     #var_n.F = sum(unmasked*fr_varF_n/sum(unmasked), na.rm = T), 
                     #var_n.M = sum(unmasked*fr_varM_n/sum(unmasked), na.rm = T)
  ) )
}

# function to make df with composite stats for homolog clusters
# takes homolog cluster list from alignment_region_dotplot-FUNCTIONS.R GetHomologs()
#  (my convention in those alignments is ref=M qry=F)
# takes DF.SEXBIAS for each sex from make_df.sexbias
make_df.HomologSexBias <- function(CLUSTERS, DF.SEXBIAS.F, DF.SEXBIAS.M){
  
  DF.HOMOLOGS.F <- sapply(CLUSTERS, function(x){ 
      ClusterStats(x$QryScafs, DF.SEXBIAS.F) 
    }) %>% t() %>% data.frame()
  
  DF.HOMOLOGS.M <- sapply(CLUSTERS, function(x){
      ClusterStats(x$RefScafs, DF.SEXBIAS.M )
    }) %>% t() %>% data.frame()
  
  DF.HOMOLOGS <- data.frame(
    F_longest = unlist(DF.HOMOLOGS.F$longest),
    F_nscafs = unlist(DF.HOMOLOGS.F$nscafs),
    F_total.bp = unlist(DF.HOMOLOGS.F$total.bp),
    F_unmasked = unlist(DF.HOMOLOGS.F$unmasked.bp),
    F_pr.length = unlist(DF.HOMOLOGS.F$pr_length),
    F_pr.unmasked = unlist(DF.HOMOLOGS.F$pr_unmasked),
    M_longest = unlist(DF.HOMOLOGS.M$longest),
    M_nscafs = unlist(DF.HOMOLOGS.M$nscafs),
    M_total.bp = unlist(DF.HOMOLOGS.M$total.bp),
    M_unmasked = unlist(DF.HOMOLOGS.M$unmasked),
    M_pr.length = unlist(DF.HOMOLOGS.M$pr_length),
    M_pr.unmasked = unlist(DF.HOMOLOGS.M$pr_unmasked),

    #F_depth.raw.F = unlist(DF.HOMOLOGS.F$depth_raw.F),
    #F_depth.raw.M = unlist(DF.HOMOLOGS.F$depth_raw.M),
    Fasm_depth.F = unlist(DF.HOMOLOGS.F$depth.F),
    Fasm_depth.M = unlist(DF.HOMOLOGS.F$depth.M),
    #M_depth.raw.F = unlist(DF.HOMOLOGS.M$depth_raw.F),
    #M_depth.raw.M = unlist(DF.HOMOLOGS.M$depth_raw.M),
    Masm_depth.F = unlist(DF.HOMOLOGS.M$depth.F),
    Masm_depth.M = unlist(DF.HOMOLOGS.M$depth.M),
    
    F_var.count.F = unlist(DF.HOMOLOGS.F$var_count.F),
    F_var.count.M = unlist(DF.HOMOLOGS.F$var_count.M),
    F_var.fr.F = unlist(DF.HOMOLOGS.F$var_fr.F),
    F_var.fr.M = unlist(DF.HOMOLOGS.F$var_fr.M),
    #Fasm_var.F = unlist(DF.HOMOLOGS.F$var_n.F),
    #Fasm_var.M = unlist(DF.HOMOLOGS.F$var_n.M),
    M_var.count.F = unlist(DF.HOMOLOGS.M$var_count.F),
    M_var.count.M = unlist(DF.HOMOLOGS.M$var_count.M),
    M_var.fr.F = unlist(DF.HOMOLOGS.M$var_fr.F),
    M_var.fr.M = unlist(DF.HOMOLOGS.M$var_fr.M),
    #Masm_var.F = unlist(DF.HOMOLOGS.M$var_n.F),
    #Masm_var.M = unlist(DF.HOMOLOGS.M$var_n.M),
    
    DepthBias_Fasm = (unlist(DF.HOMOLOGS.F$depth.F) + 0.001) /
                      (unlist(DF.HOMOLOGS.F$depth.M) + 0.001),
    DepthBias_Masm = (unlist(DF.HOMOLOGS.M$depth.F) + 0.001) / 
                      (unlist(DF.HOMOLOGS.M$depth.M) + 0.001),
    VarBias_Fasm = (unlist(DF.HOMOLOGS.F$var_fr.F) + 0.001) / 
                    (unlist(DF.HOMOLOGS.F$var_fr.M) + 0.001),
    VarBias_Masm = (unlist(DF.HOMOLOGS.M$var_fr.F) + 0.001) / 
                    (unlist(DF.HOMOLOGS.M$var_fr.M) + 0.001)
  )
  
  DF.HOMOLOGS$AvgTotalLength <- apply(
    DF.HOMOLOGS[,c('F_total.bp','M_total.bp')], 1, mean, na.rm=T )
  DF.HOMOLOGS$AvgUnmasked <- apply(
    DF.HOMOLOGS[,c('F_unmasked','M_unmasked')], 1, mean, na.rm=T )
  DF.HOMOLOGS$LengthBias <- DF.HOMOLOGS$F_total.bp / DF.HOMOLOGS$M_total.bp
  DF.HOMOLOGS$UnmaskedBias <- DF.HOMOLOGS$F_unmasked / DF.HOMOLOGS$M_unmasked
  DF.HOMOLOGS$AvgDepthBias <- apply(
    DF.HOMOLOGS[,c('DepthBias_Fasm','DepthBias_Masm')], 1, mean, na.rm=T )
  DF.HOMOLOGS$AvgVarBias <- apply(
    DF.HOMOLOGS[,c('VarBias_Fasm','VarBias_Masm')], 1, mean, na.rm=T )
  
  return(DF.HOMOLOGS[order(-DF.HOMOLOGS$AvgTotalLength),])
}


# function to link depth windows to unmasked lengths and 
#  normalize depths by unmasked-weighted average
make_df.DepthWindows <- function(INFILE_DEPTH, INFILE_NMASKED,
                                 AVG_CLIP_PCTILE = 0){
  # read in csv with total and masked bases per window, calculate unmasked per window
  DF.NMASKED <- read.csv(INFILE_NMASKED, header = F, row.names = 1, 
                          col.names = c('window', 'total', 'masked'))
  DF.NMASKED$unmasked <- DF.NMASKED$total - DF.NMASKED$masked
  
  # read in tsv with average depths per window and number sites used
  DF.DEPTH_IN <- read.csv(INFILE_DEPTH, sep = '\t', header = F,
                       col.names = c('scaffold', 'position',
                                     'window', 'nsites', 'depth'))
  row.names(DF.DEPTH_IN) = DF.DEPTH_IN$window
  # build depth df, treating unobserved windows (present in DF.NMASK) as depth 0
  DF.DEPTH <- data.frame(
    scaffold = rownames(DF.NMASKED) %>% str_split_i(pattern = ":", 1),
    position = rownames(DF.NMASKED) %>% str_split_i(pattern = ":", 2),
    window = row.names(DF.NMASKED),
    nsites = DF.DEPTH_IN[row.names(DF.NMASKED),]$nsites,
    depth = DF.DEPTH_IN[row.names(DF.NMASKED),]$depth,
    unmasked = DF.NMASKED$unmasked,
    row.names = row.names(DF.NMASKED)
  )
  # missing values to 0
  DF.DEPTH[is.na(DF.DEPTH)] = 0
  # account for false or missing 0s in depth average
  DF.DEPTH <- DF.DEPTH %>% filter(unmasked > 0) %>% 
    mutate(depth_w0s = depth * nsites/unmasked, depth)
  # length-weighted average depth
  # subset of windows for average, can clip outliers by setting AVG_CLIP_PCTILE > 0
  DF.DEPTH.NoHighLow <- DF.DEPTH %>% filter(depth_w0s > quantile(.$depth_w0s, AVG_CLIP_PCTILE),
                                            depth_w0s < quantile(.$depth_w0s, (1-AVG_CLIP_PCTILE)))
  DEPTH_AVG = sum(DF.DEPTH.NoHighLow$depth_w0s 
                  * DF.DEPTH.NoHighLow$unmasked/sum(DF.DEPTH.NoHighLow$unmasked))
  #DEPTH_AVG = sum(DF.DEPTH$depth_w0s 
  #                * DF.DEPTH$unmasked/sum(DF.DEPTH$unmasked))
  # normalize by length-weighed average
  DF.DEPTH <- DF.DEPTH %>% mutate(depth_norm = depth_w0s / DEPTH_AVG)
  
  return(DF.DEPTH)
}


# function to link variant count windows to unmasked lengths 
#  and find variant frequencies per unmasked bp
make_df.VariantWindows <- function(INFILE_VARIANT, INFILE_NMASKED){
  # read in csv with total and masked bases per window, calculate unmasked per window
  DF.NMASKED <- read.csv(INFILE_NMASKED, header = F, row.names = 1, 
                         col.names = c('window', 'total', 'masked'))
  DF.NMASKED$unmasked <- DF.NMASKED$total - DF.NMASKED$masked
  
  # read in tsv with variant site count per window 
  DF.VARIANT <- read.csv(INFILE_VARIANT, sep = '\t',
                       col.names = c('scaffold', 'position',
                                     'window', 'n_variants'))
  # add unmasked and total lengths for each window from NMASKED df
  DF.VARIANT$unmasked <- DF.NMASKED[DF.VARIANT$window,]$unmasked
  DF.VARIANT$total <- DF.NMASKED[DF.VARIANT$window,]$total
  # frequency variant sites per window
  DF.VARIANT <- DF.VARIANT %>% mutate(fr_var = n_variants / unmasked)

  return(DF.VARIANT)
}


# Functions for plotting windows data
# make the dataframe for plotting 
make_df.WindowPlots <- function(DF_WINDOWS, SEQ_LIST, MIN_UNMASKED){
  PLOTDF <- DF_WINDOWS %>% filter(scaffold %in% SEQ_LIST, unmasked >= MIN_UNMASKED)
  OFFSETS <- PLOTDF %>% group_by(scaffold) %>% 
    # highest value per scaffold, add 1 for plot spacing
    summarise(lastpos = max(as.numeric(position)) + 1) %>% arrange(match(scaffold, SEQ_LIST)) #%>% arrange(desc(lastpos)) # arrange by size order not stated order 
  OFFSETS$offset <- (c(0,OFFSETS$lastpos[1:nrow(OFFSETS)-1]) %>% cumsum())
  PLOTDF$offset <- apply(PLOTDF, 1, function(x){
    OFFSETS$offset[OFFSETS$scaffold == x[1]]
  })
  return(PLOTDF)
}

PlotDepthWindows <- function(DF_WINDOWS, SEQ_LIST, MIN_UNMASKED, DROP_OUTLIERS = T,
                             POOL_N = 1, # can pool windows to speed up plotting
                             ASSEMBLY_NAME, WINDOW_LABEL,
                             SPAN=0.4, ALPHA=0.05, 
                             Y_MIN='auto', Y_MAX='auto', FLIP_SCAF_LAB = F){
  # get dataframe for plotting
  PLOTDF <- make_df.WindowPlots(DF_WINDOWS, SEQ_LIST, MIN_UNMASKED)
  
  if( DROP_OUTLIERS ){
    PLOTDF <- PLOTDF  %>% group_by(scaffold, rep) %>% 
      filter(!(abs(depth_norm - median(depth_norm)) > 2*sd(depth_norm)) 
             | n() == 1) # keep single observations (were getting dropped since sd NA)
  }
  
  if( POOL_N > 1 ){
    PLOTDF <- PLOTDF %>% 
      mutate(pool = as.integer(as.numeric(position) / POOL_N)) %>% 
      group_by(scaffold, pool, sex, rep) %>% 
      summarise(depth_norm = weighted.mean(depth_norm, unmasked), 
                position = median(as.numeric(position)), 
                offset = max(offset) )
    WINDOW_LABEL = paste0(POOL_N, " pooled ", WINDOW_LABEL)
  }
  
  # set color scale for sex based on how many sexes represented
  N_SEX = PLOTDF$sex %>% factor() %>% levels() %>% length()
  if( N_SEX == 2 ){
    SCALE_VALUES = c('red', 'blue')
  } else {
    SCALE_VALUES=c('slateblue')
  }
  
  # set Y min and max based on PLOTDF unless set by command call
  if(Y_MIN == 'auto'){
    Y_MIN <- min(PLOTDF$depth_norm)
  }
  if(Y_MAX == 'auto'){
    Y_MAX <- max(PLOTDF$depth_norm)
  }
  
  if(FLIP_SCAF_LAB){
    ANNO_Y = Y_MAX
    ANNO_ANGLE = 90
    ANNO_HJUST = 1
    ANNO_VJUST = 1.1
  } else {
    ANNO_Y = Y_MIN
    ANNO_ANGLE = 0
    ANNO_HJUST = -0.1
    ANNO_VJUST = 1
  }
  
  # get names and x-positions for scaffold annotation
  OFFSETS <- PLOTDF %>% group_by(scaffold) %>% summarise(offset=max(offset))
  
  PLOT <- ggplot(PLOTDF, aes(x = as.numeric(position) + offset + 0.5, 
                             y = depth_norm, 
                             color = sex, 
                             shape = rep, linetype = rep, 
                             lineend = scaffold)) +
    # point cloud
    geom_point(alpha=ALPHA) + 
    # smoothed lines
    geom_smooth(se=F, method = 'loess', span=SPAN) +
    # improve legibility of legend for rep
    guides(shape = guide_legend(override.aes = list(alpha=0.6, size=3, color='black')) ) +
    # colors for sex (set above)
    scale_color_manual(values = SCALE_VALUES) +
    # lines and labels for sequences
    geom_vline(xintercept = OFFSETS$offset, color='black', lty=2, linewidth=0.8) +
    annotate("text", label=OFFSETS$scaffold, x=OFFSETS$offset, 
             y = ANNO_Y, angle = ANNO_ANGLE, 
             hjust = ANNO_HJUST, vjust = ANNO_VJUST) +
    # other plot formatting
    coord_cartesian(ylim = c(Y_MIN,Y_MAX)) + 
    theme_light() +
    labs( y = 'normalized depth', 
          x = paste0(ASSEMBLY_NAME, " - ", WINDOW_LABEL, " windows, >", 
                   MIN_UNMASKED, " unmasked bp, loess span=", SPAN) )
  
  return(PLOT)

}


PlotVariantWindows <- function(DF_WINDOWS, SEQ_LIST, MIN_UNMASKED, 
                             ASSEMBLY_NAME, WINDOW_LABEL,
                             POOL_N = 1,
                             SPAN=0.4, ALPHA=0.05, 
                             Y_MIN='auto', Y_MAX='auto', FLIP_SCAF_LAB = F){
  # get dataframe for plotting
  PLOTDF <- make_df.WindowPlots(DF_WINDOWS, SEQ_LIST, MIN_UNMASKED)
  
  if( POOL_N > 1 ){
    PLOTDF <- PLOTDF %>% 
      mutate(pool = as.integer(as.numeric(position) / POOL_N)) %>% 
      group_by(scaffold, pool, sex, rep) %>% 
      summarise(fr_var = weighted.mean(fr_var, unmasked), 
                position = median(as.numeric(position)), 
                offset = max(offset) )
    WINDOW_LABEL = paste0(POOL_N, " pooled ", WINDOW_LABEL)
  }

    # set color scale for sex based on how many sexes represented
  N_SEX = PLOTDF$sex %>% factor() %>% levels() %>% length()
  if( N_SEX == 2 ){
    SCALE_VALUES = c('red', 'blue')
  } else {
    SCALE_VALUES=c('slateblue')
  }

  # set Y min and max based on PLOTDF unless set by command call
  if(Y_MIN == 'auto'){
    Y_MIN <- min(PLOTDF$fr_var)
  }
  if(Y_MAX == 'auto'){
    Y_MAX <- max(PLOTDF$fr_var)
  }
  
  if(FLIP_SCAF_LAB){
    ANNO_Y = Y_MAX
    ANNO_ANGLE = 90
    ANNO_HJUST = 1
    ANNO_VJUST = 1.1
  } else {
    ANNO_Y = Y_MIN
    ANNO_ANGLE = 0
    ANNO_HJUST = -0.1
    ANNO_VJUST = 1
  }
  
  # get names and x-positions for scaffold annotation
  OFFSETS <- PLOTDF %>% group_by(scaffold) %>% summarise(offset=max(offset))
  
  PLOT <- ggplot(PLOTDF, aes(x = position+offset + 0.5, 
                             y = fr_var, 
                             color = sex, 
                             shape = rep, linetype = rep, 
                             lineend = scaffold)) +
    # point cloud
    geom_point(alpha=ALPHA) + 
    # smoothed lines
    geom_smooth(se=F, method = 'loess', span=SPAN) +
    # improve legibility of legend for rep
    guides(shape = guide_legend(override.aes = list(alpha=0.6, size=3, color='black')) ) +
    # colors for sex (set above)
    scale_color_manual(values = SCALE_VALUES) +
    # lines and labels for sequences
    geom_vline(xintercept = OFFSETS$offset, color='black', lty=2, linewidth=0.8) +
    annotate("text", label=OFFSETS$scaffold, x=OFFSETS$offset, 
             y = ANNO_Y, angle = ANNO_ANGLE, 
             hjust = ANNO_HJUST, vjust = ANNO_VJUST) +
    # other plot formatting
    coord_cartesian(ylim = c(Y_MIN,Y_MAX)) + 
    theme_light() +
    labs( y = 'frequency variant sites', 
          x = paste0(ASSEMBLY_NAME, " - ", WINDOW_LABEL, " windows, >", 
                     MIN_UNMASKED, " unmasked bp, loess span=", SPAN) )
  
  return(PLOT)
  
}




