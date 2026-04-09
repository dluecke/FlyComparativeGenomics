# rephase_XY-FUNCTIONS.R 

# functions used in the rephase_XY scripts

library(tidyverse)


# ALIGNMENT STATS
# USE MALE DIPLOID NON-AUTOSOME CTGS ALIGNED TO FEMALE DIPLOID SCAFFOLDS

# function to report fraction of query sequence not included in any alignment
# takes Q1 and Q2 coordinate vectors and full qry sequence length
# returns proportion [0-1] of query not covered by any alignment
p_qry_unaligned = function(vQ1, vQ2, QLEN){
  UNCOVERED = c(1:QLEN)
  for(i in seq_along(vQ1)){
    ALNi = c(vQ1[i]:vQ2[i])
    UNCOVERED = UNCOVERED[!(UNCOVERED %in% ALNi)]
  }
  p_UNCOVERED = length(UNCOVERED)/QLEN
  return(p_UNCOVERED)
}

# ALIGNMENT STAT 1
# function to get per-contig alignment-based Yscore (how much of male contig aligns to female assembly)
# takes l.coords list (with coords and seq length dfs) of alignment vs female
make_df.align.pct_absent = function(L.COORDS.F){
  df.ALIGN.Yscore = L.COORDS.F$coords %>% 
    # add qry sequence length and percent of qry covered by alignment
    mutate(qry_SeqLength = L.COORDS.F$breaksQry[qry_SeqName,2], 
           qry_cov = qry_length/qry_SeqLength) %>%
    group_by(qry_SeqName) %>% 
    summarise(n=n(), # number alignments per qry
              seqlength = qry_SeqLength[1], # qry length
              # length-weighted average of pctID across aligned regions (doesn't account for unaligned qry sequence)
              avg_pctID = weighted.mean(pctID, qry_length),
              # fraction total query length not covered by alignment
              p_qry_unaln = p_qry_unaligned(Q1, Q2, seqlength), 
              # length-weighted average of pctID, scaled by fraction of query in alignments
              # this isn't exactly right due to alignment overlap, but pretty close
              scaled_pctID = avg_pctID * (1 - p_qry_unaln),
              #total_pctID = mean(pctID * qry_cov) * n, 
              # estimate of percent bp not in female reference
              pctABSENT = (100 - scaled_pctID)) %>% 
    ungroup %>% 
    # add rows for all query sequences not in alignment, treat as fully absent
    add_row(qry_SeqName = L.COORDS.F$breaksQry$SeqName[!(L.COORDS.F$breaksQry$SeqName %in% .$qry_SeqName)],
            n = 0,
            seqlength = L.COORDS.F$breaksQry[qry_SeqName,2],
            p_qry_unaln = 1,
            avg_pctID = 0,
            pctABSENT = 100) %>%
    # Yscore_aln is just pctABSENT - alignment length-weighted average 
    # Yscore=0 means ctg at average, Yscore>0 means more diverged from Fp than average), 
    mutate(Yscore_aln = (pctABSENT - weighted.mean(pctABSENT, seqlength))
           # need to scale to var=1 before combining with other metrics, but decided to keep this version "raw" 
           # will scale variance at the metric combination step
           #, Yscore_aln = Yscore_aln/sqrt(var(Yscore_aln))
    ) %>%
    as.data.frame()
  row.names(df.ALIGN.Yscore) = df.ALIGN.Yscore$qry_SeqName
  return(df.ALIGN.Yscore)
}

# ALIGNMENT STAT 2
# scores contigs for how split/fragmented the alignment is across reference sequences
# large discontinuity between contig and reference may point at rearrangement, possible XY divergence
make_df.align.ctg_splits = function(L.COORDS, N_AUTOSOME=5, PLOIDY=2){
  CHRS = L.COORDS$breaksRef %>% arrange(desc(SeqLength)) %>% select(SeqName) %>% head(n=N_AUTOSOME*PLOIDY) %>% unlist %>% as.character()
  df.SPLITS = L.COORDS$coords %>% 
    filter( ! ref_SeqName %in% CHRS ) %>%
    mutate(qry_SeqLength = L.COORDS$breaksQry[qry_SeqName,2],
           ref_SeqLength = L.COORDS$breaksRef[ref_SeqName,2]) %>%
    group_by(qry_SeqName, ref_SeqName) %>%
    summarise(
      qry_SeqLength = qry_SeqLength[1],
      ref_SeqLength = ref_SeqLength[1],
      sum_qry_len = sum(qry_length),
      minQ = min( min(Q1), min(Q2) ),
      maxQ = max( max(Q1), max(Q2) ),
      Q_span = maxQ - minQ,
      # fraction total query length not covered by alignment
      p_qry_aln = 1 - p_qry_unaligned(Q1, Q2, qry_SeqLength), 
      minR = min(R1),
      maxR = max(R2),
      R_span = maxR - minR,
      Q_split = abs(R_span - Q_span)/qry_SeqLength,
      log_Q_splt = log2(Q_split + 0.001)
    ) %>% ungroup() %>% 
    group_by(qry_SeqName) %>%
    slice_max(p_qry_aln, n = 1, with_ties = F) %>%
    as.data.frame()
  rownames(df.SPLITS) = df.SPLITS$qry_SeqName
  return(df.SPLITS)
}

# COMBINE ALIGN STATS
# function to join df.align and df.q_split by row name (qry SeqName)
# Stats centered at 0 (X < 0, Y > 0)
#  "Yscore_aln" = pct(qry_absent) - wmean(qry_absent)
#  "log_Q_splt" = log2(Rspan/Qspan) vs female diploid scaffolds
make_df.align_yscores = function(DF.ALIGN, DF.SPLIT){
  # full join by row name
  df.ALIGN_YSCORES = full_join( DF.ALIGN, DF.SPLIT  ) %>%
    mutate(
      # scale variances to 1
      Yscore_PctAbsnt = Yscore_aln/sqrt(var(Yscore_aln, na.rm = T)),
      Yscore_logQsplt = log_Q_splt/sqrt(var(log_Q_splt, na.rm = T)),
      # sum (will scale again when combing with mapping/ymer stats)
      Yscore_align = rowSums(across(c(Yscore_PctAbsnt, Yscore_logQsplt)), na.rm = T)
    ) %>% as.data.frame() 
  row.names(df.ALIGN_YSCORES) = df.ALIGN_YSCORES$qry_SeqName
  return(df.ALIGN_YSCORES)
}


# DEPTH STATS

# function to output list with:
#  "long" (reps split)
#  "avgs" (rep normalization, averages)
#  "bias" (log2 ratio F to M normalized depth)
# input list of 2 lists of variant file: IN$female and IN$male
make_l.df.depth_by_seq = function(l.DEPTH_FILES){
  # tags used by rbind combine list (sex) and sublist (rep) names
  REPS.F = paste("female", names(l.DEPTH_FILES$female), sep = ".")
  REPS.M = paste("male", names(l.DEPTH_FILES$male), sep = ".")
  # read all input tables into combined long df
  df.long.DEPTH = lapply(
    l.DEPTH_FILES %>% unlist(), function(x){
      read.csv(x, sep='\t', header = F,
               col.names = c('sequence', 'length', 'depth'))  %>%
        mutate(norm_depth = depth/weighted.mean(depth, length))
    }
  ) %>% bind_rows(.id = "rep") %>%
    mutate(
      # assign sex based on input file ID ("rep")
      sex = case_when(
        rep %in% REPS.F ~ "female",
        rep %in% REPS.M ~ "male",
        TRUE ~ "unk"),
    )
  # convert long df into sequence/sex of reads, find average values
  df.DEPTH = df.long.DEPTH %>%
    group_by(sequence, sex) %>%
    # average number variants per sequence/sex 
    summarize(length = length[1],
              avg_norm_depth = mean(norm_depth)) %>%
    ungroup()
  # collapse df to 1 row per sequence
  df.BIAS = df.DEPTH %>% 
    group_by(sequence) %>%
    summarise(length = length[1],
              # average values from females, males
              norm_depth.F = mean( avg_norm_depth[sex == "female"] ),
              norm_depth.M = mean( avg_norm_depth[sex == "male"] )) %>%
    ungroup() %>%
    # NAs from sequences missing in sex of one species, replace with 0s
    mutate(norm_depth.F = if_else(is.na(norm_depth.F), 0, norm_depth.F),
           norm_depth.M = if_else(is.na(norm_depth.M), 0, norm_depth.M),
           # Yscore_depth is log ratio of Mval:Fval (Y chr means more depth in M reads)
           #  epsilon 0.01 to avoid issues at fr=0
           # score=0 for same in both sexes, score>0 for more female variants (y-like)
           # will need to scale to var=1 before combining with other metrics (X/sqrt(var(X)))
           Yscore_depth = log2( (norm_depth.M + 0.01) / (norm_depth.F + 0.01) ) ) %>% 
    as.data.frame()
  row.names(df.BIAS) = df.BIAS$sequence
  return(list(
    long = df.long.DEPTH,
    avgs = df.DEPTH,
    bias = df.BIAS
  ))
}


# VARIANT STATS

# function to output list with:
#  "long" (reps split)
#  "avgs" (rep avgs)
#  "bias" (log2 ratio F to M fr_var)
# inputs:
#   list of 2 lists of variant file: IN$female and IN$male
#   csv files with seqname,n_var,seqlength
make_l.df.var_by_ctg = function(l.N_VAR_FILES){
  # tags used by rbind combine list (sex) and sublist (rep) names
  REPS.F = paste("female", names(l.N_VAR_FILES$female), sep = ".")
  REPS.M = paste("male", names(l.N_VAR_FILES$male), sep = ".")
  # read all input tables into combined long df
  df.long.N_VARS = lapply(
    l.N_VAR_FILES %>% unlist(), function(x){
      read.csv(x, header = F,
               col.names = c('sequence','n_vrnts','seqlength'))      
    }
  ) %>% bind_rows(.id = "rep") %>%
    mutate(
      # replace "NA" with 0 variants
      n_vrnts = if_else(is.na(n_vrnts), 0, n_vrnts),
      # assign sex based on input file ID ("rep")
      sex = case_when(
        rep %in% REPS.F ~ "female",
        rep %in% REPS.M ~ "male",
        TRUE ~ "unk"),
      fr_vrnts = n_vrnts / seqlength
    )
  # convert long df into sequence/sex of reads, find average values
  df.VARS = df.long.N_VARS %>%
    group_by(sequence, sex) %>%
    # average number variants per sequence/sex 
    summarize(avg_vrnts = mean(n_vrnts),
              seqlength = seqlength[1]) %>%
    # frequency of variant average
    mutate(fr_vrnts = avg_vrnts / seqlength)
  # collapse df to 1 row per sequence
  df.BIAS = df.VARS %>% 
    group_by(sequence) %>%
    summarise(seqlength = seqlength[1],
              # average values from females, males
              fr_vrnts.F = mean( fr_vrnts[sex == "female"] ),
              fr_vrnts.M = mean( fr_vrnts[sex == "male"] )) %>%
    ungroup() %>%
    # NAs from sequences missing in sex of one species, replace with 0s
    mutate(fr_vrnts.F = if_else(is.na(fr_vrnts.F), 0, fr_vrnts.F),
           fr_vrnts.M = if_else(is.na(fr_vrnts.M), 0, fr_vrnts.M),
           # Yscore_vrnts is log ratio of Fval:Mval (Y chr means more varnts in F reads)
           #  epsilon 1e-5 (0.001%) to avoid issues at fr=0
           # score=0 for same in both sexes, score>0 for more female variants (y-like)
           # will need to scale to var=1 before combining with other metrics (X/sqrt(var(X)))
           Yscore_vrnts = log2( (fr_vrnts.F + 1e-5) / (fr_vrnts.M + 1e-5) ) ) %>% 
    as.data.frame()
  row.names(df.BIAS) = df.BIAS$sequence
  return(list(
    long = df.long.N_VARS,
    avgs = df.VARS,
    bias = df.BIAS
  ))
}


# YMER STATS

# function to read ymer coverage file and calculate yscore
make_df.ymer = function(COVFILE){
  # read in samtools .cov coverate file
  df.YMERS = read.csv(COVFILE, sep='\t', row.names = 1)
  # yscore: log2 ratio of value to average (weighted mean) so >0 means more Ymers than average, 1=twice average
  df.YMERS$Yscore_ymers = log2(
    (df.YMERS$meandepth + 1e-4) / # epsilon 1e-4 keeps depth=0 from separating too much from depth=very small
      (weighted.mean(df.YMERS$meandepth, df.YMERS$endpos) + 1e-4) 
  )
  return(df.YMERS)
}


# COMBINING YSCORES AND REPHASING

# function to write DF.REPHASE, row per qry contig 
# with hap and phase columns for original phase (h1/h2) and rephase (X/Y) based on Yscore vs THRESHOLD
make_df.rephase = function(DF.ALIGN, # df with alignment Yscores (percent absent, qry contig splitting, combined)
                           DF.DEPTH, # df with M:F depth bias Yscore
                           DF.VRNTS, # df with site fr all-alt-allele Yscore
                           DF.YMERS, # df with ymer coverage Yscore
                           DF.CTGS, # df with qry ctg rows "SeqName" & "SeqLength", eg l.coords$breaksQry (filter out scaffolds)
                           MIN_CTG_LEN = 0, THRESHOLD=0){
  df.REPHASE = DF.CTGS %>%
    filter(SeqLength >= MIN_CTG_LEN) %>%
    select(SeqName) %>%
    mutate(
      # hap as first 2 characters of seq name
      hap = substr(SeqName, 1, 2),
      seqlength = DF.ALIGN[SeqName,]$seqlength,
      # alignment-based components (unscaled)
      aln.Yscore_PctAbsnt = DF.ALIGN[SeqName,]$Yscore_PctAbsnt,
      aln.Yscore_logQsplt = DF.ALIGN[SeqName,]$Yscore_logQsplt ,
      # "raw" values have unscaled variance
      # align yscores combined (var=1, sum)
      raw.Yscore_align = DF.ALIGN[SeqName,]$Yscore_align,
      raw.Yscore_depth = DF.DEPTH[SeqName,]$Yscore_depth,
      raw.Yscore_vrnts = DF.VRNTS[SeqName,]$Yscore_vrnts,
      raw.Yscore_ymers = DF.YMERS[SeqName,]$Yscore_ymers,
      # scale so all variance = 1
      Yscore_align = raw.Yscore_align/sqrt(var(raw.Yscore_align, na.rm = T)),
      Yscore_depth = raw.Yscore_depth/sqrt(var(raw.Yscore_depth, na.rm = T)),
      Yscore_vrnts = raw.Yscore_vrnts/sqrt(var(raw.Yscore_vrnts, na.rm = T)),
      Yscore_ymers = raw.Yscore_ymers/sqrt(var(raw.Yscore_ymers, na.rm = T)),
      Yscore_all = rowSums(across(c(Yscore_align, Yscore_depth, Yscore_vrnts, Yscore_ymers)), na.rm = T),
      # set phase based on Yscore vs THRESHOLD
      phase = if_else(Yscore_all > THRESHOLD, "Y", "X")
    )
  return(df.REPHASE)
}



# PLOTTING

# Add Yscores to coordinates
make_df.coords_Yscores = function(L.COORDS, DF.REPHASE, 
                                  REF_SEQS = NULL, # list of reference sequences to filter for before calling GetScaffoldCooordsDF
                                  MAX_REF_LENGTH = 10000000) {
  if(length(REF_SEQS) > 0){
    L.COORDS[["coords"]] <- L.COORDS$coords %>% filter(ref_SeqName %in% REF_SEQS)
    MAX_REF_LENGTH = max(L.COORDS$breaksRef[REF_SEQS,]$SeqLength) + 1
  }
  
  df.COORDS_YSCORES = GetScaffoldCoordsDF(L.COORDS, MAX_REF_LENGTH =  MAX_REF_LENGTH) %>%
    mutate(
      hap = DF.REPHASE[qry_SeqName,]$hap,
      phase = DF.REPHASE[qry_SeqName,]$phase,
      # component alignment scores
      aln.Yscore_PctAbsnt = DF.REPHASE[qry_SeqName,]$aln.Yscore_PctAbsnt,
      aln.Yscore_logQsplt = DF.REPHASE[qry_SeqName,]$aln.Yscore_logQsplt,
      # unscaled variance Yscores
      raw.Yscore_align = DF.REPHASE[qry_SeqName,]$raw.Yscore_align,
      raw.Yscore_depth = DF.REPHASE[qry_SeqName,]$raw.Yscore_depth,
      raw.Yscore_vrnts = DF.REPHASE[qry_SeqName,]$raw.Yscore_vrnts,
      raw.Yscore_ymers = DF.REPHASE[qry_SeqName,]$raw.Yscore_ymers,
      # variance=1 scaled Yscores
      Yscore_align = DF.REPHASE[qry_SeqName,]$Yscore_align,
      Yscore_depth = DF.REPHASE[qry_SeqName,]$Yscore_depth,
      Yscore_vrnts = DF.REPHASE[qry_SeqName,]$Yscore_vrnts,
      Yscore_ymers = DF.REPHASE[qry_SeqName,]$Yscore_ymers,
      # combined score: var=1 and sum
      Yscore_all = DF.REPHASE[qry_SeqName,]$Yscore_all,
      # use combined by default, can mutate "Yscore" to above values on function call
      Yscore = Yscore_all)
  return(df.COORDS_YSCORES)
}

# plot alignment colored by Yscore, faceting by phase, alpha by qry_length
dotplot_yscore = function(DF.COORDS_YSCORE_PHASE,
                          REF_LAB = "female primary scaffold"){
  # split point between haps in qry coords
  HAP_SPLIT = DF.COORDS_YSCORE_PHASE %>%
    group_by(hap) %>% 
    summarize(Q1min = min(Q1), 
              Q2min = min(Q2), 
              Q1max = max(Q1), 
              Q2max = max(Q2)) %>% 
    select(-hap) %>% 
    unlist %>% 
    median()
  ggplot(DF.COORDS_YSCORE_PHASE, aes(
      x = R1, 
      y = Q1, 
      color = Yscore,
      alpha = qry_length
    )) +
    geom_point() + 
    geom_point(aes(x = R2, y = Q2)) + 
    geom_segment(aes(xend = R2, yend = Q2)) +
    scale_color_gradientn(colours =  c("#B2182B", "#FDB863", "#2727F5"),
                          values = scales::rescale(c(-4, 0, 4)),
                          limits = c(-4, 4),
                          oob = scales::squish) + 
    geom_hline(yintercept = HAP_SPLIT, linetype = 3) +
    facet_wrap(vars(phase), ncol = 1) +
    labs(
      y = "male contigs, hap1 and hap2",
      x = REF_LAB,
      alpha = "align bp"
    )
}


# ADJUSTING PHASES BASED ON PURGE_DUPS
# function to read purge_dups bed file output and link Yscores for haplotigs pairs
# takes 
#  list of BED filenames from X and Y contigs (names X and Y)
#  df.REPHASE 
# returns list of dataframes (tibble format) with clustered haplotig pairs and corresponding Yscores
make_l.df.purge_dups_pairs = function(l.BEDFILES, df.REPHASE){
  # contig seqnames, use for filtering purge_dups output
  CTGS = df.REPHASE %>%
    select(SeqName) %>%
    unique %>% unlist %>% as.character()
  # bed file output from purge_dups on MhX contigs
  DF.PURGE_DUPS_X = read.table(
    l.BEDFILES$X, col.names = c("ctgB", "start", "stop", "type", "ctgA"), fill = T ) %>% 
    # filter for X contigs and non-OVLP lines (overlaps have ambiguous clustering)
    filter(ctgA %in% CTGS & ctgB %in% CTGS & type != "OVLP") %>%
    # once filtered for contigs/haplotigs all ctgA lines are distinct (these are kept by purge_dups)
    group_by(ctgA) %>%
    summarize(
      # Yscore for ctgA
      Y_A = df.REPHASE[ctgA,]$Yscore_all[1], 
      # contig A length
      A_length = df.REPHASE[ctgA,]$seqlength[1], 
      # list of all ctgB paired with ctgA
      l.ctgB=list(ctgB), 
      # total length of ctgB set
      B_length = sum(df.REPHASE[ctgB,]$seqlength), 
      # average Yscore for ctgB set (mean weighted by length)
      Y_B = weighted.mean(df.REPHASE[ctgB,]$Yscore_all, df.REPHASE[ctgB,]$seqlength),
      # difference between lengths
      d = abs(A_length - B_length)/mean(c(A_length, B_length)) )
  # bed file output from purge_dups on MhY contigs
  DF.PURGE_DUPS_Y = read.table(
    l.BEDFILES$Y, col.names = c("ctgB", "start", "stop", "type", "ctgA"), fill = T ) %>% 
    # filter for Y contigs and HAPLOTIG lines (otherwise ambiguous clustering)
    filter(ctgA %in% CTGS & ctgB %in% CTGS & type != "OVLP") %>%
    # once filtered for contigs/haplotigs all ctgA lines are distinct (these are kept by purge_dups)
    group_by(ctgA) %>%
    summarize(
      # Yscore for ctgA
      Y_A = df.REPHASE[ctgA,]$Yscore_all[1], 
      # contig A length
      A_length = df.REPHASE[ctgA,]$seqlength[1], 
      # list of all ctgB paired with ctgA
      l.ctgB=list(ctgB), 
      # total length of ctgB set
      B_length = sum(df.REPHASE[ctgB,]$seqlength), 
      # average Yscore for ctgB set (mean weighted by length)
      Y_B = weighted.mean(df.REPHASE[ctgB,]$Yscore_all, df.REPHASE[ctgB,]$seqlength),
      # difference between lengths
      d = abs(A_length - B_length)/mean(c(A_length, B_length)))
 l.PURGE_DUPS = list(X = DF.PURGE_DUPS_X, Y = DF.PURGE_DUPS_Y)
 return(l.PURGE_DUPS)
}

# function to update phasing calls with purge_dups input
# set a threshold for difference in length between haplotigs, avoids small ctg forcing rephasing of much bigger pair
# set a threshold for Y score of "Y to X" changes, since high Yscore not likely in true X sequence
update_phasing = function(df.REPHASE, l.PURGE_DUPS,
                          LENGTH_DIFF_THRESHOLD = 0.5,
                          X_YSCORE_THRESHOLD = 2){
  # save original phasing
  df.REPHASE$phase_v0 = df.REPHASE$phase
  # contigs to swap from X to Y
  v.CTG_newY = c(
    (l.PURGE_DUPS$X %>%
      # filter by Y score and length difference 
      filter(Y_A > Y_B,
             d < LENGTH_DIFF_THRESHOLD) %>%
      select(ctgA) %>%
      unlist() %>% as.character()) ,
    (l.PURGE_DUPS$X %>% 
       filter(Y_A <= Y_B,
              d < LENGTH_DIFF_THRESHOLD) %>%
       select(l.ctgB) %>%
       unlist() %>% as.character())
  )
  # contigs to swap from Y to X
  v.CTG_newX = c(
    (l.PURGE_DUPS$Y %>% 
       filter(Y_A <= Y_B,
              d < LENGTH_DIFF_THRESHOLD,
              Y_A < X_YSCORE_THRESHOLD) %>%
       select(ctgA) %>%
       unlist() %>% as.character()) ,
    (l.PURGE_DUPS$Y %>%
       filter(Y_A > Y_B,
              d < LENGTH_DIFF_THRESHOLD,
              Y_A < X_YSCORE_THRESHOLD) %>%
       select(l.ctgB) %>%
       unlist() %>% as.character())
  )
  # reassign phasing column
  if(length(v.CTG_newY) > 0){
    df.REPHASE[v.CTG_newY,]$phase = "Y"
  }
  if(length(v.CTG_newX) > 0){
  df.REPHASE[v.CTG_newX,]$phase = "X"
  }
  # return modified df
  return(df.REPHASE)
}


