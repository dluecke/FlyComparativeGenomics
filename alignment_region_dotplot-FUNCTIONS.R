# alignment_region_dotplot-FUNCTIONS.R 
# takes alignment coordinates as formatted by convert_gnuplot_to_tsv.sh
# produces dotplot of alignment regions 
# FUNCTIONS file to source by other scripts

library(RColorBrewer)
library(plotly)
library(tidyverse)

# FUNCTIONS

# function to get match coordinate and scaffold ref/qry breaks dataframes, in list format
make_l.coords_breaks <- function(INFILE, INFILE_BREAKS){
  # Make df for coordinates of blast hits
  df.coords <- read.csv(INFILE, sep = "\t")
  df.coords <- df.coords[order(df.coords$R1),]
  # use the match midpoint for sequence ID assignment to avoid boundary issues
  df.coords$Rmid <- rowMeans( subset(df.coords, select = c(R1, R2)) )
  df.coords$Qmid <- rowMeans( subset(df.coords, select = c(Q1, Q2)) )
  # get lengths of hits in ref and qry space, and match orientation
  df.coords$ref_length <- df.coords$R2 - df.coords$R1
  df.coords$qry_length <- df.coords$Q2 - df.coords$Q1
  df.coords$orientation <- apply(df.coords, 1, FUN = function(x){
    if( sign(x[7]) == sign(x[8]) ){return("forward")}
    if( sign(x[7]) != sign(x[8]) ){return("reverse")}
  })
  
  # Make df for scaffold positions in reference and query alignment sequences
  df.breaks <- read.csv(INFILE_BREAKS, sep = "\t")
  # using code line from .gp file to find break between Reference and Query scaffolds (gap between lines)
  # find the row with a gap in the code line, will be last line of Reference rows
  breaks.gpLineJumpRow <- which(diff(df.breaks$gpLine)>1)
  df.breaksRef <- df.breaks[1:breaks.gpLineJumpRow,]
  df.breaksQry <- df.breaks[(breaks.gpLineJumpRow+1):nrow(df.breaks),]
  # assign rownames for later value access
  rownames(df.breaksRef) <- df.breaksRef$SeqName
  rownames(df.breaksQry) <- df.breaksQry$SeqName
  # add length for each sequence
  df.breaksRef$SeqLength <- c( diff(df.breaksRef$Position) + 1, NaN )
  df.breaksQry$SeqLength <- c( diff(df.breaksQry$Position) + 1, NaN )
  
  # Assign Ref and Qry scaffold IDs for each alignment line in df.coords based on position in total alignment
  df.coords$ref_SeqName <- sapply(df.coords$Rmid, function(x){ 
    df.breaksRef$SeqName[ which( df.breaksRef$Position > x )[1] - 1]
  })
  df.coords$qry_SeqName <- sapply(df.coords$Qmid, function(x){ 
    df.breaksQry$SeqName[ which( df.breaksQry$Position > x )[1] - 1]
  })
  
  # drop last row, no Scaffold name (need to keep for previous step assigning ref/qry scaf IDs)
  df.breaksRef <- df.breaksRef[-nrow(df.breaksRef),]
  df.breaksQry <- df.breaksQry[-nrow(df.breaksQry),]
  
  return( list(coords = df.coords, breaksRef = df.breaksRef, breaksQry = df.breaksQry, filename = INFILE) )
}

# Function to build a df with subset of coordinates from given sequence names in order given, by default returns df.coords
GetScaffoldCoordsDF <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, 
                                DROP_EMPTY_SCAFFOLDS=T,
                                MIN_REF_LENGTH=0, MIN_QRY_LENGTH=0,
                                MIN_MATCH_LENGTH=0){
  
  COORDS <- L.COORDS$coords
  BREAKS_REF <- L.COORDS$breaksRef
  BREAKS_QRY <- L.COORDS$breaksQry
  
  # Set default REFERENCE and QUERY based on L.COORDS breaks df 
  if( is.null(REFERENCE) ){
    REFERENCE <- BREAKS_REF$SeqName
  }
  if( is.null(QUERY) ){
    QUERY <- BREAKS_QRY$SeqName
  }
  
  # filter qry and ref by length
  REFERENCE <- REFERENCE[BREAKS_REF[REFERENCE,]$SeqLength > MIN_REF_LENGTH]
  QUERY <- QUERY[BREAKS_QRY[QUERY,]$SeqLength > MIN_QRY_LENGTH]
  
  # get matches for qry and reference
  DF.COORDS_SUBSET <- COORDS[ (COORDS$ref_SeqName %in% REFERENCE) & (COORDS$qry_SeqName %in% QUERY) ,]
  
  # filter based on minimum match length
  DF.COORDS_SUBSET <- DF.COORDS_SUBSET[ (abs(DF.COORDS_SUBSET$ref_length) > MIN_MATCH_LENGTH) &
                                          (abs(DF.COORDS_SUBSET$qry_length) > MIN_MATCH_LENGTH),]
  
  if( DROP_EMPTY_SCAFFOLDS ){
    REFERENCE <- REFERENCE[ REFERENCE %in% DF.COORDS_SUBSET$ref_SeqName ]
    QUERY <- QUERY[ QUERY %in% DF.COORDS_SUBSET$qry_SeqName ]
  }
  
  # lengths of remaining qry/scaffs
  LENGTHS_REF <- sapply(REFERENCE, function(x) { as.numeric(BREAKS_REF[x,]$SeqLength) }) 
  LENGTHS_QRY <- sapply(QUERY, function(x) { as.numeric(BREAKS_QRY[x,]$SeqLength) })
  
  STARTPOS_REF <- c(1, 1+cumsum(LENGTHS_REF)[-length(LENGTHS_REF)] )
  names(STARTPOS_REF) <- REFERENCE
  STARTPOS_QRY <- c(1, 1+cumsum(LENGTHS_QRY)[-length(LENGTHS_QRY)] )
  names(STARTPOS_QRY) <- QUERY
  
  
  DF.COORDS_SUBSET$ref_StartPos <- sapply(DF.COORDS_SUBSET$ref_SeqName, function(x) STARTPOS_REF[x])
  DF.COORDS_SUBSET$qry_StartPos <- sapply(DF.COORDS_SUBSET$qry_SeqName, function(x) STARTPOS_QRY[x])
  
  DF.COORDS_SUBSET$ref_SeqLength <- sapply(DF.COORDS_SUBSET$ref_SeqName, function(x) LENGTHS_REF[x])
  DF.COORDS_SUBSET$qry_SeqLength <- sapply(DF.COORDS_SUBSET$qry_SeqName, function(x) LENGTHS_QRY[x])
  
  
  DF.COORDS_SUBSET[,c('R1','R2')] <- DF.COORDS_SUBSET[,c('R1','R2')] - 
    BREAKS_REF[DF.COORDS_SUBSET$ref_SeqName,]$Position + 
    DF.COORDS_SUBSET$ref_StartPos 
  
  DF.COORDS_SUBSET[,c('Q1','Q2')] <- DF.COORDS_SUBSET[,c('Q1','Q2')] - 
    BREAKS_QRY[DF.COORDS_SUBSET$qry_SeqName,]$Position + 
    DF.COORDS_SUBSET$qry_StartPos
  
  return(DF.COORDS_SUBSET)
  
}


PlotDFCoords <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, REORDER_QRY = F,  # reorder query names by ref position
                         MIN_REF_LENGTH = 0, MIN_QRY_LENGTH = 0, MIN_MATCH_LENGTH = 0, 
                         TIC_LABELS = FALSE, TIC_COUNT = 40, 
                         DROP_EMPTY_SCAFFOLDS = T,
                         ALPHA = 0.25, POINTSIZE = 0.8, COORD_OFFSET = 0.02,
                         LAB_REF = "", LAB_QRY = ""){
  
  COORDS <- L.COORDS$coords
  BREAKS_REF <- L.COORDS$breaksRef
  BREAKS_QRY <- L.COORDS$breaksQry
  
  # Set default REFERENCE and QUERY based on L.COORDS breaks df 
  if( is.null(REFERENCE) ){
    REFERENCE <- BREAKS_REF$SeqName
  }
  if( is.null(QUERY) ){
    QUERY <- BREAKS_QRY$SeqName
  }

  # filter qry and ref by length
  REFERENCE <- REFERENCE[BREAKS_REF[REFERENCE,]$SeqLength > MIN_REF_LENGTH]
  QUERY <- QUERY[BREAKS_QRY[QUERY,]$SeqLength > MIN_QRY_LENGTH]
  COORDS <- COORDS[((COORDS$ref_SeqName %in% REFERENCE) & (COORDS$qry_SeqName %in% QUERY)),]
  
  # reorder the query sequences to best show contiguity across contigs
  if( REORDER_QRY ){
    # longest reference match for each scaffold pair
    topR_coords <- COORDS %>% dplyr::group_by(ref_SeqName, qry_SeqName) %>% top_n(1, ref_length)
    topQ_coords <- topR_coords %>% dplyr::group_by(qry_SeqName) %>% top_n(1, ref_length)
    ORDERED_QRY_ALL <- topQ_coords$qry_SeqName %>% unique() # COORDS[order(COORDS$R1),]$qry_SeqName %>% unique()
    QUERY <- ORDERED_QRY_ALL[ ORDERED_QRY_ALL %in% QUERY ]
  }
  
  DF_COORDS <- GetScaffoldCoordsDF(L.COORDS = L.COORDS, REFERENCE = REFERENCE, QUERY = QUERY, 
                                   DROP_EMPTY_SCAFFOLDS = DROP_EMPTY_SCAFFOLDS,
                                   MIN_REF_LENGTH = MIN_REF_LENGTH, MIN_QRY_LENGTH = MIN_QRY_LENGTH,
                                   MIN_MATCH_LENGTH = MIN_MATCH_LENGTH)
  
  # By default only plot scaffolds with matches
  if( DROP_EMPTY_SCAFFOLDS ){  
    
    REF_STARTPOS <- unique(DF_COORDS$ref_StartPos)
    REF_SEQNAMES <- unique(DF_COORDS$ref_SeqName)
    QRY_STARTPOS <- unique(DF_COORDS$qry_StartPos)
    QRY_SEQNAMES <- unique(DF_COORDS$qry_SeqName)
    REF_MAX = max(DF_COORDS$ref_StartPos + DF_COORDS$ref_SeqLength)
    QRY_MAX = max(DF_COORDS$qry_StartPos + DF_COORDS$qry_SeqLength)
    
  } else {
    # otherwise need to replicate some of GetScaffoldCoordsDF() but without filtering for matches
    DF_BREAKSREF <- L.COORDS$breaksRef[REFERENCE,]
    DF_BREAKSQRY <- L.COORDS$breaksQry[QUERY,]
    
    REF_STARTPOS <- c(1, 1+cumsum(DF_BREAKSREF$SeqLength)[-length(REFERENCE)])
    REF_SEQNAMES <- REFERENCE
    QRY_STARTPOS <- c(1, 1+cumsum(DF_BREAKSQRY$SeqLength)[-length(QUERY)])
    QRY_SEQNAMES <- QUERY
    REF_MAX <- DF_BREAKSREF$Position[nrow(DF_BREAKSREF)] + DF_BREAKSREF$SeqLength[nrow(DF_BREAKSREF)]
    QRY_MAX <- DF_BREAKSQRY$Position[nrow(DF_BREAKSQRY)] + DF_BREAKSQRY$SeqLength[nrow(DF_BREAKSQRY)]
    
  }
  
  # find subset of scaffold positions for axis marking, if requested
  if( TIC_COUNT > 0 ){
    REF_TIC_SPACING = as.integer( REF_MAX / TIC_COUNT )
    QRY_TIC_SPACING = as.integer( QRY_MAX / TIC_COUNT )
    
    GetAxisTics <- function(POSITIONS, SPACING) {
      TICGUIDE <- c( seq( min(POSITIONS), max(POSITIONS), by = SPACING ), max(POSITIONS) )
      TICINDEX <- sapply(TICGUIDE, function(guidepoint) which.min( abs(POSITIONS - guidepoint) ) ) %>% unique()
      TICPOSITIONS <- POSITIONS[TICINDEX][unique( c(1, which( diff(POSITIONS[TICINDEX[-length(TICINDEX)]]) > SPACING ), length(TICINDEX)) )] 
      return(TICPOSITIONS[ diff(TICPOSITIONS) > SPACING ])
    }
    
    l.REF_TIC_POSITIONS <- lapply(REF_SEQNAMES, function(seq){ 
      GetAxisTics( sort(c( DF_COORDS[DF_COORDS$ref_SeqName==seq,]$R1, DF_COORDS[DF_COORDS$ref_SeqName==seq,]$R2 )) , REF_TIC_SPACING)
    })
    names(l.REF_TIC_POSITIONS) <- REF_SEQNAMES
    l.REF_TIC_LABELS <- lapply(REF_SEQNAMES, function(seq){ 
      l.REF_TIC_POSITIONS[[seq]] - DF_COORDS[DF_COORDS$ref_SeqName==seq,]$ref_StartPos[1] + 1
    })
    names(l.REF_TIC_LABELS) <- REF_SEQNAMES
    df.REF_TICS <- data.frame(position = unlist(l.REF_TIC_POSITIONS), label = as.character(unlist(l.REF_TIC_LABELS)))
    
    l.QRY_TIC_POSITIONS <- lapply(QRY_SEQNAMES, function(seq){ 
      GetAxisTics(sort(c(DF_COORDS[DF_COORDS$qry_SeqName==seq,]$Q1, DF_COORDS[DF_COORDS$qry_SeqName==seq,]$Q2)), QRY_TIC_SPACING)
    })
    names(l.QRY_TIC_POSITIONS) <- QRY_SEQNAMES
    l.QRY_TIC_LABELS <- lapply(QRY_SEQNAMES, function(seq){ 
      l.QRY_TIC_POSITIONS[[seq]] - DF_COORDS[DF_COORDS$qry_SeqName==seq,]$qry_StartPos[1] + 1
    })
    names(l.QRY_TIC_LABELS) <- QRY_SEQNAMES
    df.QRY_TICS <- data.frame(position = unlist(l.QRY_TIC_POSITIONS), label = as.character(unlist(l.QRY_TIC_LABELS)))
  } else {
    df.REF_TICS <- data.frame(position = c(REF_MAX))
    df.QRY_TICS <- data.frame(position = c(QRY_MAX))
  }
  
  p <- ggplot(DF_COORDS, aes(x=R1, y=Q1))
  p2 <- p + 
    geom_segment(aes(xend=R2, yend=Q2, color=orientation)) + 
    geom_point(aes(color=orientation), alpha=ALPHA, size=POINTSIZE) +
    geom_point(aes(x=R2, y=Q2, color=orientation), alpha=ALPHA, size=POINTSIZE) +
    geom_segment(data = data.frame(x = c(REF_STARTPOS, REF_MAX), xend=c(REF_STARTPOS, REF_MAX), y=0, yend=QRY_MAX), 
                 aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', linewidth=0.25) + 
    geom_segment(data = data.frame(x=0, xend=REF_MAX, y = c(QRY_STARTPOS, QRY_MAX), yend = c(QRY_STARTPOS, QRY_MAX)), 
                 aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', linewidth=0.25) +
    labs(x=LAB_REF, y=LAB_QRY) +
    scale_x_continuous(breaks = c(REF_STARTPOS, df.REF_TICS$position),
                       labels = c(REF_SEQNAMES, rep("", nrow(df.REF_TICS))),
                       minor_breaks = NULL) + #, expand = c(0.05,0), limits = c(-REF_MAX*0.1,REF_MAX*1.1)) +
    scale_y_continuous(breaks = c(QRY_STARTPOS, df.QRY_TICS$position),
                       labels = c(QRY_SEQNAMES, rep("", nrow(df.QRY_TICS))), position = "right" ,
                       minor_breaks = NULL ) + #, expand = c(0.05,0), limits = c(-QRY_MAX*1.1, QRY_MAX)) +
    # just mapping orientation to color via aes and scale_color_brewer means only-reverse plots use the otherwise "forward" color
    # so set manually to keep consistent, but still use the nice palette
    scale_color_manual(values = c("forward" = brewer.pal(3,"Dark2")[1],
                                  "reverse" = brewer.pal(3,"Dark2")[2])) + 
    theme_minimal() +
    theme(plot.margin = margin(10,0,5,10), axis.title.y.right = element_text(angle = 90),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0, color = "slateblue", face = "bold"), 
          axis.text.y = element_text(angle = 45, hjust = 0, vjust = 0, color = "slateblue", face = "bold") )
  
  if( TIC_LABELS ){
    p2 <- p2 +
      geom_text(data = df.REF_TICS, aes(x = position, y = -1, label = label), angle = -45, hjust = 0, vjust = 1, size = 2) +
      geom_text(data = df.QRY_TICS, aes(x = REF_MAX + 1, y = position, label = label), angle = 45, hjust = 0, vjust = 1, size = 2) +
      coord_cartesian(xlim = c(COORD_OFFSET*REF_MAX, (1+COORD_OFFSET)*REF_MAX), ylim = c(-COORD_OFFSET*QRY_MAX, (1+COORD_OFFSET)*QRY_MAX)) 
  } else {
    p2 <- p2 +
      coord_cartesian(xlim = c(COORD_OFFSET*REF_MAX, REF_MAX*(1-COORD_OFFSET)), ylim = c(COORD_OFFSET*QRY_MAX, (1-COORD_OFFSET)*QRY_MAX)) 
  }
  
  return(p2)
}


# Function to find which reference and query sequences have lowest match percent
# takes l.coords, returns list with dataframes of ref and qry seqs sorted by match coverage
GetMatchStats <- function(L.COORDS){
  
  # DF for reference sequence match statistics
  # get stats by reference sequence
  REF_MATCH_STATS <- L.COORDS$coords  %>% group_by(ref_SeqName) %>% 
    summarise( total_match = sum(ref_length), longest_match = max(ref_length), n_matches = length(ref_length) )
  # ref sequences with no matches (in breaks but not coords)
  REF_NO_MATCH <- rownames(L.COORDS$breaksRef)[ ! (rownames(L.COORDS$breaksRef) %in% REF_MATCH_STATS$ref_SeqName)]
  # add to ref match stats df with match stats = 0
  REF_MATCH_STATS <- rbind(REF_MATCH_STATS, 
                        data.frame(ref_SeqName = REF_NO_MATCH,
                                   total_match = 0, longest_match = 0, n_matches = 0))
  # get lengths for ref seqs
  REF_MATCH_STATS$SeqLength <- L.COORDS$breaksRef[REF_MATCH_STATS$ref_SeqName,]$SeqLength
  # calculate length-normalized stats
  REF_MATCH_STATS$match_pct <- REF_MATCH_STATS$total_match / REF_MATCH_STATS$SeqLength
  REF_MATCH_STATS$longest_pct <- REF_MATCH_STATS$longest_match / REF_MATCH_STATS$SeqLength
  # sort by percent sequence with match
  REF_MATCH_STATS <- REF_MATCH_STATS[order(REF_MATCH_STATS$match_pct),]
  
  # DF for query sequence match statistics, same process
  QRY_MATCH_STATS <- L.COORDS$coords  %>% group_by(qry_SeqName) %>% 
    summarise( total_match = sum(ref_length), longest_match = max(ref_length), n_matches = length(ref_length) )
  QRY_NO_MATCH <- rownames(L.COORDS$breaksQry)[ ! (rownames(L.COORDS$breaksQry) %in% QRY_MATCH_STATS$qry_SeqName)] 
  QRY_MATCH_STATS <- rbind(QRY_MATCH_STATS, 
                           data.frame(qry_SeqName = QRY_NO_MATCH,
                                      total_match = 0, longest_match = 0, n_matches = 0))
  QRY_MATCH_STATS$SeqLength <- L.COORDS$breaksQry[QRY_MATCH_STATS$qry_SeqName,]$SeqLength
  QRY_MATCH_STATS$match_pct <- QRY_MATCH_STATS$total_match / QRY_MATCH_STATS$SeqLength
  QRY_MATCH_STATS$longest_pct <- QRY_MATCH_STATS$longest_match / QRY_MATCH_STATS$SeqLength
  QRY_MATCH_STATS <- QRY_MATCH_STATS[order(QRY_MATCH_STATS$match_pct),]
  

  return(list( RefStats = REF_MATCH_STATS,
               QryStats = QRY_MATCH_STATS ))
}


# Function to get scaffold clusters of mutual homology based on top pairwise sum of match lengths
# written to be embedded in GetHomologs() below
# takes main reference and query scaffold as definted in GetHomologs() df.HOMOLOGS
# takes full pairwise list of match lengths sums as GetHomologs() df.MATCHES
GetHomologCluster <- function(MainRef, MainQry, DF.MATCHES){
  
  MatchedRefs <- DF.MATCHES %>% 
    filter(qry_SeqName == MainQry, matchsum > 0) %>% select(ref_SeqName)
  ClusterRefs <- DF.MATCHES %>% 
    filter(ref_SeqName %in% MatchedRefs$ref_SeqName) %>% 
    group_by(ref_SeqName) %>% slice(which.max(matchsum)) %>% 
    filter(qry_SeqName == MainQry) %>% select(ref_SeqName) %>% 
    arrange(nchar(ref_SeqName))
  
  MatchedQrys <- DF.MATCHES %>% 
    filter(ref_SeqName == MainRef, matchsum > 0) %>% select(qry_SeqName)
  ClusterQrys <- DF.MATCHES %>% 
    filter(qry_SeqName %in% MatchedQrys$qry_SeqName) %>% 
    group_by(qry_SeqName) %>% slice(which.max(matchsum)) %>% 
    filter(ref_SeqName == MainRef) %>% select(qry_SeqName) %>% 
    arrange(nchar(qry_SeqName))
  
  return(list(RefScafs = ClusterRefs$ref_SeqName,
              QryScafs = ClusterQrys$qry_SeqName))
}

# Function to get unambiguous homologous scaffolds across female/male assemblies
# takes L.COORDS from make_l.coords_breaks()
# outputs list with unambiguous list and all matches sorted by sum of match lengths
GetHomologs <- function(L.COORDS){
  
  # sum of match lengths for each pair
  df.MATCHES <- L.COORDS$coords %>% group_by(ref_SeqName, qry_SeqName) %>% 
    summarize(matchsum=sum(ref_length)) %>% arrange(desc(matchsum)) %>% ungroup
  
  # top partner for each scaffold in both assemblies
  df.HOMOLOGS <- rbind(df.MATCHES[!duplicated(df.MATCHES$ref_SeqName),], 
                       df.MATCHES[!duplicated(df.MATCHES$qry_SeqName),]) %>% 
    # keep reciprocal pairs (present in both lists)
    group_by(ref_SeqName, qry_SeqName) %>% count %>% filter(n == 2) %>% 
    # keep scaffold name columns and order by name length
    select(ref_SeqName, qry_SeqName) %>% arrange(nchar(ref_SeqName)) %>% ungroup
  
  # find homology clusters based on unambiguous pairs in df.HOMOLOGS
  HomologsIndex <- paste("Homolog", sprintf("%03d", seq_len(nrow(df.HOMOLOGS))), sep = "_")
  l.CLUSTERS <- by(df.HOMOLOGS, HomologsIndex , function(x) 
    GetHomologCluster(x$ref_SeqName, x$qry_SeqName, df.MATCHES) )
  # check that clusters are non-overlapping
  AllClusterRefs <- lapply(l.CLUSTERS, function(x){ x$RefScafs }) %>% unlist
  if( sum(duplicated(AllClusterRefs) ) > 0) {
    warning( paste0(c("Homology clusters have overlapping reference scaffolds:",
                     AllClusterRefs[duplicated(AllClusterRefs)]), 
                   sep = " " ) )
  }
  AllClusterQrys <- lapply(l.CLUSTERS, function(x){ x$QryScafs }) %>% unlist
  if( sum(duplicated(AllClusterQrys) ) > 0) {
    warning( paste0(c("Homology clusters have overlapping query scaffolds:",
                      AllClusterQrys[duplicated(AllClusterQrys)]), 
                    sep = " " ) )
  }
  
  return(list(homologs = df.HOMOLOGS,
              matches = df.MATCHES,
              clusters = l.CLUSTERS))
}

