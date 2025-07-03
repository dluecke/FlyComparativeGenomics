# alignment_region_dotplot-FUNCTIONS.R 
# takes alignment coordinates as formatted by convert_gnuplot_to_tsv.sh
# produces dotplot of alignment regions 
# FUNCTIONS file to source by other scripts

library(RColorBrewer)
#library(plotly)
library(tidyverse)

# FUNCTIONS

# new version using mummer's coord output format, along with fai index files for seq lengths
#   FAI inputs need sequence ID and length columns, 
#   can use/reuse any version so long as all aligned sequences represented
# function to list dfs of match coordinates and sequence ref/qry lengths ("breaks")
# copies format of older version l.coord header names etc
make_l.coords <- function(COORDSFILE, REF_FAI, QRY_FAI){
  # read coords format, skip first 5 header lines
  df.coords <- read.table(COORDSFILE, stringsAsFactors = F, skip = 5)
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
  
  # "breaks" dfs for tracking sequence lengths
  df.breaksRef <- read.table(REF_FAI)
  df.breaksRef <- df.breaksRef[,c(1,2)] 
  colnames(df.breaksRef) <- c("SeqName", "SeqLength")
  rownames(df.breaksRef) <- df.breaksRef$SeqName
  
  df.breaksQry <- read.table(QRY_FAI)
  df.breaksQry <- df.breaksQry[,c(1,2)]
  colnames(df.breaksQry) <- c("SeqName", "SeqLength")
  rownames(df.breaksQry) <- df.breaksQry$SeqName
  
  return( list(coords = df.coords, breaksRef = df.breaksRef, breaksQry = df.breaksQry, filename = COORDSFILE) )
}

# old version which takes plot.tsv and breaks.tsv custom outputs produced from mummer .delta and .gp
# function to get match coordinate and scaffold ref/qry breaks dataframes, in list format
make_l.coords_breaks <- function(INFILE, INFILE_BREAKS){
  # Make df for coordinates of blast hits
  df.coords <- read.csv(INFILE, sep = "\t")
  # not necessary for plotting but makes visually scanning df easier
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
  
  # convert coordinates from alignment-space to sequence-space (shift back to aligment position when plotting)
  df.coords[,c('R1','R2')] <- df.coords[,c('R1','R2')] - 
    df.breaksRef[df.coords$ref_SeqName,]$Position
  
  df.coords[,c('Q1','Q2')] <- df.coords[,c('Q1','Q2')] - 
    df.breaksQry[df.coords$qry_SeqName,]$Position
  
  return( list(coords = df.coords, breaksRef = df.breaksRef, breaksQry = df.breaksQry, filename = INFILE) )
}

# Function to build a df with subset of coordinates from given sequence names in order given, by default returns df.coords
GetScaffoldCoordsDF <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, 
                                DROP_EMPTY_SCAFFOLDS=T,
                                MIN_REF_LENGTH=0, MIN_QRY_LENGTH=0,
                                MAX_REF_LENGTH=F, MAX_QRY_LENGTH=F,
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
  # if filtering for MAX length (by default skip)
  if( is.numeric(MAX_REF_LENGTH) ){
    REFERENCE <- REFERENCE[BREAKS_REF[REFERENCE,]$SeqLength < MAX_REF_LENGTH]
  }
  if( is.numeric(MAX_QRY_LENGTH) ){
    QUERY <- QUERY[BREAKS_QRY[QUERY,]$SeqLength < MAX_QRY_LENGTH]
  }
  
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
  
  DF.COORDS_SUBSET$ref_StartPos <- STARTPOS_REF[DF.COORDS_SUBSET$ref_SeqName]
  DF.COORDS_SUBSET$qry_StartPos <- STARTPOS_QRY[DF.COORDS_SUBSET$qry_SeqName]
  
  DF.COORDS_SUBSET$ref_SeqLength <- LENGTHS_REF[DF.COORDS_SUBSET$ref_SeqName]
  DF.COORDS_SUBSET$qry_SeqLength <- LENGTHS_QRY[DF.COORDS_SUBSET$qry_SeqName]
  
  # as originally written, needlessly complicated
  # DF.COORDS_SUBSET$ref_StartPos <- sapply(DF.COORDS_SUBSET$ref_SeqName, function(x) STARTPOS_REF[x])
  # DF.COORDS_SUBSET$qry_StartPos <- sapply(DF.COORDS_SUBSET$qry_SeqName, function(x) STARTPOS_QRY[x])
  # 
  # DF.COORDS_SUBSET$ref_SeqLength <- sapply(DF.COORDS_SUBSET$ref_SeqName, function(x) LENGTHS_REF[x])
  # DF.COORDS_SUBSET$qry_SeqLength <- sapply(DF.COORDS_SUBSET$qry_SeqName, function(x) LENGTHS_QRY[x])
  
  DF.COORDS_SUBSET[,c('R1','R2')] <- DF.COORDS_SUBSET[,c('R1','R2')] + #- 
#    BREAKS_REF[DF.COORDS_SUBSET$ref_SeqName,]$Position + # now doing this in make_l.coord_breaks()
    DF.COORDS_SUBSET$ref_StartPos 
  
  DF.COORDS_SUBSET[,c('Q1','Q2')] <- DF.COORDS_SUBSET[,c('Q1','Q2')] + #- 
#    BREAKS_QRY[DF.COORDS_SUBSET$qry_SeqName,]$Position + 
    DF.COORDS_SUBSET$qry_StartPos
  
  # Add pctID column if doesn't exist (default 100)
  if( ! "pctID" %in% colnames(DF.COORDS_SUBSET) ){
    DF.COORDS_SUBSET$pctID = 100
  }
  
  return(DF.COORDS_SUBSET)
  
}


PlotDFCoords <- function(L.COORDS, REFERENCE=NULL, QUERY=NULL, REORDER_QRY = F,  # reorder query names by ref position
                         MIN_REF_LENGTH = 0, MIN_QRY_LENGTH = 0, MIN_MATCH_LENGTH = 0, 
                         MAX_REF_LENGTH = F, MAX_QRY_LENGTH = F,
                         REF_LIM = NULL, QRY_LIM = NULL, # akin to xlim, ylim - takes c(START, STOP) coordinates
                         TIC_LABELS = FALSE, TIC_COUNT = 40, 
                         DROP_EMPTY_SCAFFOLDS = T, FLIP_AXES = F, COORD_OFFSET = 0.02,
                         LINEWIDTH = 1, ALPHA = 0.25, POINTSIZE = 0.8,
                         SEQLABANGLE = 45,
                         MAR_T = 25, MAR_R = 5, MAR_B = 5, MAR_L = 10,
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
  # if filtering for MAX length (by default skip)
  if( is.numeric(MAX_REF_LENGTH) ){
    REFERENCE <- REFERENCE[BREAKS_REF[REFERENCE,]$SeqLength < MAX_REF_LENGTH]
  }
  if( is.numeric(MAX_QRY_LENGTH) ){
    QUERY <- QUERY[BREAKS_QRY[QUERY,]$SeqLength < MAX_QRY_LENGTH]
  }
  
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
                                   MAX_REF_LENGTH = MAX_REF_LENGTH, MAX_QRY_LENGTH = MAX_QRY_LENGTH,
                                   MIN_MATCH_LENGTH = MIN_MATCH_LENGTH)
  
  if( nrow(DF_COORDS) > 0 ){  
    
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
      REF_MAX <- REF_STARTPOS[length(REF_STARTPOS)] + DF_BREAKSREF$SeqLength[nrow(DF_BREAKSREF)]
      QRY_MAX <- QRY_STARTPOS[length(QRY_STARTPOS)] + DF_BREAKSQRY$SeqLength[nrow(DF_BREAKSQRY)]
    }
    
    # adjust boundaries for REF_LIM and QRY_LIM
    if( !is.null(REF_LIM) ){
      REF_SEQNAMES = REF_SEQNAMES[intersect(which(c(REF_STARTPOS,REF_MAX) > REF_LIM[1]) - 1, 
                                            which(REF_STARTPOS < REF_LIM[2]))]
      REF_STARTPOS = c(REF_LIM[1], 
                       REF_STARTPOS[(REF_STARTPOS > REF_LIM[1]) & (REF_STARTPOS < REF_LIM[2])])
      REF_MIN = REF_LIM[1]
      REF_MAX = REF_LIM[2]
    } else {
      REF_MIN = 0
    }
    if( !is.null(QRY_LIM) ){
      QRY_SEQNAMES = QRY_SEQNAMES[intersect(which(c(QRY_STARTPOS,QRY_MAX) > QRY_LIM[1]) - 1, 
                                            which(QRY_STARTPOS < QRY_LIM[2]))]
      QRY_STARTPOS = c(QRY_LIM[1], 
                       QRY_STARTPOS[(QRY_STARTPOS > QRY_LIM[1]) & (QRY_STARTPOS < QRY_LIM[2])])
      QRY_MIN = QRY_LIM[1]
      QRY_MAX = QRY_LIM[2]
    } else {
      QRY_MIN = 0
    }
    
    # find subset of scaffold positions for axis marking, if requested (doing this with empty scaffolds caused problems, wouldn't need tics for empties anyway)
    if( TIC_COUNT > 0 & DROP_EMPTY_SCAFFOLDS ){
      REF_TIC_SPACING = as.integer( (REF_MAX-REF_MIN) / TIC_COUNT )
      QRY_TIC_SPACING = as.integer( (QRY_MAX-QRY_MIN) / TIC_COUNT )
      
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
      df.REF_TICS <- data.frame(position = c(REF_MAX), label = c(""))
      df.QRY_TICS <- data.frame(position = c(QRY_MAX), label = c(""))
    }
    
    if( FLIP_AXES ){
      XPOS = "top"
      YPOS = "left"
    } else {
      XPOS = "bottom"
      YPOS = "right"
    }
    
    # only want match % to map to linewidth if pctID is variable (taken from newer .coords input path)
    if(min(DF_COORDS$pctID) < 100){
      p <- ggplot(DF_COORDS, aes(x=R1, y=Q1)) +
        geom_segment(aes(xend=R2, yend=Q2, color = orientation, linewidth = pctID)) +
        scale_linewidth_continuous(range = c(0.4, LINEWIDTH)) +
        labs(x=LAB_REF, y=LAB_QRY, linewidth = "match %")
    } else {
      p <- ggplot(DF_COORDS, aes(x=R1, y=Q1)) +
        geom_segment(aes(xend=R2, yend=Q2, color = orientation), linewidth = LINEWIDTH) +
        labs(x=LAB_REF, y=LAB_QRY)
    }
    
    p2 <- p + 
      geom_point(aes(color=orientation), alpha=ALPHA, size=POINTSIZE) +
      geom_point(aes(x=R2, y=Q2, color=orientation), alpha=ALPHA, size=POINTSIZE) +
      geom_segment(data = data.frame(x = c(REF_STARTPOS, REF_MAX), xend=c(REF_STARTPOS, REF_MAX), y=QRY_MIN, yend=QRY_MAX), 
                   aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', linewidth=0.25) + 
      geom_segment(data = data.frame(x=REF_MIN, xend=REF_MAX, y = c(QRY_STARTPOS, QRY_MAX), yend = c(QRY_STARTPOS, QRY_MAX)), 
                   aes(x=x, xend=xend, y=y, yend=yend), color = "slateblue", linetype='solid', linewidth=0.25) +
      scale_x_continuous(breaks = c(REF_STARTPOS, df.REF_TICS$position), 
                         labels = c(REF_SEQNAMES, rep("", nrow(df.REF_TICS))), position = XPOS,
                         minor_breaks = NULL) +
      scale_y_continuous(breaks = c(QRY_STARTPOS, df.QRY_TICS$position),
                         labels = c(QRY_SEQNAMES, rep("", nrow(df.QRY_TICS))), position = YPOS,
                         minor_breaks = NULL ) + 
      guides(color = guide_legend(override.aes = list(linewidth = LINEWIDTH))) +
      # just mapping orientation to color via aes and scale_color_brewer means only-reverse plots use the otherwise "forward" color
      # so set manually to keep consistent, but still use the nice palette
      scale_color_manual(values = c("forward" = brewer.pal(3,"Dark2")[1],
                                    "reverse" = brewer.pal(3,"Dark2")[2])) + 
      theme_minimal() 
    
    if( TIC_LABELS ){
      if( FLIP_AXES ){
        p2 <- p2 +
          geom_text(data = df.REF_TICS %>% filter(position > REF_MIN, position < REF_MAX), 
                    aes(x = position, y = QRY_MAX + 1, label = label), angle = 45, hjust = 0, vjust = 1, size = 2) +
          geom_text(data = df.QRY_TICS %>% filter(position > QRY_MIN, position < QRY_MAX), 
                    aes(x = REF_MIN-1, y = position, label = label), angle = -45, hjust = 0, vjust = 1, size = 2) +
          coord_flip(xlim = c( (REF_MIN - COORD_OFFSET*(REF_MAX-REF_MIN)), 
                               (REF_MAX + COORD_OFFSET*(REF_MAX-REF_MIN)) ), 
                     ylim = c( (QRY_MIN + COORD_OFFSET*(QRY_MAX-QRY_MIN)), 
                               (QRY_MAX + COORD_OFFSET*(QRY_MAX-QRY_MIN)) ))
      } else {
        p2 <- p2 +
          geom_text(data = df.REF_TICS %>% filter(position > REF_MIN, position < REF_MAX), 
                    aes(x = position, y = QRY_MIN-1, label = label), angle = -45, hjust = 0, vjust = 1, size = 2) +
          geom_text(data = df.QRY_TICS %>% filter(position > QRY_MIN, position < QRY_MAX), 
                    aes(x = REF_MAX + 1, y = position, label = label), angle = 45, hjust = 0, vjust = 1, size = 2) +
          coord_cartesian(xlim = c( (REF_MIN + COORD_OFFSET*(REF_MAX-REF_MIN)), 
                                    (REF_MAX + COORD_OFFSET*(REF_MAX-REF_MIN)) ), 
                          ylim = c( (QRY_MIN - COORD_OFFSET*(QRY_MAX-QRY_MIN)), 
                                    (QRY_MAX + COORD_OFFSET*(QRY_MAX-QRY_MIN)) )) 
      }  
    } else {
      if( FLIP_AXES ){
        p2 <- p2 +
          coord_flip(xlim = c( (REF_MIN + COORD_OFFSET*(REF_MAX-REF_MIN)), 
                               (REF_MAX - COORD_OFFSET*(REF_MAX-REF_MIN)) ), 
                     ylim = c( (QRY_MIN + COORD_OFFSET*(QRY_MAX-QRY_MIN)), 
                               (QRY_MAX - COORD_OFFSET*(QRY_MAX-QRY_MIN)) )) 
      } else {
        p2 <- p2 +
          coord_cartesian(xlim = c( (REF_MIN + COORD_OFFSET*(REF_MAX-REF_MIN)), 
                                    (REF_MAX - COORD_OFFSET*(REF_MAX-REF_MIN)) ), 
                          ylim = c( (QRY_MIN + COORD_OFFSET*(QRY_MAX-QRY_MIN)), 
                                    (QRY_MAX - COORD_OFFSET*(QRY_MAX-QRY_MIN)) )) 
        
      }
    }
    
    return(p2 +
             theme(plot.margin = margin(MAR_T, MAR_R, MAR_B, MAR_L), 
                   axis.title.y.right = element_text(angle = 90),
                   axis.text.x = element_text(angle = -SEQLABANGLE, 
                                              hjust = 0, vjust = 0, 
                                              color = "slateblue", 
                                              face = "bold"), 
                   axis.text.y = element_text(angle = 90-SEQLABANGLE, 
                                              hjust = 0, vjust = 0, 
                                              color = "slateblue", 
                                              face = "bold") )
           
    )
  } else {
    return(NULL)
  }
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
    summarize(matchsum=sum(ref_length)) %>% arrange(desc(matchsum)) %>% ungroup %>%
    mutate(ref_SeqLength = L.COORDS$breaksRef[ref_SeqName,]$SeqLength, 
           qry_SeqLength = L.COORDS$breaksQry[qry_SeqName,]$SeqLength, 
           ref_Prop = matchsum/ref_SeqLength, 
           qry_Prop = matchsum/qry_SeqLength)
  
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


# function to run make_l.coords() and GetHomologs() on all files in IN_FILE list:
# IN_FILE = list(aln1 = c("alignment1_data.coords", "reference1.fa.fai", "query1.fa.fai"),
#                aln2 = c("alignment2_data.coords", "reference2.fa.fai", "query2.fa.fai"))
# returns list (per alignment) of list(l.coord, l.homolog): 
make_l.l_coords.l_homologs = function(L.FILES){
  l.LCoords.LHomologs <- lapply(L.FILES, function(ALN){
    l.coords = make_l.coords(ALN[1], ALN[2], ALN[3])
    l.homologs = GetHomologs(l.coords)
    RETURN = list(l.coords = l.coords,
                  l.homologs = l.homologs)
    return(RETURN)
  })
  return(l.LCoords.LHomologs)
}

# function to take list structure from male_l.l_coords.l_homologs 
# and produce list (per alignment) of standard plots:
#   Chromosomes All vs All
#   Unplaced Scaffolds All vs All
#   Single Chromosome Pairs (as list per Chr)
#   Homolog clusters (as list per cluster)
#   Unplaced scaffolds against matched Chromosome, as list per target Chr:
#     v.IntoRef: names of unplaced query scaffolds to place into target Chr
#     p.IntoRef: plot with alignment matches for placing query scaffolds
#     v.IntoQry: names of unplaced reference scaffolds to place into targte Chr
#     p.IntoQry: plot with alignment matches for placing reference scaffolds
# Also takes dataframe with plot axis labels:
#   rows named for alignment names
#   RefLab and QryLab columns
# Min Length cutoff parameters for Chromosome scaffolds and unplaced scaffolds to consider
plot_l.l_coords.l_homologs = function(L.L.COORD.HOM, df.Labels, 
                                      MIN_CHR = 10000000, 
                                      MIN_MATCH_CHR = 10000,
                                      MIN_UNPLACED = 100000,
                                      MIN_MATCH_UNPLACED = 10000,
                                      MIN_MATCH_PCT_UNPLACED = 0.8){
  
  l.plots <- lapply(L.L.COORD.HOM %>% names, function(ALN.name){ 
    ALN = L.L.COORD.HOM[[ALN.name]]
    # All Chromosomes
    p.Chr = PlotDFCoords(ALN$l.coords,
                         LAB_REF = paste(df.Labels[ALN.name,]$RefLab,
                                         "chromosomes"),
                         LAB_QRY = paste(df.Labels[ALN.name,]$QryLab,
                                         "chromosomes"),
                         MIN_REF_LENGTH = MIN_CHR,
                         MIN_QRY_LENGTH = MIN_CHR,
                         MIN_MATCH_LENGTH = MIN_MATCH_CHR)
    # Unplaced scaffolds (or microchromosomes)
    p.nonChr = PlotDFCoords(ALN$l.coords,
                            LAB_REF = paste(df.Labels[ALN.name,]$RefLab,
                                            "unplaced scaffolds"),
                            LAB_QRY = paste(df.Labels[ALN.name,]$QryLab,
                                            "unplaced scaffolds"),
                            MIN_MATCH_LENGTH = MIN_MATCH_UNPLACED,
                            MAX_REF_LENGTH = MIN_CHR,
                            MAX_QRY_LENGTH = MIN_CHR,
                            MIN_REF_LENGTH = MIN_UNPLACED,
                            MIN_QRY_LENGTH = MIN_UNPLACED)
    # List of each chromosome individually
    l.p.Chr = lapply(ALN$l.coords$breaksRef %>%
                       filter(SeqLength > MIN_CHR) %>%
                       rownames,
                     function(scaffold){
                       qry_scaffold = ALN$l.homologs$homologs %>%
                         filter(ref_SeqName == scaffold) %>%
                         select(qry_SeqName) %>% unlist %>% as.character()
                       p.Chr = PlotDFCoords(ALN$l.coords,
                                    REFERENCE = scaffold,
                                    QUERY = qry_scaffold,
                                    LAB_REF = df.Labels[ALN.name,]$RefLab,
                                    LAB_QRY = df.Labels[ALN.name,]$QryLab,
                                    TIC_LABELS = T)
                       return(p.Chr)
                     })
    names(l.p.Chr) = ALN$l.coords$breaksRef %>%
      filter(SeqLength > MIN_CHR) %>% rownames
    # each homology cluster
    l.p.clusters = lapply(ALN$l.homologs$clusters %>% names, function(n){
      PlotDFCoords(ALN$l.coords,
                   REFERENCE = ALN$l.homologs$clusters[[n]]$RefScafs,
                   QUERY = ALN$l.homologs$clusters[[n]]$QryScafs,
                   LAB_REF = paste(df.Labels[ALN.name,]$RefLab, n),
                   LAB_QRY = paste(df.Labels[ALN.name,]$QryLab, n) )
    })
    names(l.p.clusters) = ALN$l.homologs$clusters %>% names
    # unplaced scaffolds with good match to other assembly, for chromosome placement
    l.p.ChrPlacement = lapply(
      ALN$l.coords$breaksRef %>%
        filter(SeqLength > MIN_CHR) %>% rownames,
      function(scaffold){
        v.IntoQry = ALN$l.homologs$matches %>%
          filter(ref_SeqName == scaffold,
                 # restrict to appropriate homolog cluster, ensures no duplicates
                 qry_SeqName %in% ALN$l.homologs$clusters[[
                   which(sapply(ALN$l.homologs$clusters, 
                                function(homolog) scaffold %in% homolog$RefScaf
                                ))]]$QryScaf,
                 qry_SeqLength > MIN_UNPLACED,
                 qry_SeqLength < MIN_CHR,
                 qry_Prop > MIN_MATCH_PCT_UNPLACED) %>%
          select(qry_SeqName) %>% unlist %>% as.character()
        if( length(v.IntoQry) > 0 ){
          p.IntoQry = PlotDFCoords(ALN$l.coords,
                                   REFERENCE = scaffold,
                                   QUERY = v.IntoQry,
                                   LAB_REF = df.Labels[ALN.name,]$RefLab,
                                   LAB_QRY = df.Labels[ALN.name,]$QryLab,
                                   TIC_LABELS = T, TIC_COUNT = 100)
        } else { p.IntoQry = NULL }
        v.IntoRef = ALN$l.homologs$matches %>%
          filter(qry_SeqName == scaffold,
                 # restrict to appropriate homolog cluster, ensures no duplicates
                 ref_SeqName %in% ALN$l.homologs$clusters[[
                   which(sapply(ALN$l.homologs$clusters, 
                                function(homolog) scaffold %in% homolog$RefScaf
                   ))]]$RefScaf,
                 ref_SeqLength > MIN_UNPLACED,
                 ref_SeqLength < MIN_CHR,
                 ref_Prop > MIN_MATCH_PCT_UNPLACED) %>%
          select(ref_SeqName) %>% unlist %>% as.character()
        if( length(v.IntoRef) > 0 ){
          p.IntoRef = PlotDFCoords(ALN$l.coords,
                                   QUERY = scaffold,
                                   REFERENCE = v.IntoRef,
                                   LAB_REF = df.Labels[ALN.name,]$RefLab,
                                   LAB_QRY = df.Labels[ALN.name,]$QryLab,
                                   TIC_LABELS = T, TIC_COUNT = 100)
        } else { p.IntoRef = NULL }
        RETURN = list(v.IntoQry = v.IntoQry,
                      p.IntoQry = p.IntoQry,
                      v.IntoRef = v.IntoRef,
                      p.IntoRef = p.IntoRef)
        return(RETURN)
      })
    names(l.p.ChrPlacement) = ALN$l.coords$breaksRef %>%
      filter(SeqLength > MIN_CHR) %>% rownames
    
    RETURN = list(p.Chr = p.Chr,
                  p.nonChr = p.nonChr,
                  l.p.Chr = l.p.Chr,
                  l.p.clusters = l.p.clusters,
                  l.p.ChrPlacement = l.p.ChrPlacement)
    return(RETURN)
  })
  names(l.plots) = names(L.L.COORD.HOM)
  
  return(l.plots)
}

# function to compile all unplaced scaffold for each chromosome in all assemblies
# takes list of plots from plot_l.l_coords.l_homologs 
#   (reads l.p.ChrPlacement $v.IntoRef and $v.IntoQry)
# takes lists of positions in L.PLOTS alignments where  
# returns list (RefPri, QryPri) of lists (per Chr) of scaffolds to place
# lists for primaries are combined between all alignments
# lists for haps are taken from L.PLOTS l.p.ChrPlacement
make_l.unplaced = function(L.PLOTS,
                                        MPri_AS_REF = c(1,4,5),
                                        MPri_AS_QRY = NULL,
                                        FPri_AS_REF = c(2,3),
                                        FPri_AS_QRY = c(1)){
  CHRS = L.PLOTS[[1]]$l.p.Chr %>% names
  # for each scaffold combine unplaced with chr match in all alignemnts
  l.MpriAsRef.unplaced = sapply(CHRS, function(CHR){
                            sapply(MPri_AS_REF, 
                                   function(ALNi){
                                     ALN = L.PLOTS[[ALNi]]
                                     unlist(ALN$l.p.ChrPlacement[[CHR]]$v.IntoRef)
                                   } %>% unlist %>% as.character)
                          } %>% unlist %>% unique )
  l.MpriAsQry.unplaced = sapply(CHRS, function(CHR){
                                 sapply(MPri_AS_QRY, 
                                        function(ALNi){
                                          ALN = L.PLOTS[[ALNi]]
                                          unlist(ALN$l.p.ChrPlacement[[CHR]]$v.IntoQry)
                                        } %>% unlist %>% as.character)
                               } %>% unlist %>% unique )
  
  l.FpriAsRef.unplaced = sapply(CHRS, function(CHR){
                                 sapply(FPri_AS_REF, 
                                        function(ALNi){
                                          ALN = L.PLOTS[[ALNi]]
                                          unlist(ALN$l.p.ChrPlacement[[CHR]]$v.IntoRef)
                                        } %>% unlist %>% as.character)
                               } %>% unlist %>% unique )
  l.FpriAsQry.unplaced = sapply(CHRS, function(CHR){
                                 sapply(FPri_AS_QRY, 
                                        function(ALNi){
                                          ALN = L.PLOTS[[ALNi]]
                                          unlist(ALN$l.p.ChrPlacement[[CHR]]$v.IntoQry)
                                        } %>% unlist %>% as.character)
                               } %>% unlist %>% unique )
  # remove duplicates, order by scaffold name
  l.Mpri.unplaced = sapply(CHRS, function(CHR){
    CHR.unplaced = c(l.MpriAsRef.unplaced[[CHR]], l.MpriAsQry.unplaced[[CHR]])
    CHR.unplaced = CHR.unplaced[order(nchar(CHR.unplaced), CHR.unplaced)] %>% unique
    return(CHR.unplaced)
  })
  l.Fpri.unplaced = sapply(CHRS, function(CHR){
    CHR.unplaced = c(l.FpriAsRef.unplaced[[CHR]], l.FpriAsQry.unplaced[[CHR]])
    CHR.unplaced = CHR.unplaced[order(nchar(CHR.unplaced), CHR.unplaced)] %>% unique
    return(CHR.unplaced)
  })
  
  # names of haps from L.PLOT, only primary vs primary represented twice in ALNi sets
  PRIi = c(MPri_AS_QRY, MPri_AS_REF, FPri_AS_QRY, FPri_AS_REF)[
    duplicated(c(MPri_AS_QRY, MPri_AS_REF, FPri_AS_QRY, FPri_AS_REF))
  ]
  # list for each hap alignment, always Qry against relevant primary
  l.haps.unplaced = lapply(names(L.PLOTS)[-PRIi],
                           function(ALNi){
                             ALN = L.PLOTS[[ALNi]]
                             sapply(CHRS, function(CHR){
                               hap.unplaced = ALN$l.p.ChrPlacement[[CHR]]$v.IntoQry
                               hap.unplaced = hap.unplaced[order(nchar(hap.unplaced),
                                                                 hap.unplaced)]
                             })
                           })
  names(l.haps.unplaced) = names(L.PLOTS)[-PRIi]
  l.unplaced = c(list(Mpri = l.Mpri.unplaced,
                      Fpri = l.Fpri.unplaced),
                 l.haps.unplaced)
  return(l.unplaced)
}


# plot_l.unplaced_self() takes l.unplaced (from make_l.unplaced_in_primaries) 
# and l.self_align (from make_l.l_coords.l_homologs on self alignments)
#   names for l.self_align primaries need to be "Fpri" and "Mpri" to match l.unplaced
# produces list (per primary Fpri/Mpri) of plot lists (per scaffold)
# of unplaced scaffolds against chromosomes from self alignment
plot_l.unplaced_self = function(L.UNPLACED, L.SELF, DUP_THRESHOLD = 0.5){
  l.unplaced.self.plots = lapply(
    L.UNPLACED %>% names, 
    function(ASM){
      p.asm = lapply(L.UNPLACED[[ASM]] %>% names,
                     function(SCAF){
                       if( length(L.UNPLACED[[ASM]][[SCAF]]) == 0 ){
                         return(NULL)
                       } else {
                       p = PlotDFCoords(L.SELF[[ASM]]$l.coords,
                                    REFERENCE = SCAF,
                                    QUERY = L.UNPLACED[[ASM]][[SCAF]],
                                    TIC_LABELS = T, TIC_COUNT = 100,
                                    LAB_REF = ASM, LAB_QRY = ASM)
                       d = L.SELF[[ASM]]$l.homologs$matches %>%
                         filter(qry_SeqName %in% L.UNPLACED[[ASM]][[SCAF]],
                                qry_Prop > DUP_THRESHOLD) %>%
                         select(qry_SeqName) %>% unlist
                       return(list(plot = p, dups = d))
                       }   
                     })
      names(p.asm) = L.UNPLACED[[ASM]] %>% names
      return(p.asm)
    })
  names(l.unplaced.self.plots) = L.UNPLACED %>% names
  return(l.unplaced.self.plots)
}


# function to convert lists into dataframe with unplaced and duplicate scaffolds 
#  one row for each assembly/chromosome
#  $unplaced and $dups columns: single comma sep character string of non-Chr scaffolds
# takes l.unplaced and l.unplaced.self.plots
make_df.unplaced_dups = function(L.UNPLACED, L.SELF.PLOTS){
  data.frame(
    asm = rep(names(L.UNPLACED), each = length(L.UNPLACED[[1]])),
    seq = rep(names(L.UNPLACED[[1]]), length(L.UNPLACED)),
    unplaced = lapply(L.UNPLACED, function(ASM){
      sapply(ASM, paste, collapse = ",") %>% unlist %>% as.character()
    }) %>% unlist,
    dups = lapply(L.SELF.PLOTS, function(ASM){
      lapply(ASM, function(SCAF){
        paste(SCAF$dups, collapse = ",") %>% unlist %>% as.character()
      })
    }) %>% unlist,
    row.names = NULL
    # remove duplicated scaffolds from unplaced list
  ) %>% mutate(unplaced = str_remove(unplaced, 
  # throws warnings if search string empty so adding "filler" to avoid warning
                                     str_replace_all(paste("filler", 
                                                           dups, sep = ","), 
                                                     ",", "|")), 
  # clean up leftover commas
               unplaced = str_replace_all(unplaced, ",,", ","), 
               unplaced = str_remove(unplaced, "^,|,$"))
  
}

