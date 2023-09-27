run_DiffExp <- function( inputData, Var1, Var2 = NULL, mode = "1 VS rest", 
                         m1v1fg = NULL, m1v1bg = NULL, min.pct = 0.25, minLFC = 0.5, 
                         minCells = 20, rmMitoG =TRUE, rmRiboG = TRUE
){
  # Main function to calculate diferential expression
  ## ins:
  # -  inputData: input dataset, Seurat Object
  # - Var1: metadata variable separating cells (cluster)
  # - Var2: optional, for conditional DE (i.e, DE due to treatment in each cell type )
  # - mode: mode of analysis, one of ['1 VS res'|'1 VS 1'|'Conditional']
  #   - m1v1fg (only in 1v1 mode): foreground group
  #   - m1v1bg (only in 1v1 mode): background group
  # - min.pct: minimum percentage of cells in foreground or background expressing
  #     the gene.
  # - minLFC: minimum value of log2(Fold Change) to consider differential expression
  #     relevant.
  # - minCells: minimum number of cells in either foreground or background, otherwise 
  #     the comparison is discarded.
  # - rmMitoG and rmRiboG: remove or keep mitocondrial and ribosomal genes
  ## outs:
  # - list of:
  #   - Parameters: list of parameters used.
  #   - Tables: differential expression data.frames
  #   - Plots: differential expression plots:
  #     - Heatmap: top 5 up-reg genes in each group
  
  # Set undefined variables to null
  if( Var2 == "None" ){ Var2 <- NULL }
  if(!is.null(m1v1fg)){if( m1v1fg == "None" ){ m1v1fg <- NULL }}
  if(!is.null(m1v1bg)){if( m1v1bg == "None" ){ m1v1bg <- NULL }}
  if( min.pct > 1 ){ min.pct <- min.pct/100 }
  
  # Set cell identity
  
  # drop empty levels of the grouping variables
  if(! is.null ( Var1 ) ){
    if( is.factor( inputData[[ Var1 ]][,1] )){
      inputData[[ Var1 ]][,1] <- droplevels( inputData[[ Var1 ]][,1]  )  } }
  if(! is.null ( Var2 ) ){
    if( is.factor( inputData[[ Var2 ]][,1]) ) {
      inputData[[ Var2 ]][,1] <- droplevels( inputData[[ Var2 ]][,1]  )   }}
  
  inputData <- DiffExpr_cellIdent( inputData, Var1, Var2 )
  groups <- sort( unique( inputData$groups ) )
  
  # Switch active assay to normalized counts (CPM)
  inputData <- DiffExpr_CPMnorm( inputData )
  defAssay <- "RNA"
  
  # set features of interest
  feats <- rownames( inputData )
  
  if( rmMitoG | rmRiboG ){ feats <- DiffExpr_rmFeats( feats, rmMitoG, rmRiboG ) }
  
  groupList <- DiffExpr_setComparisons( inputData, mode, Var1, Var2, m1v1fg, m1v1bg )
  
  markers <- lapply(
    groupList,
    function(i, ptc, cellDistr ){
      groups <- as.character( unlist(i) )
      
      if( all( cellDistr[ groups ] > minCells ) ){
        mark <- FindMarkers( inputData, ident.1= i$fg , ident.2=i$bg, min.pct=min.pct,
                             features = feats, logfc.threshold =minLFC, assay = defAssay  )
        
      } else { NULL }
      
    }, ptc=min.pct, cellDistr = table(inputData$groups)
  )
  names(markers) <- sapply( groupList, function(i){
    if( is.null( i$bg) ){ i$bg <- "rest"}
    paste0( i$fg, " VS ", i$bg ) } )
  
  # Make parameter table
  if( is.null(Var2) ){ Var2 <- "None" }
  if( is.null(m1v1fg) ){ m1v1fg <- "None" }
  if( is.null(m1v1bg) ){ m1v1bg <- "None" }
  params <- data.frame(
    Value=c( inputData@project.name, Var1, Var2, mode, m1v1fg, m1v1bg, min.pct, minCells, minLFC,
             rmMitoG, rmRiboG),
    row.names = 
      c("Dataset", "Grouping Variable", "Secondary grouping variable", "Mode of comparison",
        "Foreground group", "Background group", "Minimum % of cells expressing a gene",
        "Minimum of cells per group", "Minimum Log2FC", "Discard mitocondrial genes",
        "Discard ribosomal genes")
  )
  params <- params[ params$Value != "None", , drop=FALSE]
  
  # Make summary table
  sumTab <- DiffExpr_summary( 
    inputData, markers, Var1, Var2, mode, minCells, minLFC, m1v1fg, m1v1bg, groupList )
  
  g <- unique( unlist( groupList, recursive = TRUE ) )
  
  if(mode=="1 VS rest" ){
    plotHM <- DiffExpr_plotHM( inputData, markers, groups=g, minLFC=minLFC )
    gc()
  } else {
    plotHM <- NULL
  }
  
  return( list( Parameters=params, Tables=markers, Summary=sumTab,
                Plots = list( Heatmap=plotHM ) ) )
}

##### auxiliary functions #####

DiffExpr_rmFeats <- function( feats, rmMitoG, rmRiboG, rmOtherG=NULL ){
  # Returns the array of gene IDs 'genes' after removing selected features
  ## ins:
  # - feats: array of gene IDs
  # - rmMitoG (T/F): remove mitocondrial genes
  # - rmRiboG (T/F): remove ribosomal genes
  # - rmOtherG: array of other IDs to remove 
  ## outs:
  # list of features (gnene IDs) facter removing ribosomal/mitocondrial/other
  # gene IDs.
  
  spec <- DiffExpr_genomeSpecies( feats )
  
  spec[1] <- gsub("\\. ", "", spec[1] )
  avblOrgs <- unique( gsub( "\\.txt", "", gsub( "^[A-Za-z]*_", "", dir("data/QC.Genes/") ) ) )
  
  
  if( spec[1] %in% avblOrgs ){  # skip if species is "Other" or uniavailable in the data folder
    if( spec[ 2 ] != "ENSEMBL Gene" ){  # gene symbols
      GeneList <- gen_loadGeneLists( IDformat = "symbols", section = "Rb.Mt.NRY"  )
    } else { # ENSEMBL gene IDs
      GeneList <- gen_loadGeneLists( IDformat = "ENSids", section = "Rb.Mt.NRY"  )
    }
    
    GeneList <- GeneList[[ spec[ 1 ] ]]
    
    if( rmMitoG ){ feats <- feats[ ! feats %in% GeneList[[ "mitoG" ]] ] }
    if( rmRiboG ){ feats <- feats[ ! feats %in% GeneList[[ "riboG" ]] ] }
    
  }
  if( !is.null(rmOtherG) ){ feats <- feats[ ! feats %in% rmOtherG ] }
  
  
  return(feats)
}

DiffExpr_cellIdent <- function( inputData, Var1, Var2 = NULL, m1v1fg = NULL ){
  # In a Seurat object, set cell identity according to variables 'Var1' & 'Var2'.
  ## ins: 
  #  - inputData: scRNA-seq dataset, Seurat object.
  #  - Var1 and Var2, names of columns in inputData@meta.data
  ## outs: 
  #  - Seurat object with updated cell identities and removed cells with NA
  #    values for 'Var1' and 'Var2'
  
  if( !is.null(Var2) ){ if( is.na(Var2) ){ Var2 <- NULL } }
  
  if( ! is.null( Var2 ) ){
    keep <- !is.na( inputData@meta.data[[ Var1 ]] ) & !is.na( inputData@meta.data[[ Var2 ]] )
    inputData <- subset( inputData, cells=colnames(inputData)[ keep ] )
    inputData[["groups"]] <- paste0( 
      inputData@meta.data[[ Var1 ]], " - ", inputData@meta.data[[ Var2 ]] )
    
  } else { 
    keep <- !is.na( inputData[[Var1]] )
    inputData <- subset( inputData, cells=colnames( inputData )[ keep ] )
    inputData[["groups"]] <- inputData@meta.data[[ Var1 ]]
    groups <- unique( inputData@meta.data[[ Var1 ]] ) 
  }
  Idents( inputData ) <- "groups"
  return( inputData )
}

DiffExpr_sortGroups <- function( x ){
  # funtion to return sorted levels of x -in case x is a char, but starts
  # with numbers and so on, to avoid ['1-...','10-...','2-...' ]
  
  if (is.factor(x)){ x <- as.character(x) }
  x <- unique(x)
  
  if( all( !is.na( as.numeric( x ) ) ) ){ # if group names are numbers
    x_levels <- sort( as.numeric( x ) )
  } else {
    if( all( grepl(" - ", x ) ) ){ 
      if( all( !is.na( as.numeric( gsub(" -.*", "", x ) ) ) ) ){ # if group names start with numbers
        x_levels <- x[ order( gsub(".*- ", "", x ) ) ][
          order(x[ order(  gsub(".*- ", "", x ) )  ]) ]
        
      } else if( all( !is.na( as.numeric( gsub(".*- ", "", x ) ) ) )  ){ # if group names end with numbers
        x_levels <-  x[ order( gsub(" -.*", "", x ) ) ][
          order( as.numeric( gsub( ".*- ", "",  x[ order( gsub(" -.*", "", x ) ) ] ) ) ) ]
      } else {
        x_levels <- sort( x ) 
      }
    } else {  x_levels <- sort( x ) }
  }
  return( x_levels )  
}

DiffExpr_setComparisons <- function( tmp, mode, Var1, Var2, m1v1fg, m1v1bg ){
  # Function define comparisons for differential expression analysis according to the
  # grouping variables set by the user.
  ## ins:
  # tmp: Seurat object
  # mode: mode of analysis, one of ['1 VS res'|'1 VS 1'|'Conditional']
  # Var1 and Var2: grouping variables, column names in tmp@meta.data
  # m1v1fg: in mode '1 VS 1', foreground group
  # m1v1bg: in mode '1 VS 1', background group
  
  ## outs:
  # list of comparisons. Each element of the list contains the name of
  # a foreground and background group (if mode is '1 VS rest', background is set to NULL)
  
  groups <- DiffExpr_sortGroups( tmp$groups )
  
  if( mode == "1 VS rest" ){
    groupList <- lapply(
      groups, function(i){
        list( fg = i, bg=NULL )
      } )
    
  }else if( mode=="Conditional" ){
    
    v2groups <- DiffExpr_sortGroups( tmp[[Var2]][,1] )
    v1groups <- sort( unique( tmp[[Var1]][,1] ) )
    if( !is.null( m1v1fg ) ){
      v1groups <- c( m1v1fg, v1groups[ v1groups != m1v1fg ] )
    }
    
    groupList <- lapply(
      v2groups,
      function(i, v1groups ){
        subGroups <- groups[ groups %in% tmp$groups[ tmp[[Var2]][,1]==i ] ]
        dimNa <- subGroups[ match( v1groups, gsub( " -.*", "", subGroups)  ) ]
        m <- matrix( T, length(subGroups), length(subGroups), dimnames = list( dimNa,dimNa ) )
        m[ lower.tri(m, T) ] <- F
        
        com <- lapply( seq_len( nrow( m ) ), function(i){
          com <- lapply( seq_len( ncol( m ) ), function(j){
            if( m[i, j] ){ list( fg=rownames(m)[i], bg=rownames(m)[j] ) } 
          })
          com[ sapply( com, is.null) ]<- NULL
          return(  com  )
        })
        com <- unlist( com, recursive = FALSE )
        return( com )
      }, v1groups = v1groups
    )
    groupList <- unlist( groupList, recursive=FALSE )
    groupList[ sapply( groupList, is.null) ]<- NULL
    
  }else if( mode == "1 VS 1" ){
    if(is.null( m1v1fg) & is.null( m1v1bg ) ){ # make all combinations
      m <- matrix( T, length(groups), length(groups), dimnames = list(groups,groups) )
      m[ lower.tri(m, T) ] <- F
      groupList <- lapply( seq_len( nrow( m ) ), function(i){
        com <- lapply( seq_len( ncol( m ) ), function(j){
          if( m[i, j] ){ list( fg=rownames(m)[i], bg=rownames(m)[j] ) }
        })
        com[ sapply( com, is.null) ]<- NULL
        return(com)
      })
      groupList <- unlist( groupList, recursive=FALSE )
      
    } else if( !is.null( m1v1fg ) & is.null( m1v1bg ) ){
      groups <- groups[ groups != m1v1fg ]
      groupList <- lapply( groups, function(i){ list( fg= m1v1fg, bg = i ) })
      
    } else {
      groupList <- list( list( fg = m1v1fg, bg = m1v1bg ) )
    }
  }
  return( groupList )
  
}

DiffExpr_genomeSpecies <- function( x ){
  # Check the feature ids to determine the species
  ## ins:
  # x: character vector - gene IDs
  ## outs: 
  # character vector containing dataset species and ID type [symbol | ensemble GeneID]
  
  if( all( grepl("^ENS", x ) ) ){ # ensembl IDs
    genePatts <- gen_loadGeneLists( "ENSids", "patterns" )
    matches <- sapply(
      genePatts, function(i, x){ sum( grepl( i, x ) ) } ,x=x )
    return( c(
      Species = names(genePatts)[ which( matches == max(matches) ) ],
      IDtype = "ENSEMBL Gene"
    ))
  } else {
    genePatts <- gen_loadGeneLists( "symbols", "patterns" )
    symbolOrgs = c( "HGNC", "MGI", "RGNC", "ZFIN" )
    matches <- sapply(
      genePatts,
      function(i, x){
        sum( sum( grepl( i[1], x ) ) + sum( grepl( i[2], x ) ) )
      }, x=x
    )
    return( c(
      Species = names(genePatts)[ which( matches == max(matches) ) ],
      IDtype = paste0( symbolOrgs[which( matches == max(matches) )], " gene symbols")
    ))
  }
  
}

DiffExpr_CPMnorm <- function( inputData ){
  # Return normalized (CPM) RNA counts as the default assay of inputData.
  ## ins: 
  #     - inputData: scRNA-seq dataset, Seurat object.
  ## outs: 
  #     - Seurat object with DefaultAssay set to "RNA" and counts
  #       normalized to CPM values.
  
  DefaultAssay( inputData ) <- "RNA"
  NormalizeData( 
    inputData,
    normalization.method = "RC",
    scale.factor = 1e6 )
}

DiffExpr_plotHM <- function( inputData, markers, groups=NULL, Var1=NULL, Var2=NULL, n=5, minLFC=0.5 ){
  # Plots heatmap of top DEGs for every cell group. If the DE results are loaded
  # from a file, cell groups will be generated from the variables Var1 and Var2.
  
  ## ins:
  # - inputData: Seurat object
  # - markers: list of diferential expression results
  # - groups: character vector, combinations of values or the grouping variables
  # - Var1 and Var2: grouping variables, column names in inputData@meta.data
  # - n: number of top genes to plot per cell group
  # - minLFC: minimum value of log2(Fold Change) to consider differential expression
  #     relevant.
  
  if( is.null( groups ) ){ # if loading DE results from file
    inputData <- DiffExpr_cellIdent( inputData, Var1, Var2 )
    groups <- levels( inputData$groups )
  }
  
  cells <- colnames(inputData)[ inputData$groups %in% groups ]
  inputData <- subset( inputData, cells = cells)
  
  inputData$groups <- factor( as.character( inputData$groups ), levels=groups )
  
  if(length(markers)>1){
    top <- unlist( lapply(
      markers,
      function(i){
        i <- i[ i$avg_log2FC > minLFC, ]
        if(nrow(i) >= n ){ rownames( i )[ 1:n ]
        } else { rownames(i)  }  } ))
  } else {
    top <- c(
      rownames(markers[[1]])[ markers[[1]]$avg_log2FC>minLFC ][1:n],
      rownames(markers[[1]])[ markers[[1]]$avg_log2FC< -minLFC ][1:n]
    )
  }
  
  if( DefaultAssay( inputData ) != "RNA" ){
    DefaultAssay( inputData ) <- "RNA"
    cat("Switching default assay to 'RNA' to generate DEG plots.\n")
  }
  inputData <- NormalizeData( inputData, normalization.method = "LogNormalize",
                              scale.factor = 10000 )
  inputData <- ScaleData( inputData, features = rownames(inputData) )
  
  DoHeatmap( inputData, features = top, size = 4, raster = TRUE, group.by = "groups"  ) + NoLegend()
}

DiffExpr_summary <- function( inputData, markers, Var1=NULL, Var2=NULL,
                              mode="1 VS rest", minCells=20,  minLFC=0.5,
                              m1v1fg=NULL, m1v1bg=NULL, groupList=NULL ){
  # Function to create a summary data frame containing details of the DE analysis.
  ## ins:
  # - inputData: Seurat object
  # - markers: list of diferential expression results
  # - Var1 and Var2: grouping variables, column names in inputData@meta.data
  # - mode: mode of analysis, one of ['1 VS res'|'1 VS 1'|'Conditional']
  # - minLFC: minimum value of log2(Fold Change) to consider differential expression
  #     relevant.
  # - minCells: minimum number of cells in either foreground or background, otherwise 
  #     the comparison is discarded.
  # - m1v1fg (only in 1v1 mode): foreground group
  # - m1v1bg (only in 1v1 mode): background group
  # - groupList: list of comparisons to make, output of DiffExpr_setComparisons()
  
  ## outs: data frame containing details on the DE results
  
  # If this function is run outside of 'run_diffExpr()':
  if( !is.null(Var2) ){ if( is.na(Var2) ){ Var2 <- NULL } }
  if( !is.null(m1v1fg) ){ if( is.na(m1v1fg) ){ m1v1fg <- NULL  } }
  if( !is.null(m1v1bg) ){ if( is.na(m1v1bg) ){ m1v1bg <- NULL  } }
  if( !is.numeric( minLFC ) ){ minLFC <- as.numeric( minLFC ) }
  
  
  if( ! "groups" %in% colnames( inputData@meta.data )){
    inputData <- DiffExpr_cellIdent( inputData, Var1, Var2 )
  }
  
  cellsPerGroup <- table( inputData@meta.data$groups )
  
  if( is.null(groupList) ){
    groupList <- DiffExpr_setComparisons( inputData, mode, Var1, Var2, m1v1fg, m1v1bg )
  }
  
  if( length( groupList ) != length( markers ) ){ # loading DE from file, if some comparison didn't meet cell # minimum
    expectedM <- sapply( groupList, paste0, collapse=" VS " )
    missingM <- expectedM[ ! expectedM %in% names(markers) ]
    markers <- c( markers, list( missingM ))
    markers[ missingM ] <- data.frame( NA )
  }
  
  sumTab <- data.frame(
    "Comparison" = names(markers),
    "Cells in Foreground" = as.integer( cellsPerGroup[ as.character(sapply(groupList, "[[", "fg" ) )]),
    "Cells in Background" = rep(NA, length(markers)),
    "Total DEG" = sapply( markers, function(i){sum(i$p_val_adj<0.05)}),
    "Up-regulated Genes" =  sapply( markers, function(i){sum(i$p_val_adj<0.05 & i$avg_log2FC>minLFC) } ),
    "Down-regulated Genes" =  sapply( markers, function(i){sum(i$p_val_adj<0.05 & i$avg_log2FC< -minLFC)})
  )
  
  if( mode =="1 VS rest" ){
    sumTab$Cells.in.Background <- as.integer( ncol(inputData) - sumTab$Cells.in.Foreground )
  } else {
    sumTab$Cells.in.Background <-  as.integer( cellsPerGroup[ sapply(groupList, "[[", "bg" ) ] )
  }
  
  sumTab[( sumTab$Cells.in.Foreground < minCells |
             sumTab$Cells.in.Background < minCells), 4:6 ]  <- NA
  
  return( sumTab )
}


