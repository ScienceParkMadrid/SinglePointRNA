run_paths <- function( inputData, DEtabs, pathDB, DEparams, minDE=5 ){
  # Function to run enrichment analysis of differentially expressed genes on known metabolic pathways
  ## ins:
  # - inputData: Seurat object
  # - DEtabs: list of differential gene expression result data frames, output of 
  #     run_DEtabs().
  # - pathDB: list of gene sets (char. arrays) containing genes involved in metabolic pathways.
  # - DEparams: parameters used to calculate differential expression.
  # - minDE: minimum number of DE genes in a differential expression analysis to proceed with
  #     pathway analysis.
  
 
  # set background array of expressed genes.
  bg_genes <- paths_getBG( inputData, DEparams, DEtabs )
  minLFC <- as.numeric( DEparams[ "Minimum Log2FC", ] )
  
  
  if( DEparams[ "Mode of comparison", ] %in% c("1 VS 1", "Conditional") ){ 
    pathEnrich <- bplapply( 
      seq_along(DEtabs),
      function(i){
        
        # filter out pathway gene sets with less than 10 genes in either foregrounds or backgounds.
        allG <- unique( unlist( c( rownames(DEtabs[[i]]), bg_genes[[i]] )  ) )
        pathDB <- pathDB[ sapply(  
          pathDB, function(i){ sum(i %in% allG )>=10 }  ) ]
        
        tab <- paths_Fisher( DEtabs[[i]], pathDB = pathDB, bg = bg_genes[[ i ]], LFCthreshold = minLFC )
    }, BPPARAM = MulticoreParam(workers= 4 ) )
    
  } else { # 1 VS rest DE analysis
    allG <- unique( unlist( c( lapply(DEtabs,rownames), bg_genes )  ) )
    pathDB <- pathDB[ sapply(  
      pathDB, function(i){ sum(i %in% allG )>=10 }  ) ]
    
    pathEnrich <- lapply( 
      DEtabs, paths_Fisher, pathDB = pathDB, bg = bg_genes, LFCthreshold=minLFC  )
  }
  names( pathEnrich ) <- names( DEtabs )
  
  enrSum <- paths_summary(  inputData, DEtabs, DEparams, pathEnrich, bg_genes )
  
  return( list( Tables = pathEnrich, Summary=enrSum, backGrounds = bg_genes ) )
}

#### Auxiliary functions ####

paths_loadDBs <- function(  iDsummary, collections ){
  # function to select and load a patway genese database
  # based on the gene ID type and species of the dataset
  ## ins:
  # iDsummary: data frame summarizing characteristics of the input dataset.
  #   output of load_getInputReport()
  ## outs:
  # 
  
  if( grepl( "^ENSEMBL", iDsummary[1, "Feature ID type"] ) ){
    mSig_Paths <- gen_loadGeneLists("ENSids", "pathways"  )
  } else {
    mSig_Paths <- gen_loadGeneLists("symbols", "pathways"  )
  }
  
  mSig_Paths <- mSig_Paths[ grep( gsub("\\. ", "", iDsummary$Species[1] ), names(mSig_Paths) ) ]
  names( mSig_Paths) <- gsub( "_[a-zA-Z]*", "", names(mSig_Paths) )
  
  if( collections !="All"){
    mSig_Paths <- mSig_Paths[ collections ]
  }
  
  mSig_Paths <- unlist( mSig_Paths, recursive = FALSE )
  names(mSig_Paths) <- gsub("^.*\\.", "", names(mSig_Paths))
  
  return( mSig_Paths )
}

paths_selectDEresult <- function( DEGtables, DEGselected ){
  # function to subset the DE results to analyze.
  ## ins:
  # DEGtables: list of data frames containing the differential expression results.
  #   Part of the output of run_DiffExp()
  ## outs: list of data frames. Subset of DEGtables selected by the user.
  
  if( DEGselected != "All"  ){
    return( DEGtables[ names( DEGtables) %in% DEGselected ] )
  } else{ 
    return( DEGtables )
  }
}

paths_getBG <- function( inputData, DEparams, DEtabs, laxBG=TRUE ){
  # function to gen genes expressed in at least min.ptc% of cells of at least
  # one cluster.
  ## ins:
  # - inputData: Seurat object
  # - DEparams: data frame of parameters used for DEG analysis
  # - DEtabs: differential expression results
  # - laxBG: BG with any genes that are expressed, if FALSE, limit to genes 
  #   expressed in at least a % of the cells of each group.
  ## outs: character vector or list of character vectors of background genes.
  
  # get DE parameters
  deMode <- DEparams[ "Mode of comparison", ]
  min.ptc <- as.numeric( DEparams["Minimum % of cells expressing a gene",] )
  Var1 <- DEparams["Grouping Variable",]
  Var2 <- DEparams["Secondary grouping variable",]
  
  DefaultAssay( inputData ) <- "RNA"
  
  if( !"groups" %in% colnames(inputData@meta.data) ){
    inputData <- paths_cellIdent( inputData, Var1, Var2 )
  }
  
  # genes in DE results will be accepted automatically
  # if minLFC > 0, examine the rest of the genes
  if( deMode =="1 VS rest"){
    acceptedG <- unique( unlist(lapply( DEtabs, rownames ) ) )
    acceptedG <- acceptedG[ !is.na(acceptedG ) ]
    groupList <- sort( unique( inputData$groups ) )  
  } else if( deMode %in% c( "1 VS 1","Conditional") ){
    acceptedG <- lapply( DEtabs, rownames )
    names( acceptedG ) <- gsub( ".*- ","", names(acceptedG) )
    groupList <- lapply( names( DEtabs ), function(i){ unlist( strsplit( i, " VS ") ) } )
  }
  
  if( is.factor( inputData$groups ) ){
     levels( inputData$groups ) <- c( levels(inputData$groups), "Other" )
  }
  
  if( DEparams[ "Minimum Log2FC", ] > 0 & laxBG ){
    
    expressedGenes <- lapply(
      seq_along(groupList),
      function( i, inputData, genesNoTest, deMode ){
        
    
        if( deMode=="1 VS rest" ){
          i <- levels( groupList )[ i ]
          
          inputData$groups[ inputData$groups !=i ] <- "Other"
          g <- rownames(inputData)[ ! rownames(inputData) %in% genesNoTest ]
          g <- AverageExpression(
            inputData, assays = "RNA", slot="counts", features = g, group.by = "groups" )[[1]] 
          g <-  rownames( g )[ apply( g, 1, function( j ){ any( j > 0 ) } ) ]
          gc()
          
        } else if ( deMode %in% c("1 VS 1","Conditional") ){
          groupList_i <- groupList[[ i ]]
          gnt_i <- genesNoTest[[ i ]]
          
          cellIDs <- colnames(inputData)[ inputData$groups %in% groupList_i ]
          g <- rownames(inputData)[ ! rownames(inputData) %in% gnt_i ]
          
          g <- AverageExpression( subset(inputData,cells = cellIDs), 
            assays = "RNA", slot="counts", features = g, group.by = "groups" )[[1]]
          g <- rownames( g )[ apply( g, 1, function( j ){ any( j > 0 ) } ) ] 
          gc()
        }
        return( g )
      }, 
      inputData = inputData, genesNoTest = acceptedG, deMode=deMode
    ) 
    
    
  } else if( DEparams[ "Minimum Log2FC", ] > 0 & ! laxBG ) {
    
    expressedGenes <- lapply(
      seq_along(groupList),
      function( i, inputData, genesNoTest, deMode, min.ptc ){
        
        if( deMode=="1 VS rest" ){
          
          fg <- colnames(inputData)[ inputData$groups == groupList[ i ] ]
          bg <- colnames(inputData)[ !colnames(inputData) %in% fg ]
          
          g <- rownames(inputData)[ ! rownames(inputData) %in% genesNoTest ]
          minCell = min.ptc * length( fg )
          
          fg <- FetchData( inputData, cells = cfg, vars = g, slot = "counts" ) > 0
          fg <-  colnames( fg )[ colSums( fg ) > minCell ]
          bg <- FetchData( inputData, cells = cbg, vars = g[!g %in% fg ], slot = "counts" ) > 0
          bg <-  colnames( bg )[ colSums( bg ) > minCell ]
          g <- c(fg,bg)
          gc()
          
        } else if ( deMode %in% c("1 VS 1","Conditional") ){
          
          gnt_i <- genesNoTest[[ i ]]
          
          fg <- colnames(inputData)[ inputData$groups == groupList[ i ][ 1 ] ]
          bg <- colnames(inputData)[ inputData$groups == groupList[ i ][ 2 ] ]
          
          fg <- FetchData( inputData, cells = cfg, vars = g, slot = "counts" ) > 0
          fg <-  colnames( fg )[ colSums( fg ) > minCell ]
          bg <- FetchData( inputData, cells = cbg, vars = g[!g %in% fg ], slot = "counts" ) > 0
          bg <-  colnames( bg )[ colSums( bg ) > minCell ]
          g <- c(fg,bg)
          
          gc()
        }
        return( g )
      }, 
      inputData = inputData, genesNoTest = acceptedG, deMode=deMode, min.ptc=min.ptc
    )
    
  }
  
  bgList <- lapply(
    seq_along(expressedGenes), function(i){
      if( is.list( acceptedG ) ){c( acceptedG[[i]], expressedGenes[[i]])
      } else { c( acceptedG, expressedGenes[[i]]) }
    })
  names( bgList ) <- names( DEtabs )
  return( bgList )
}

paths_cellIdent <- function( inputData, Var1, Var2 = NULL ){
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

paths_sortGroups <- function( x ){
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

paths_Fisher <- function(deg, pathDB, bg, pVthreshold=0.05, LFCthreshold=0.5,
                         upOnly=F, downOnly=F ){
  # Function to calculate enrichment of DEG genes in a given list of genesets.
  ## ins: 
  # - deg: data frame of differential expression results.
  # - pathDB: list of character vectors, relevant genesets
  # - bg: background genes
  # - pVthreshold:
  # - LFCthreshold:
  # - upOnly and DownOnly: for "1 VS 1" differential expression, reporte results of
  #     genes DE in group 1 (upOnly) or group 2 (donwOnly).
  
  
  if( upOnly ){ deg_g <- rownames( deg )[ deg$avg_log2FC > 0 &deg$p_val_adj < pVthreshold ] 
  } else if ( downOnly ){ 
    deg_g <- rownames( deg )[ deg$avg_log2FC < 0 & deg$p_val_adj < pVthreshold ]  
  } else { deg_g <- rownames( deg )[ deg$p_val_adj < pVthreshold ]
  }
  
  if( length(deg_g) > 0 ){
    bg <- bg[ !bg %in% deg_g ]
    
    enrich_fisher <- do.call( rbind, lapply(
      pathDB,
      function(i, deg_g, bg){
        
        crossTab <- rbind(
          c( sum( deg_g %in% i ), sum( !deg_g %in% i )),
          c( sum( bg %in% i ), sum( !bg %in% i ))
        )
        
        fish <- fisher.test( crossTab )
        data.frame(
          OddsRatio = fish[["estimate"]],
          pval = fish[["p.value"]],
          nDEG = crossTab[1,1],
          DEG = paste0( deg_g[ deg_g %in% i], collapse = ", "),
          list_length= length(i),
          list_present = sum( i %in% c(bg,deg_g) )
          
        )
      }, deg_g = deg_g, bg=bg
    ) )
    
    enrich_fisher$list <- names( pathDB )
    enrich_fisher <- enrich_fisher[ enrich_fisher$DEG > 0, ]
    enrich_fisher$p_val_adj <- p.adjust( enrich_fisher$pval, "fdr" )
    enrich_fisher <- enrich_fisher[ order(enrich_fisher$pval ), ]
    rownames( enrich_fisher ) <- NULL
    enrich_fisher
    
  } else { NULL }
  
}

paths_summary <- function( inputData, DEtabs, DEparams, pathwayEnrichment, bg_genes,
                           pVthreshold=0.05 ){
  # Function so summarize results of DEG enrichment in metabolic pathways
  ## ins:
  # - inputData: Seurat object
  # - DEtabs: list of data frames, results of diferential analysis
  # - DEpartams: data frame, parameters used in DE analysis
  # - pathwayEnrichment: list of data frames, results of paths_Fisher()
  # - bg_genes: character vector or list of expressed genes
  # - pVthreshold: maximum adjusted p-value to consider a gene DE and a pathway
  #     enriched in DEGs.
  ## outs: data frame containing some details of the dataset, differential expression,
  #    and DEG enrichment in molecular pathways.
  
  # get DE parameters
  deMode <- DEparams[ "Mode of comparison", ]
  min.ptc <- as.numeric( DEparams["Minimum % of cells expressing a gene",] )
  Var1 <- DEparams["Grouping Variable",]
  Var2 <- DEparams["Secondary grouping variable",]
  
  
  inputData <- paths_cellIdent( inputData, Var1, Var2 )  
  groups <- sort( unique( inputData$groups ) )    
  groupSize <- table(inputData$groups)
 
  groupList <- strsplit( names( DEtabs ), " VS "  )
  
  groupData <- do.call(rbind, lapply(
    groupList,
    function(i, deMode ){
      comp <- paste0( i, collapse = " VS "  )
      
      if(is.list( bg_genes )){ exprG <- length( bg_genes[[ comp ]]  )
      } else { exprG <- length( bg_genes ) }
      
      deg <- sum( DEtabs[[ comp ]]$p_val_adj <= pVthreshold )
      
      if( !is.null( pathwayEnrichment[[comp]] )){
        totalPaths <- nrow( pathwayEnrichment[[ comp ]] )
        signifPaths <- sum( 
        pathwayEnrichment[[comp]]$p_val_adj <= pVthreshold &
          pathwayEnrichment[[ comp ]]$OddsRatio > 1 )
      } else {
        totalPaths <- NA
        signifPaths <- NA
      }
      
      data.frame( row.names = NULL,
        Comparison = comp,
        Mode.of.DE.analysis = deMode,
        Cells.in.foreground = groupSize[ i[1] ],
        Cells.in.background = groupSize[ i[2] ],
        Total.expressed.genes = exprG,
        Differentially.expressed.genes = deg,
        Pathways.examined = totalPaths,
        Pathways.enriched.in.DEGs = signifPaths
      )
    }, deMode=deMode
  ))
  if(deMode == "1 VS rest"){
    groupData$Cells.in.background <- sum( groupData$Cells.in.foreground ) - groupData$Cells.in.foreground
  } 
  
  rownames( groupData ) <- NULL
  return( groupData )
}

paths_likertPlots <- function( DEtabs, DEparams, bg_genes, pathwayEnrichment, 
                               plotComp = NULL, plotPath = NULL, pVthreshold=0.05 ){
  # Function to generate Likert type plots for DE genes in metabolic pathway gene sets. 
  # This function plot either:
  # * Several gene sets for one DE analysis.
  # * Several DE analysis for one gene set.
  ## ins:
  # - DEtabs:
  # - DEparams:
  # - bg_genes:
  # - pathwayEnrichment:
  # - plotComp: comparisons (DE analysis results) to use for the plot.
  # - plotPath: pathway gene sets to use for the plot.
  ## outs:
  # Likert type plot displayig the state ( up-regulated, down-regulated, not DE) of the genes
  # expressed in the subset of cells compared that are part of a given metabolic pathway.
  
  deMode <- DEparams[ "Mode of comparison", ]
  min.ptc <- 100* as.numeric( DEparams["Minimum % of cells expressing a gene",] )
  minLFC <- as.numeric( DEparams[ "Minimum Log2FC", ]) 
  if( plotComp=="All" ){ plotComp <- names( pathwayEnrichment ) }
  
  if( length(plotPath) == 1 & length(plotComp) > 1  ){ # plot one gene set for several DE analysis
    
    pwRes <- do.call( rbind, lapply(
      pathwayEnrichment, function(i){ i[ i$list == plotPath, ] }
    ))
    
    degs <- lapply( 
      seq_along( DEtabs),
      function( i_name ){
        i <- DEtabs[[ i_name ]]
        i <- i[ unlist(strsplit( pwRes[ i_name, ]$DEG, ", ")), ]
        
        c( sum( i$p_val_adj <= pVthreshold & i$avg_log2FC >= minLFC ),
           pwRes[ i_name, ]$list_present - pwRes[ i_name, ]$nDEG,
           sum( i$p_val_adj <= pVthreshold & i$avg_log2FC <= - minLFC ))
      })

    if( deMode %in% c( "Conditional", "1 VS 1" ) ){
      pldf <- data.frame(
        Comparison = names( DEtabs ), 
        "Up_regulated in" = sapply( degs, "[[", 1),
        "Not Differentially expressed" = sapply( degs, "[[", 2),
        "Up_regulated in " = sapply( degs, "[[", 3)
      )
      colnames(pldf) <- gsub( "\\.", " ", colnames(pldf) )
      colnames(pldf) <- gsub( "_", "-", colnames(pldf) )
      
      gr <- strsplit(plotComp, split = " VS ")[[1]]
      colnames(pldf)[2] <- paste0( colnames(pldf)[2], " ", gr[1] )
      colnames(pldf)[4] <- paste0( colnames(pldf)[4], " ", gr[2] )
      
      HH::likert( Comparison ~ . , data=pldf, ylab=NULL, as.percent=FALSE,
        main = list( paste0( "Differential expression in ", gsub("_", " ", plotPath ) ) ),
        sub= list( paste0( "Only genes expressed in >",min.ptc,"% of cells are counted"))
      )
    } else { # DE in "1 VS rest" mode
      
      pldf <- data.frame(
        Comparison = names( DEtabs ), 
        "Up_regulated in foreground" = sapply( degs, "[[", 1),
        "Not Differentially expressed" = sapply( degs, "[[", 2),
        "Down regulated in foreground" = sapply( degs, "[[", 3)
      )
      colnames(pldf) <- gsub( "\\.", " ", colnames(pldf) )
      colnames(pldf) <- gsub( "_", "-", colnames(pldf) )
      
      gr <- strsplit(plotComp, split = " VS ")[[1]]
      
      colnames(pldf)[2] <- paste0( colnames(pldf)[2], " ", gr[1] )
      colnames(pldf)[4] <- paste0( colnames(pldf)[4], " ", gr[1] )
      
      HH::likert( Comparison ~ . , data=pldf, ylab=NULL, as.percent=FALSE,
        main = list( paste0( "DE in ", gsub("_", " ", plotPath ) ) ),
        sub= list( paste0( "Only genes expressed in >",min.ptc,"% of cells are counted"))
      )
    }
  } else if( length(plotComp) == 1 ){ # plot several gene sets for one DE analysis

    pwRes <- pathwayEnrichment[[ plotComp ]]
    if( is.null( plotPath )){
      pwRes <- pwRes[ 1:10, ] # id paths aren't specified, select top10
      plotPath <- pwRes$list
    } else {
      pwRes <- pwRes[ match( plotPath, pwRes$list ), ]
    }
    DEtab <- DEtabs[[ plotComp ]]

    degs <- lapply( 
      plotPath,
      function( i ){
        tabRow <- which(  pwRes$list == i  )
        setGenes <- unlist( strsplit( pwRes$DEG[ tabRow ], ", " ) )
        DEtab <- DEtab[ rownames(DEtab) %in% setGenes, ]
        
        c( sum( DEtab$p_val_adj <= pVthreshold & DEtab$avg_log2FC >= minLFC ),
           pwRes[ tabRow, ]$list_present - pwRes[ tabRow, ]$nDEG,
           sum( DEtab$p_val_adj <= pVthreshold & DEtab$avg_log2FC <= - minLFC ))
      })

    if( deMode %in% c( "Conditional", "1 VS 1" ) ){
      
      pldf <- data.frame(
        Pathway = paths_wrapLabels( plotPath ), 
        "Up_regulated in" = sapply( degs, "[[", 1),
        "Not Differentially expressed" = sapply( degs, "[[", 2),
        "Up_regulated in " = sapply( degs, "[[", 3)
      )
      colnames(pldf) <- gsub( "\\.", " ", colnames(pldf) )
      colnames(pldf) <- gsub( "_", "-", colnames(pldf) )
      
      gr <- strsplit(plotComp, split = " VS ")[[1]]
      
      colnames(pldf)[2] <- paste0( colnames(pldf)[2]," ", gr[1] )
      colnames(pldf)[4] <- paste0( colnames(pldf)[4], " ",  gr[2] )
      
      HH::likert( Pathway ~ . , data=pldf, ylab=NULL, as.percent=FALSE,
                  main = list( paste0( "Altered pathways in ", gsub("_", " ", plotComp ) ) ),
                  sub= list( paste0( "Only genes expressed in >",min.ptc,"% of cells are counted"))
      )
    } else { # DE in "1 VS rest" mode

      pldf <- data.frame(
        Pathway = paths_wrapLabels( plotPath ), 
        "Up_regulated in " = sapply( degs, "[[", 1),
        "Not Differentially expressed" = sapply( degs, "[[", 2),
        "Down regulated in " = sapply( degs, "[[", 3)
      )
      colnames(pldf) <- gsub( "\\.", " ", colnames(pldf) )
      colnames(pldf) <- gsub( "_", "-", colnames(pldf) )
      
      gr <- strsplit(plotComp, split = " VS ")[[1]]
      
      colnames(pldf)[2] <- paste0( colnames(pldf)[2], " ", gr[1] )
      colnames(pldf)[4] <- paste0( colnames(pldf)[4], " ", gr[1] )
      
      HH::likert( Pathway ~ . , data=pldf, ylab=NULL, as.percent=FALSE,
                  main = list( paste0( "Altered pathways in ", gsub("_", " ", plotComp ) ) ),
                  sub= list( paste0( "Only genes expressed in >",min.ptc,"% of cells are counted"))
      )
    }
  }
}

paths_wrapLabels <- function( x, len=40 ){
  # Function to format pathway names for display in likert plots
  sapply(
    x,
    function( i ){
      if( nchar(i) > 50 ){
        str_wrap( gsub( "_"," ",  i), width=len )
      } else {
        gsub( "_", " ", i)
      }
    }, USE.NAMES = FALSE
  )
  
}






