run_QCanalysis <- function( inputData, cellGroups="None", scoreCellCycle="No", iD_summary ){
  # Function to calculate and plot several cell-level features linked to the 
  # quality of the extraction / sequencing / mapping processes.
  ## ins:
  # - inputData: Seurat object
  # - cellGroups (optional): variable grouping the cells (i.e.: 'cluster' or 'sample')
  # - scoreCellCycle (T/F): wheter to estimate the cell cycle phase in which
  #   each cell is.
  # - iD_summary: data frame summarizing characteristics of the input dataset.
  #   output of load_getInputReport() or Merge_getInputReport()
  ## outs, list of;
  # - rawPlot & filteredPlot: arrangement of plots depicting cell features (read count,
  #   number of genes, etc) before (raw) and after (filtered) applying default filters
  #   to the dataset.
  # - CCscoringPlot (optional): plot of a PCA of cell cycle genes with cells colored
  #   by their assigned cell cycle phase.
  # - reportTable: table summarizing cell features before and after applying default
  #   filters
  # - defaultFilters: table containing the default filters applied.
  # - cellMD: data frame containing all cell meta data.

  
  #### plots and default filtering ####
  
  # Check which species the feature IDs belong to
  
  spec <- unlist( iD_summary[, 1:2] )
  spec[1] <- gsub("\\. ", "", spec[1] )
  avblOrgs <- unique( gsub( "\\.txt", "", gsub( "^[A-Za-z]*_", "", dir("data/QC.Genes/") ) ) )

  
  if( spec[1] %in% avblOrgs ){  # skip if species is "Other" or uniavailable in the data folder
    # calculate mitocondrial / ribosomal / ChrY %
    
    if( spec[ 2 ] != "ENSEMBL Gene" ){  # gene symbols
      GeneList <- gen_loadGeneLists( IDformat = "symbols", section = "Rb.Mt.NRY"  )
    } else { # ENSEMBL gene IDs
      GeneList <- gen_loadGeneLists( IDformat = "symbols", section = "Rb.Mt.NRY"  )
    }
    GeneList <- GeneList[[ spec[ 1 ] ]]
    allFeats <- inputData@assays$RNA@counts@Dimnames[[ 1 ]]
    GeneList <- lapply( GeneList, function(i){ i[i %in% allFeats ] } ) #filter out genes not in the dataset
    
    inputData[["percent.mt"]] <- PercentageFeatureSet(
      inputData, features = GeneList[["mitoG"]], assay = "RNA" )
    inputData[["percent.ribo"]] <- PercentageFeatureSet( 
      inputData, features = GeneList[["riboG"]], assay="RNA" )
    if( spec[ 1 ] != "D. rerio" ){
      inputData[["percent.chrY"]] <- PercentageFeatureSet( 
        inputData, features = GeneList[["NRY"]], assay="RNA" )
    } else {
      inputData[["percent.chrY"]] <- NA
    }
    
    inputData[["percent.mt"]] <- round( inputData[["percent.mt"]], digits = 3 )
    inputData[["percent.ribo"]] <- round( inputData[["percent.ribo"]], digits = 3 )
    inputData[["percent.chrY"]] <- round( inputData[["percent.chrY"]], digits = 3 )
    
  } else {
    inputData[["percent.mt"]] <- rep(0, nrow( inputData@meta.data) )
    inputData[["percent.ribo"]] <- rep(0, nrow( inputData@meta.data) )
    inputData[["percent.chrY"]] <- rep(0, nrow( inputData@meta.data) )
  }
  
  defaultFilters <- list(
    Counts = c(1000, Inf),
    Features = c(500, Inf),
    maxMtRNA = c(-Inf,20),
    maxRiboRNA = c(-Inf,20)
  )
  
  rawPlot <- QC_QCplots( inputData, cellGroups )
  
  rawDataReport <- QC_getInputReport( inputData, iD_summary )
  
  cellSub <- QC_inRange( inputData$nFeature_RNA, defaultFilters[["Features"]] ) &
    QC_inRange( inputData$nCount_RNA, defaultFilters[["Counts"]] ) &
    QC_inRange( inputData$percent.mt, defaultFilters[["maxMtRNA"]] ) &
    QC_inRange( inputData$percent.ribo, defaultFilters[["maxRiboRNA"]] )
  gc()

    
  #### cell cycle scoring #### 
  
  if( scoreCellCycle=="Yes" ){

    if( grepl("^ENS",spec[2] ) ){
      ccGenes <- gen_loadGeneLists("ENSids", "CC" )[[spec[1]]]
    } else {
      ccGenes <- gen_loadGeneLists("symbols", "CC" )[[spec[1]]] }
    
    oldDefAssay <-  DefaultAssay( inputData )
    DefaultAssay( inputData ) <- "RNA"
    
    inputData <- NormalizeData(
      inputData, normalization.method = "LogNormalize", scale.factor = 10000)
    inputData <- FindVariableFeatures(
      inputData, selection.method = "vst", nfeatures = 2000)
    inputData <- ScaleData( inputData, features = rownames( inputData ) )
    inputData <- QC_CCscore( inputData, ccGenes=ccGenes, set.ident = FALSE )
    
    inputData <- RunPCA( 
      inputData, features = c( ccGenes$S, ccGenes$G2M),  npcs = 2)
    
    ccPlot <- DimPlot( inputData, group.by = "Phase", reduction = "pca" ) +
      ggtitle("Cell cycle scoring")
    
    DefaultAssay( inputData ) <- oldDefAssay
    
  } else { ccPlot <- NULL }
  
  rawMD <- inputData@meta.data
  
  #### default filtering ####
  inputData <- subset(
    inputData, cells=colnames(inputData)[cellSub] )
  
  filteredPlot <- QC_QCplots( inputData, cellGroups )
  
  filteredDataReport <- QC_getInputReport( inputData, iD_summary )
  
  reportTab <- rbind( rawDataReport, filteredDataReport )
  rownames(reportTab) <- c("Raw data", "Filtered data" )
  colnames(reportTab) <- gsub("\\.", " ", colnames(reportTab))
  
  defaultFilters <- do.call(rbind, defaultFilters)
  colnames( defaultFilters ) <- c("Minimum value", "Maximum value")
  rownames( defaultFilters ) <- c(
    "Read counts per Cell", "Features per Cell",
    "Maximum mitocondrial RNA (%)", "Maximum ribosomal RNA (%)" )
  
  #### returns ####
  
  return( list(
    rawPlot=rawPlot, filteredPlot=filteredPlot, CCscoringPlot=ccPlot,
    reportTable = reportTab, defaultFilters=defaultFilters,
    cellMD = rawMD
  ) )
}


#### auxiliary functions ####

QC_getInputReport <- function( inputData, iD_summary ){
  # Generate a data frame containing a basic summary about the 
  # dataset loaded.
  ## Inputs:
  # - inputData: input dataset, Seurat Object
  # - iD_summary: data frame summarizing characteristics of the input dataset.
  #   output of load_getInputReport() or Merge_getInputReport()
  ## Outputs: data frame specifying species, gene ID type, etc.
  
  if( is.list( inputData ) ){
    # Multiple samples
    inputDF <- do.call(rbind, lapply( inputData, QC_getInputReport, iD_summary=iD_summary ) )
    inputDF$'Sample name' <- names(inputData)
    return( inputDF )
    
  } else {
    if( DefaultAssay( inputData ) != "RNA" ){
      DefaultAssay( inputData ) <- "RNA"
    }
    # One sample
    spec <- unlist( iD_summary[, 1:2] )
    spec[1] <- gsub("\\. ", "", spec[1] )
    
    inputDF <- data.frame(
      "Species" = spec[1],
      "Feature ID type" = spec[2],
      "nCells" = ncol(inputData),
      "Genes detected" = sum( rowSums( inputData ) > 0 ),
      "Median genes per cell" =  median( colSums( 
        GetAssayData(object = inputData, slot = "counts") > 0  ) ),
      "Median counts per cell" = median( colSums( 
        GetAssayData(object = inputData, slot = "counts") ) )
    )
    return( inputDF )
  }
}

QC_get_legend <- function( myggplot ){
  # Function to extract the legend of an arrangement of ggplots
  ## ins: myggplot, ggplot object or gtable
  tmp <- ggplot_gtable( ggplot_build( myggplot ) )
  leg <- which( sapply( tmp$grobs, function( x ) x$name) == "guide-box" )
  legend <- tmp$grobs[[ leg ]]
  return( legend )
}

QC_QCplots <- function( inputData, cellGroups=NULL ){
  # Function to plot QC features of a Seurat object (#genes, #reads, etc)
  ## ins:
  #  - inputData: Seurat object
  #  - cellGroups (optional): variable grouping the cells (i.e.: 'cluster' or 'sample')
  ## outs:
  #  - gtable, arrangement of ggplots.
  
  nGroups <- nrow(unique(inputData[[ cellGroups ]]))
  
  if( !is.null( cellGroups ) ){
    Plot1 <- c( VlnPlot(
      inputData, ncol = 4, group.by = cellGroups , combine = F, 
      log=TRUE, cols = scales::hue_pal()( nGroups ), 
      features = c("nCount_RNA", "nFeature_RNA" ) ),
      VlnPlot(
        inputData, ncol = 4, group.by = cellGroups , combine = F, 
        cols = scales::hue_pal()( nGroups ),
        features = c("percent.mt", "percent.ribo", "percent.chrY"  )  ) )
    
  } else {
    Plot1 <- c( VlnPlot(
      inputData, ncol = 4, combine = F, log=TRUE, 
      features = c("nCount_RNA", "nFeature_RNA" ) ),
      VlnPlot(
        inputData, ncol = 4, combine = F, 
        features = c("percent.mt", "percent.ribo", "percent.chrY"  )  ) )
  }
  
  # adjust log10 break ticks
  Plot1[1:2] <- lapply( Plot1[1:2], function( i ){
    i <- i + annotation_logticks(scaled = TRUE, sides="l")
    i
  } )
  Plot2 <- FeatureScatter( inputData, "nCount_RNA", "nFeature_RNA", raster = T, shuffle = T  ) + 
    scale_y_log10() + scale_x_log10() + annotation_logticks(scaled = TRUE, sides="lb")
  Plot3 <- FeatureScatter( inputData, "nFeature_RNA", "percent.mt", raster = T, shuffle = T  )  +
    scale_x_log10( )+annotation_logticks(scaled = TRUE, sides="b")
  Plot4 <- FeatureScatter( inputData, "nFeature_RNA", "percent.ribo",raster = T, shuffle = T  )  +
    scale_x_log10( )+annotation_logticks(scaled = TRUE, sides="b")
  
  if( ncol(inputData) > 5000 ){
    Plot1 <- lapply( Plot1, function(i){
      i$layers <- i$layers[ -2 ]
      i$layers <- i$layers[ -2 ]
      i <- i + annotation_logticks(scaled = TRUE, sides="l")
      i
    })
    Plot2$layers[[1]]$aes_params$alpha <- 0.3
    Plot2$layers[[1]]$aes_params$size <- 0.5
    Plot3$layers[[1]]$aes_params$alpha <- 0.3
    Plot3$layers[[1]]$aes_params$size <- 0.5
    Plot4$layers[[1]]$aes_params$alpha <- 0.3
    Plot4$layers[[1]]$aes_params$size <- 0.5
  }
  
  leg1 <- QC_get_legend( Plot1[[1]] )
  Plot1 <- lapply( Plot1, function(i){
    i +  theme(legend.position="none")
  })
  leg2 <- QC_get_legend( Plot2 )
  Plot2 <- Plot2 + theme(legend.position="none") +
    ggtitle( "Counts VS Genes" )
  Plot3 <- Plot3 + theme(legend.position="none") +
    ggtitle( "Genes VS % mtRNA" )
  Plot4 <- Plot4 + theme(legend.position="none") +
    ggtitle( "% Ribosomal RNA VS Genes" )
  
  Plot <- arrangeGrob( 
    arrangeGrob( 
      do.call( "arrangeGrob", c( c(Plot1 ), ncol=5) ), 
      arrangeGrob( Plot2, Plot3, Plot4, ncol=3 ),
      ncol=1 ),
    arrangeGrob( leg1, leg2, ncol=1 ),
    ncol=2, widths=c(0.9,0.1) )
  
  return( Plot )
}

QC_inRange <- function( x, limits=NULL, minVal= -Inf, maxVal=Inf ){
  # return TRUE for elements of x are whithin the boundaries of 'limits'
  # Inputs:
  # x: numerical vector
  # limits: numerical vector of length 2 - c( minVal, maxVal )
  # minVal and maxVal are an alternative way to define boundaries.
  
  if( is.null( limits ) ){ limits <- c( minVal, maxVal )}
  x_bool <- (x >= limits[ 1 ] & x <= limits[ 2 ]  )
  x_bool[ is.na( x_bool ) | is.null( x_bool ) ] <- FALSE
  return( x_bool )
}

QC_CCscore <- function( inputData, ccGenes, set.ident = FALSE ){
  # Function to assign a cell cycle phase to each cell in the experiment
  ## ins:
  # - inputData: Seurat object
  # - ccGenes: list of cell cycle genes for S and G2M phases
  # - set.ident: whether to set cell cycle phase as cell identity
  s.genes <- ccGenes$S
  g2m.genes <- ccGenes$G2M
  feats <- inputData@assays$RNA@counts@Dimnames[[1]]
  
  s.genes <- s.genes[ s.genes %in% feats ]
  g2m.genes <- g2m.genes[ g2m.genes %in% feats ]
  
  inputData <- CellCycleScoring(
    inputData, s.features = s.genes,
    g2m.features = g2m.genes, set.ident = set.ident )
  
  return( inputData )
}


