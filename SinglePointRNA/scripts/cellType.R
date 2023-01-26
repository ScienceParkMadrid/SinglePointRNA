run_cellType <- function( inputData, mode=1, plotMarkerDE=TRUE, groupVar = NULL,
                          m1_species=NULL, m1_tissue=NULL,
                          m2_markers=NULL, m2_DEparams=NULL, m2_diffExpr=NULL
                          ){
  # Main function to annotate cells with a celltype based on the gene expression of a set
  # of markers. 
  ## ins:
  # - input dataset (Seurat Object)
  # - mode: [1|2]: 1=automatic annotation with scCATCH, 2=Manual annotation
  # - plotMarkerDE [TRUE | FALSE ]: generate and return plots displaying the annotation
  # - m1_species: species to select markers from in the cellMatch database.
  # - m1_tissue: tissue to select markers from in the cellMatch database.
  # - m2_markers: data frame of cell types and their markers
  # - m2_DEparams: parameters used to calculate differential expression
  # - m2_diffExpr: list of differential gene expression result data frames, output of 
  #     run_diffExpr().
  
  ## output, a list with 4 elements:
  # - cell.metaData: dataset metadata after cell type imputation
  # - plots: graphical displays of the annotation
  # - imputation: data frame with information on the cell types assigned
  #     to each group in the dataset
  # - markerInfo: data frame with information on the gene markers used for
  #     annotation
  
  if( mode == 1 ){ # automatic annotation
    
    DefaultAssay( inputData ) <- "RNA"
    inputData <- NormalizeData( 
       inputData, normalization.method = "LogNormalize",
       scale.factor = 10000 )

      cellType_res <- cType_scCATCH( inputData,  groupVar, m1_species, m1_tissue, plotMarkerDE )
    
      gc()
    
    return( cellType_res )

  } else if( mode == 2){ # manual annotation
  
    if( ! is.factor( m2_markers$Cell.Type ) ){
      m2_markers$Cell.Type <- as.factor(m2_markers$Cell.Type )
    }
    
    markerExpr <- cType_getMarkerExpr( inputData, m2_markers, m2_diffExpr )
    
    if( m2_DEparams[ "Mode of comparison", ] =="1 VS rest" ){
      
      ctype_annot <- cType_getCtype( markerExpr, minLFC =  as.numeric( m2_DEparams["Minimum Log2FC",] ) )  
      

      foregrounds <- gsub( " VS.*", "", names( ctype_annot )  )
      
      inputData$CellType <- ctype_annot[
        match( as.character( inputData@meta.data[[ groupVar ]] ), foregrounds ) ]
      
      clusterIdents <-  paste0( foregrounds, "-", ctype_annot )
      
      if( all( !is.na( as.numeric( inputData@meta.data[[ groupVar ]] ) ) ) ){ # if group names are numbers
        clusterIdents_levels <- clusterIdents[ order(
          as.numeric( gsub("-.*", "", clusterIdents ) ) ) ]
      } else {
        clusterIdents_levels <- clusterIdents[ order(clusterIdents ) ]
      }
     
      
      matches <- match( as.character( inputData[[groupVar]][,1] ), foregrounds )
      inputData$clusterIdents <- factor( 
        clusterIdents[ matches ], levels = clusterIdents_levels  )
      
      AnnotSumm <- cType_summary( inputData, markerExpr, m2_markers, m2_DEparams, m2_diffExpr )

      
      if( plotMarkerDE == TRUE ){
        ringPlots <- cType_ringPlot( inputData,
          markerExpr, groupVar, names( ctype_annot ), ctype_annot )
        
        if( any( c( "tsne", "umap") %in% Reductions(inputData)) ){
          red <- c( "tsne", "umap")[ c( "tsne", "umap") %in% Reductions(inputData ) ][1]
          
          ringPlots <- cType_arrangePlots(
            inputData, red, groupVar, foregrounds,
            ringPlots, names(ctype_annot) )
          names(ringPlots) <- names(ctype_annot)
        } else {
          names( ringPlots ) <- names( ctype_annot )
        }
      } else { ringPlots <- NULL }
    } else { ringPlots <- NULL }
    
    return( list( cell.metaData = inputData@meta.data, plots = ringPlots,
                  imputation=AnnotSumm$CellSummary, markerInfo=AnnotSumm$MarkerSummary ) ) 
  }
  
  
}



##### auxiliary functions #####
cType_checkGenes <- function( inputData, manual_markers=NULL, auto_markers=NULL, iDsum=NULL ){
  # Check if the cell type markers in manual_markers or auto_markers are
  # present in the dataset. Return notifications if there are missing genes.
  
  ## ins:
  # - inputData: seurat object
  # - manual_markers: data frame of cell types an associated genes
  # - auto_markers: name of the tissue in cellmatch database
  # - iDsum: data frame, summary od dataset characteristics, with a variable named "Species"
  
  ## outs:
  # Shiny pop up notification if:
  #   - there are less than 10 markers in the dataset in auto mode.
  #   - any markers introduced are missing in the dataset.
  
  if(  !is.null( auto_markers )){ # scCATCH
    species <- grep(substr( iDsum$Species,1,1), c("Human","Mouse"), value = T) 
    
    cellmatch <- scCATCH::cellmatch
    cellmatch <- cellmatch[
      cellmatch$tissue == auto_markers & cellmatch$species == species, ]
    
    g <- inputData@assays$RNA@counts@Dimnames[[1]]
    
    if( sum( g %in% cellmatch$gene ) < 10 ){
      showNotification("Less than 10 marker genes from this reference are expressed in the dataset",
                       duration = 1000, type="error", id="cType_gNoti" )   
      return("scCATCH_lowMarker")
    } else { 
      removeNotification( id="cType_gNoti")
    } 
    
  } else if( !is.null( manual_markers )){ # manual annotation
    
    g <- inputData@assays$RNA@counts@Dimnames[[1]]
    allMarks <-  unique( unlist( manual_markers$Markers ) )
    
    marks_miss <- allMarks[ ! allMarks %in% g ]
    
    if(length( marks_miss ) > 0 ){
      showNotification(
        HTML(paste0( "The following genes are not found in the dataset: <br>",
                paste0( marks_miss, collapse = ", ") ) ),
        duration = 20, type="warning", id="cType_gNoti" )    
    } else { 
      removeNotification( id="cType_gNoti")
    } 
    
    
  }
}

cType_getMarkerExpr <- function( inputData, m2_markers, m2_diffExpr ){
  # Generates a list of data frames with the expression of the cell type markers in each
  # cluster.
  
  ## ins:
  # - inputData: seurat object
  # - m2_markers: data frame of cell types an associated genes,
  # - m2_diffExpr: list of differential gene expression results, output of run_diffExpr().
  
  ## outs:
  # - markerExpr: list of dataframes containing the differential expression results
  #    of marker genes in each cluster.
  
  
  markerExpr <- lapply( 
    1:length( m2_diffExpr ),
    function(i){
      tab <- m2_diffExpr[[i]]

      data.frame(
        Group= names(m2_diffExpr)[i],
        Gene= factor( unlist(m2_markers$Markers) ), 
        cellType = rep( m2_markers$Cell.Type, sapply(m2_markers$Markers, length) ) ,
        LFC = tab$avg_log2FC[match( unlist(m2_markers$Markers), rownames(tab) )],
        log10.pval = -log10( tab$p_val_adj[match( unlist(m2_markers$Markers), rownames(tab) )] )
      )
    })
  names( markerExpr ) <- names(m2_diffExpr)
  return( markerExpr)
}

cType_getCtype <- function( markerExpr, minLFC ){
  # Funtion to determine the putative cell type of a group of cells by the 
  # differential expression of selected markers.
  
  ## ins: 
  # markerExpr: data frame, output from cType_getMarkerExpr(), contains 
  #   gene ID, celltype, LFC, and adjusted p-values for selected markers
  # minLFC: minimum absolute LFC threshold to consider differential expression
  #   relevant ( | LFC | > 0.5 by default)
  
  ## outs: 
  # Name of the cell type whose markers are significantly up-regulated and
  # have higher average LFC. If no markers are significantly up, 'Undetermined'
  # is returned
  
  
  unlist(sapply(
    markerExpr, 
    function(i){

        if( sum(i$log10.pval > 1.3, na.rm = T )>=1 ){
        
        i$LFC[is.na( i$LFC )] <- 0
        
        cType_means <- unlist( sapply(
          levels(i$cellType),
          function(i, j){
            
            i_ <- i[ i$cellType == j, ]

            if( sum( i_$log10.pval > 1.3, na.rm = T ) > 0 ){             # if any marker genes are significantly DE:
              sum( i_$LFC[ i_$log10.pval > 1.3 ], na.rm = T ) / nrow(i_) # score= sum(LFC of DE markers)/total number of markers
            } else if( sum( i_$log10.pval > 1.3, na.rm = T ) == 1 ) {
              i_$LFC[ i_$log10.pval > 1.3 & !is.na( i_$log10.pval )] / nrow(i_)  
            } else { NA }
          }, i = i ) )
        if( any( cType_means > minLFC, na.rm = T ) ){
          
          if( sum(cType_means == max( cType_means, na.rm = T ), na.rm = T)>1 ){ # ties
            paste0( levels(i$cellType)[
              which( cType_means == max( cType_means, na.rm = T )  ) ], collapse = "/")
          } else {
            levels(i$cellType)[
              which( cType_means == max( cType_means, na.rm = T )  ) ]
          }
        } else{
          "Undetermined"
        }
      } else { "Undetermined"} 
    }) )
}

cType_updInput <- function( newCT, newMark, currentDF=NULL ){
  # Function to update the inpur marker data frame. Each time a new pair of cell type
  # and gene markers are introduced, it checks if the cell type is already there
  # and replaces the associated markers if appropiate, or appends the ney cell type
  # to the data frame.
  
  ## ins:
  # - newCT: character, cell type name.
  # - newMark: character, gene IDs associated to 'newCT', separated by commas.
  # - currentDF: current marker data frame.
  
  ## outs:
  # dataframe with a new row added or the markers in an existing row replaced.
  
  
  if( is.null( currentDF ) ){
    
    currentDF <-  data.frame( Cell.Type = trimws( newCT ))
    
    currentDF$Markers = list( cType_splitMark( newMark ) )
    
    return( currentDF )
  } else {

    if( newCT %in% currentDF$Cell.Type ){ # overwrite an existing marker list
      
      
      currentDF$Markers[ currentDF$Cell.Type == newCT  ] <- list( cType_splitMark( newMark ) )
    } else {
      
      extraDF <-  data.frame( Cell.Type = trimws( newCT ))
      extraDF$Markers = list( cType_splitMark( newMark ) )
      
      currentDF <- rbind( currentDF, extraDF )
    }
    return( currentDF )
    
  }
}

cType_splitMark <- function( newMark ){
  # Function to split a vector of gene IDs separated by commas.
  newMark <- strsplit( newMark, "," )
  newMark <- unlist( lapply( newMark, trimws ) )
  return(newMark)
  
}

cType_displayMarkers <- function( currentDF ){
  # Function to format the reference marker data frame for GUI display
  ## ins:
  # - currentDF: dataframe with two variables:
  #    - Cell type: cell type name
  #    - Markers: list of markers
  
  ## outs
  # data frame where the 'Markers' variable is switched from list to
  # character vector.
  
  currentDF$Markers <- sapply(
    currentDF$Markers, function(i){
      paste0( i, collapse = ", ")
    }
  )
  return( currentDF )
  
  
}

cType_ringPlot <- function(inputData, markerExpr, groupVar, titles, subtitles ){
  # Function to plot the differential expression values of given markers in a cell group.
  ## ins:
  # iputData: Seurat object
  # markerExpr: data frame, output from cType_getMarkerExpr(), contains 
  #   gene ID, celltype, LFC, and adjusted p-values for selected markers
  # groupVar: grouping variable -usually assigned cluster-.
  # titles: char vector, names of the cell groups
  # subtitles: char vector, putative cell types of the cell groups.
  
  ## outs:
  # list of plots of LFC values of marker genes by cell type in a radial arrangement, 
  # with genes significantly DE in bright colors and genes not significantly DE in 
  # muted colors.
  
  markerExpr <- lapply( 
    markerExpr,
    function(i){
      
      i$LFC[is.na( i$LFC )] <- 0
      i$log10.pval[ is.na( i$log10.pval ) ] <- 0
      
      
      empty_bar <- 2
      to_add <- data.frame( matrix(NA, empty_bar*nlevels( i$cellType ), ncol( i )) )
      colnames(to_add) <- colnames( i )
      to_add$cellType <- rep(levels( i$cellType ), each=empty_bar )
      i <- rbind(i, to_add)
      i <- i %>% arrange( cellType )
      i$id <- seq(1, nrow( i ))
      
      i$signif <- factor( i$log10.pval>3, levels= c(F,T), labels=c("No","Yes") )
      
      return( i )
    })
  
  plotList<- lapply(
    1:length( unique( inputData[[ groupVar ]][,1]) ),
    function(i){
      
      tab <- markerExpr[[i]]
      
      # separator between cell types
      empty_bar <- 2
      
      # Gene labels 
      label_data <- tab
      number_of_bar <- nrow(label_data)
      angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
      label_data$hjust <- ifelse( angle < -90, 1, 0)
      label_data$angle <- ifelse(angle < -90, angle+180, angle)
      label_data$y <- label_data$LFC + 0.5
      label_data$y[ label_data$y < 0 ] <- 0.5
      
      # cell type labels
      base_data <- tab %>% 
        group_by( cellType ) %>% 
        summarize( start=min(id), end=max(id) - empty_bar) %>% 
        rowwise() %>% 
        mutate( title=mean( c( start, end ) ) )
      angle <- 5 - 360 * (base_data$title)/number_of_bar  
      base_data$labelAngle <- ifelse(angle < -90 & angle > -270, angle+180, angle)
      base_data$start <- base_data$start -0.3
      base_data$end <- base_data$end + 0.3
      base_data$labCol <- scales::hue_pal()(nrow(base_data))
      base_data$labY <- 5
      base_data$lineY <- 4.7
      
      
      # grid labels
      grid_scale <- base_data
      grid_scale$end <- grid_scale$end[ c( nrow(grid_scale), 1:nrow(grid_scale)-1)] + 0.5
      grid_scale$start <- grid_scale$start - 0.5
      grid_scale <- grid_scale[-1,]
      
      # plot Y limits for exterior labels
      if( max(tab$LFC, na.rm = T) < 2.5 ){ 
        plotlimY <- c( min(tab$LFC, na.rm = T)-0.5, 5 )
      } else { 
        plotlimY <- c( min(tab$LFC, na.rm = T)-0.5, max(tab$LFC, na.rm = T) * 2 )
        base_data$labY <- max(tab$LFC, na.rm = T) * 2
        base_data$lineY <- max(tab$LFC, na.rm = T) * 1.8
      } 
      
      
      grid_ticks <- c(
        0, min(tab$LFC, na.rm = T), mean(tab$LFC, na.rm = T), max(tab$LFC, na.rm = T)
      )
      
      ggplot(tab, aes(x=as.factor(id), y=LFC, fill= cellType)) + 
        
        geom_bar(aes(alpha=signif ), stat="identity") +
        geom_segment(data=grid_scale, aes(x = end, y = grid_ticks[1], xend = start, yend = grid_ticks[1]),
                     colour = "black", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_scale, aes(x = end, y = grid_ticks[2], xend = start, yend = grid_ticks[2]),
                     colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_scale, aes(x = end, y = grid_ticks[3], xend = start, yend = grid_ticks[3]),
                     colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_segment(data=grid_scale, aes(x = end, y = grid_ticks[4], xend = start, yend = grid_ticks[4]),
                     colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
        geom_bar(aes(alpha=signif ), stat="identity") +
        annotate("text", x = rep(max(tab$id),4), y = grid_ticks, label = round(grid_ticks,1) ,
                 color="grey40", size=3 , angle=0, fontface="bold", hjust=1) +
        # cell type labels
        geom_segment( data=base_data,
                      aes(x = start, y = lineY, xend = end, yend = lineY ),
                      colour = base_data$labCol,
                      alpha=0.8, size=2 , inherit.aes = FALSE ) +
        geom_text( data=base_data,
                   aes(x = title, y = labY, label=cellType),
                   colour = "black", alpha=0.8, size=4, fontface="bold",
                   angle=base_data$labelAngle,inherit.aes = FALSE) +
        # gene labels
        geom_text(data=label_data, 
                  aes(x=id, y=y, label=Gene, hjust=hjust),
                  color="black", fontface="bold",alpha=0.6, 
                  size=3, angle= label_data$angle, inherit.aes = FALSE ) +
        
        geom_bar(aes(alpha=signif ), stat="identity") +
        ylim( plotlimY[1], plotlimY[2] ) +
        coord_polar() + theme_minimal() +
        theme(
          legend.position = "none",
          axis.text = element_blank(),
          axis.title = element_blank() )  +
        labs( title=paste0( "Diff. expression (", groupVar, ") - ",titles[ i ]),
              subtitle = paste0("Cell type - ", subtitles[ i ]) )
    }
  )
  
  return( plotList )
}

cType_arrangePlots <- function( inputData, red, groupVar, foregrounds, ringPlots, titles ){
  # Function to add 2D reduction plots to a list of ring plots if 2D embeddings (UMAP or tSNE)
  # are available in the Seurat object.
  
  ## ins:
  # - inputData: Seurat object.
  # - red ['tsne'|'umap']: reduction name
  # - groupVar: grouping variable -usually clustering-.
  # - foregrounds: char vector, foreground groups for every comparison
  # - ringPlots: list of ggplot plots, output of cType_ringPlot().
  # - titles: char vector, names of the cell groups.
  
  ## outs: list of gtables, each with a 2D embedding plot higlighting each foreground
  # cell group and a ringplot of cell type markers' differential expression.
  
  inputData$groupID <- paste0( inputData[[ groupVar ]][,1], " - ", inputData$CellType  )
  
  Idents( inputData ) <- "groupID"
  
  dimplots <- lapply(
    seq_along( levels( Idents( inputData ) ) ),
    function(i){
      
      highlight <- Idents( inputData )[  inputData[[ groupVar ]][,1] == foregrounds[i] ][1]
      highlight <- which( levels( Idents( inputData )) == highlight )
      
      nGroups <- length( unique( inputData[[ groupVar ]][,1 ] ) )
      colList <- scales::hue_pal(c = 10)(nGroups)
      colList[ highlight ] <-  scales::hue_pal()( nGroups )[ highlight ]
      DimPlot( inputData, reduction = red,  cols = colList, label=T ) + NoLegend()
    }
  )
  
  plotList <- lapply(
    seq_along( levels( Idents( inputData ) ) ),
    function(i){
      ttl <- textGrob( paste0( titles[i] ), gp = gpar(fontsize = 13, fontface = 'bold')) 
      
      arrangeGrob(
        ttl,
        arrangeGrob( dimplots[[i]], ringPlots[[i]], ncol=2 ),
        ncol=1, heights = c(0.1,0.9)
      )
    }
  )
}

cType_summary <- function( inputData, markerExpr, m2_markers, m2_DEparams, m2_diffExpr ){
  # Function to make a brief summary data frame of the cell type imputation results.
  ## ins:
  # - inputData: Seurat object
  # - markerExpr: data frame, output from cType_getMarkerExpr(), contains 
  #   gene ID, celltype, LFC, and adjusted p-values for selected markers.
  # - m2_markers: data frame of cell types an associated genes.
  # - m2_DEparams: parameters used to calculate differential expression.
  # - m2_diffExpr: list of differential gene expression results, output of run_diffExpr().
  
  ## outs:
  # Two data frames, one for cells and another for marker genes:
  #   - CellSummary: Data frame containing cell group names, group sizes, details of the
  #       differential expression analysis used as input, selected markrs that are
  #       differentially expressed in each group, and putative cell type
  #   - MarkerSummary: data frame containing gene ID, cell type they're linked to, and
  #       expression levels in the experiment.
  
  
  # making sure the order of dfs in markerExpr & diffExpr is the same as the groups in
  # the table
  groupVar <- m2_DEparams["Grouping Variable",]
  
  groups <- sort( unique( inputData[[ groupVar ]][,1] ) )
  markerExpr <- markerExpr[ match( groups, gsub(" VS.*","", names( markerExpr ) ) ) ]
  m2_diffExpr <- m2_diffExpr[ match( groups, gsub(" VS.*","", names( m2_diffExpr ) ) ) ]
  
  ct <- unique( inputData[[c( groupVar,"CellType" )]] )
  ct <- ct[ match( groups, ct[, 1] ), ]
  
  minFC <- as.numeric( m2_DEparams["Minimum Log2FC",] )
  
  
  sumCells <- data.frame(
    Group = ct[ , groupVar ],
    Size = as.numeric( table( inputData[[ groupVar ]][,1 ] ) ),
    DE.analysis.mode = m2_DEparams[ "Mode of comparison", ],
    Num.of.DEGs = sapply( m2_diffExpr, function(i){ sum( i$p_val_adj < 0.05, na.rm = TRUE ) } ),
    Num.of.up.genes = sapply( m2_diffExpr, function(i){ 
      sum( i$avg_log2FC> minFC & i$p_val_adj < 0.05, na.rm = TRUE ) } ),
    Num.of.down.genes = sapply( m2_diffExpr, function(i){
      sum( i$avg_log2FC < -minFC & i$p_val_adj < 0.05, na.rm = TRUE ) } )
    
  )
  sumCells$Markers.UpReg <- sapply( markerExpr, function(i){
    g <- i$Gene[ i$LFC> minFC & i$log10.pval>1.3 ]
    paste0( as.character( g[ !is.na(g) ]), collapse = ", ") } )
  sumCells$Markers.DownReg <- sapply( markerExpr, function(i){
    g <- i$Gene[ i$LFC < -minFC & i$log10.pval>1.3 ]
    paste0( as.character( g[ !is.na(g) ]), collapse = ", ") } )
  
  sumCells$Putative.cell.type <- ct$CellType
  
  inputData <- cType_CPMnorm( inputData )
  
  sumMarkers <- data.frame(
    Gene = unlist( m2_markers$Markers ),
    Cell.type = rep( m2_markers$Cell.Type, sapply(m2_markers$Markers, length ) )
  )
  
  sumMarkers$Pct.cells.expressing <- round(digits = 2, sapply(
    sumMarkers$Gene, function(i){
      100*sum( inputData@assays$RNA@counts[ i,  ]  > 0 ) / ncol( inputData ) }  ) )
  sumMarkers$mean.CPM.in.expressing.cells <- round(digits = 2, sapply(
    sumMarkers$Gene, function(i){
      c <- inputData@assays$RNA@data[ i,  ]
      mean( c[c>0]) }  ) )
  
  return( list(
    CellSummary = sumCells, MarkerSummary = sumMarkers
  ))
  
}

cType_CPMnorm <- function( inputData ){
  # Return normalized (CPM) RNA counts as the default assay of inputData.
  ## ins: 
  #   - inputData: scRNA-seq dataset, Seurat object.
  ## outs: 
  #   - Seurat object with DefaultAssay set to "RNA" and counts
  #     normalized to CPM values.
  
  DefaultAssay( inputData ) <- "RNA"
  NormalizeData( 
    inputData,
    normalization.method = "RC",
    scale.factor = 1e6 )
}

cType_scCATCH <- function( inputData, groupVar, species, tissue, plotMarkerDE ){
  # Function to run scCATCH's cell type imputation analysis.
  
  ## ins:
  # - inputData: Seurat object
  # - groupVar: grouping variable -usually assigned cluster-.
  # - species ['Hsapiens'|'Mmusculus']: species of the experiment sample -in the app
  #     is autodetected
  # - tissue: selected tissue name.
  # - plotMarkerDE [ TRUE | FALSE ]: whether to make an annotated 2D plot -if the
  #     UMAP/tSNE reduction is already available.
  
  ## outs: list of
  # - cell.metaData = data frame of cell meta data to update the object outside the
  #     main function.
  # - scCATCH.pars: input parameters used in scCATCH
  # - imputation: data frame with details on asigned cell type and DE markers.
  # - markerInfo: details on the markers used, refrences to literature linking
  #     each marker to a particular cell type.
  # - plots: if plotMarkerDE = TRUE, 2D reduction plot with cell type annotation.
  
  if( any( is.na(  inputData[[ groupVar]][,1] )) ){  # in case there's missing info
    inputData <- subset( inputData, cells= colnames(inputData)[
      !is.na(  inputData[[ groupVar]][,1] )] )
  }
  
  groups <- unique( inputData[[ groupVar]][,1] )
  
  if( species == "H. sapiens" ){ spec <- "Human"
  } else if( species=="M. musculus"){spec <- "Mouse"}
  
  
  scC <- createscCATCH( 
    inputData@assays$RNA@data, cluster = as.character( inputData[[groupVar]][,1] ) )
  scC <- findmarkergene( object = scC, 
    species = spec, marker = cellmatch, tissue = tissue, use_method = 2,
    cell_min_pct = 0.25, logfc = 0.5 )
  scC <- findcelltype( object = scC )
  
  if(nrow(scC@celltype)< length(groups)){ #if any cluster doesn't have any DE cell type makers
    mis <- as.character( groups[ ! groups %in% scC@celltype$cluster ] )
    mis <- as.data.frame( do.call(rbind, lapply(
      mis, function(i){ c( i, "No DE cell type markers found", "Unknown", NA,NA,NA) }
    )))
    colnames(mis) <- colnames(scC@celltype)
    scC@celltype <- rbind( scC@celltype, mis ) 
  }
  
  
  clusterIdents <- paste0(
    scC@celltype$cluster, "-", scC@celltype$cell_type )

  if( all( !is.na( as.numeric( groups ) ) ) ){ # if group names are numbers
    clusterIdents_levels <- clusterIdents[ order(
      as.numeric( gsub("-.*", "", clusterIdents ) ) ) ]
  } else {
    clusterIdents_levels <- clusterIdents[ order(clusterIdents ) ]
  }
  
  matches <- match( as.character( inputData[[groupVar]][,1] ), scC@celltype$cluster)
  
  inputData$CellType <- scC@celltype$cell_type[ matches ]
  
  inputData$CellType_score <- scC@celltype$celltype_score[ matches ]
  
  inputData$clusterIdents <- factor( 
    clusterIdents[ matches ], levels = clusterIdents_levels  )
  
  
  if( any( c( "tsne", "umap") %in% Reductions(inputData)) & plotMarkerDE ){
    red <- c( "tsne", "umap")[ c( "tsne", "umap") %in% Reductions(inputData ) ][1]
    dimplot <- DimPlot( 
      inputData, group.by = "clusterIdents", label = T, reduction = red ) + 
      labs(title="Automatic cell type annotation",
           subtitle = paste0( "scCATCH - ", spec, ", ", tissue, "."))
  } else { dimplot = NULL }
  
  return( list(
    cell.metaData = inputData@meta.data, 
    scCATCH.pars = scC@para, imputation=scC@celltype, markerInfo=scC@marker,
    plots = dimplot
  ))
}





