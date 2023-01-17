run_Cluster <- function( inputData=NULL, parList ){
  # Main function to cluster cells in a Seurat object
  ## ins:
  # - inpuData: Seurat Object
  # - parlist: list of clustering parameters, output of Cluster_checkParams().
  ## outs:
  # - Seurat object with extra columns added to the metadata values -cluster asignment
  #     for each resolution in the parList, also clustering uncertainty scores if
  #     the 'measure uncertainty' option is activated.

  # Perform PCA
  if( !'pca' %in% Reductions( inputData ) ){
    inputData <- Cluster_inputPCA( inputData )
  } else {
    if( inputData@reductions[["pca"]]@assay.used != DefaultAssay( inputData )){
      inputData <- Cluster_inputPCA( inputData )
    }
  }
  
  inputData <- FindNeighbors( inputData, dims = 1:parList$nPCs )

  if( parList$uncertainty$measureUncertainty ){
    
    seeds <- runif( parList$uncertainty$iterations, min = 0, max=100 )
    
    if( length( parList$resolution ) > 1 ){
      clUncertainty <- do.call( cbind, lapply( 
        parList$resolution,
        function(i){
          tab <- Cluster_clUncertainty( 
            inputData, clResolution = i, seeds = seeds, ret="array" )
          colnames(tab)<- paste0( colnames(tab), "_res_", i )
          tab
        } ) )
      colnames(clUncertainty)<- gsub( "^[a-z]{3,4}\\.", "", colnames(clUncertainty ), perl = T)
        
      inputData@meta.data <- cbind( inputData@meta.data, clUncertainty )
      
    } else {
      inputData <- Cluster_clUncertainty( 
        inputData, clResolution = parList$resolution, seeds = seeds, ret="seurat" )
    }
  } else {
    inputData <- Cluster_clCells( inputData, parList )
  }
  return( inputData )
}

##### auxiliary functions #####

Cluster_checkParams <- function(
  inputData, nPCs, resolution, makePlot=FALSE, embedding=NULL, clUcertainty=FALSE, clUncIters=NA ){
  # Check the input options for the clustering
  ## ins:
  # - input dataset (Seurat Object)
  # - nPCs: number of principal components from a PCA to use
  # - resolution: clustering resolutions
  # - makePlot [T|F]: whether to generaste 2D reduction plots displaying the clustering
  #     results
  # - embedding ['tsne'|'umap'|'both'|'none']: if makePlot=TRUE, what reduction algorithm to use 
  # - clUncertainty [T|F]: wheter to test clustering uncertainty
  # - clUncIters: if clUncertainty=TRUE, how manky iterations of the clustering algorithms to use
  #     to calculate uncertainty.
  
  ## outs: named list of parameters.

  if ( !is.na(  suppressWarnings( as.integer( nPCs ) ) ) ){
    nPCs <- as.integer( nPCs )
  } else { nPCs <- Cluster_selectNPC( inputData ) }

  # If a precise resolution value has not been defined:
  if( is.na( suppressWarnings( as.numeric( resolution ) ) ) ){ 
    resolution <- trimws( tolower( resolution ) )
    if( ncol( inputData ) > 10000  ){     # Adapt to the size of the dataset.
      resList <- c( low=0.2, mid=0.5, high=1 )
    } else if( ncol(inputData) > 5000 ){           
      resList <- c( low=0.3, mid=0.8, high=1.2 )
    } else{
      resList <- c( low=0.4, mid=1, high=1.5) }
    
    if ( resolution %in% c('low', 'mid', 'high') ){
      resolution <- resList[ resolution ]
    } else {
      resolution <- resList
    }
  } else {
    resolution <- as.numeric( resolution )
  }

  # check if input resolutions have already been used. If they have, results will be overwritten
  
  resInData <-  colnames( inputData@meta.data )
  resInData <- resInData[ grep("Clustering_res_", resInData) ]
  repd <- resolution[ resolution %in% gsub("Clustering_res_", "", resInData) ]
  
  if( is.null(embedding) & makePlot ){
    plotPars <- list( make2Dplots=T, embedding="both")
  } else if(  is.null(embedding) & ! makePlot  ) {
    plotPars <- list( make2Dplots=F, embedding="none")
  } else {
    plotPars <- list( make2Dplots=T, embedding= tolower( embedding ) )
  }

  if( clUcertainty ){
    clUncertainty <- list( measureUncertainty=T, iterations=clUncIters )
  } else {
    clUncertainty <- list( measureUncertainty=F, iterations=NA )
  }
  
  return( list( nPCs=nPCs, resolution=resolution, make2Dplots = plotPars,
                uncertainty=clUncertainty, repeats=repd ) )
}

Cluster_selectNPC <- function( inputData, minDrop = 0.05 ){
  # Function to estimate optimal number of principal components based on
  # standard deviation drop off.
  
  ## ins:
  # - inputData: Seurat ojects
  # - minDrop: minimum decrease in standard deviation.
  
  # outs: recommended number of PCs 
  
  if( !'pca' %in% Reductions( inputData ) ){
    inputData <- Cluster_inputPCA( inputData )
  }
  maxComp <- Cluster_compLimit( inputData, minDrop )
  return( maxComp )
}

Cluster_compLimit <- function( inputData, minDrop=0.05 ){
  # Find optimal number of components of a PCA based on standard deviation
  # drop off.
  pcaStd <- Stdev( inputData[["pca"]] )
  or <- list( 1:length(pcaStd)-1 , 2:length(pcaStd) )
  pcaStd_red <- ( pcaStd[ or[[1]] ] - pcaStd[ or[[2]] ] ) / pcaStd[ or[[1]] ]
  which( zoo::rollsum( pcaStd_red < minDrop, 3, align = "right") > 2 )[1]
}

Cluster_inputPCA <- function( inputData, nFeats=2000 ){
  # Perform Principal component analysis on a Seurat object dataset.
  ## ins: 
  # - inputData: seurat object
  # - nFeats: number of variable features to use -if assay is SCT,set to 3K-.
  ## outs:
  # - seurat object with added PCA in the reductions slot.
  
  if( DefaultAssay(inputData)=="SCT" ){ nFeats <- 3000 }
    
  if( length( VariableFeatures( inputData ) ) == 0 ){
    inputData <- FindVariableFeatures(
      inputData, selection.method = "vst", nfeatures = nFeats )
    inputData <- ScaleData( inputData )
  }
  
  inputData <- RunPCA( inputData, features = VariableFeatures( inputData ) )
  return( inputData )
}

Cluster_clCells <- function( inputData, clusterParams ){
  # Function to cluster cells according to the variables nPC and 
  # resolution in 'clusterParams'.
  ## ins:
  # inputData: Seurat object
  # clusterParams: list of clustering parameters, output of Cluster_checkParams().
  
  ## outs: 
  # A seurat object with an extra metadata column containing the cluster 
  # assigned to each cell.
  
  if( length( clusterParams$resolution ) >1 ){
    for(i in seq_along( clusterParams$resolution ) ){
      res <- clusterParams$resolution[ i ]
      n <- ncol( inputData@meta.data )
      inputData <- FindClusters(
        inputData, 
        resolution = res,
        group.singletons = T)
      colnames( inputData@meta.data )[ n+1 ] <- paste0("Clustering_res_", res )
      inputData@meta.data <- inputData@meta.data[ 
        , colnames(inputData@meta.data) != "seurat_clusters" ]
    }
  } else {
    res <- clusterParams$resolution
    n <- ncol( inputData@meta.data )
    inputData <- FindClusters(
      inputData, 
      resolution = res,
      group.singletons = T)
    colnames( inputData@meta.data )[ n+1 ] <- paste0("Clustering_res_", res )
    inputData@meta.data <- inputData@meta.data[ 
      , colnames(inputData@meta.data) != "seurat_clusters" ]
  }
  
  return( inputData )
}

Cluster_plot <- function( inputData, clusterParams ){
  # Generate tSNE/UMAP plots highlighting cell clusters.
  ##ins:
  # - inputData: seurat object,
  # - clusterParams: list of clustering parameters, output of Cluster_checkParams().
  ## outs:
  # Returns a list of ggplot plots
  
  
  # Select dataset embedding
  emb <- list(tsne="tsne", umap="umap", both=c("tsne","umap") )
  
  emb <- emb[[ clusterParams$make2Dplots$embedding ]]
  
  p <- lapply( 
    emb,
    function( i ){
      if( i == "tsne" ){ 
        inputData <- RunTSNE( inputData, dims = 1:clusterParams$nPCs ) 
      } else if( i == "umap" ){ 
        inputData <- RunUMAP( inputData, dims = 1:clusterParams$nPCs )  
      }
      
      p <- lapply(
        seq_along( clusterParams$resolution ),
        function(j){
          res <- clusterParams$resolution[ j ]
          cl <-  paste0("Clustering_res_", res)
          Idents( inputData ) <- cl
          
          p <- DimPlot( inputData, reduction = i, label = TRUE ) + 
            NoLegend() 
          
          if( !is.null( names( clusterParams$resolution ) ) ){
            p$labels$title <- paste0( 
              tools::toTitleCase( names(clusterParams$resolution)[ j ]), " resolution" )            
          } else {
            p$labels$title <- paste0( "Clustering, resolution = ", clusterParams$resolution[ j ] )
          }
          return(p)
        }
      )
      return(p)
    } )
  return( unlist( p, recursive = F ) )
}

Cluster_summary <- function( inputData, clusterParams ){
  # Function to generate data frames and graphs summarizing clustering results and
  # parameters used.
  ## ins:
  # - inputData: seurat object -with the clustering results in the metadata columns
  # - clusterParams: list of clustering parameters, output of Cluster_checkParams().
  
  ## outs: list of:
  # - clusterDistribution: distribution of cells in each cluster
  # - clusterParams: data frame storing th paramenters used
  # - summaryPlots: list of ggplot plots, output of Cluster_summaryPlots().

  clCols <-  paste0( "Clustering_res_", clusterParams$resolution )
  
  clDistr <- lapply(
    clCols,
    function(i){
      data.frame(
        Cluster= levels( inputData@meta.data[, i] ),
        Cells=as.integer(table( inputData@meta.data[, i] ) ) ) 
    }
  )
  names( clDistr ) <- paste0( "Resolution=",  clusterParams$resolution )
  
  summPlots <- Cluster_summaryPlots( inputData, clusterParams, clDistr )
  
  clusterParams <- clusterParams[ names(clusterParams)!= "repeats" ]
  
  if( clusterParams$make2Dplots$embedding == "both" ){
    clusterParams$make2Dplots$embedding <- "tSNE & UMAP"
  }
  clusterParams <- as.data.frame( clusterParams, row.names = NULL )
  colnames(clusterParams)<- c(
    "Number of PCs", "Clustering resolution","Make 2D plots", "Embedding",
    "Measure Uncertainty", "Uncertainty.iterations"
  )
  
  
  return( list( clusterDistribution = clDistr, clusterParams=clusterParams,
                summaryPlots = summPlots) )
}

Cluster_summaryPlots <- function( inputData, clusterParams, clDistr ){
  # Function to generate graphs summarizing clustering results
  ## ins:
  # - inputData: Seurat object
  # - clusterParams: list of clustering parameters, output of Cluster_checkParams().
  # - clDistr: data frame of cluster sizes, part of the output of clusterSummary()
  ## outs:
  # list of two plots, one for cluster sizes, another for uncertainty scores if calculated.
  maxX <- max( as.numeric( do.call(rbind, clDistr )[,"Cluster"] ) )+0.5
  maxY <- max(do.call(rbind, clDistr )[,"Cells"] )+0.5
  
  clDistr_plots <- lapply( seq_along(clDistr),
    function( i, clDistr, maxX, maxY ){
      
      fillCol <- scales::hue_pal()(length(clDistr))[i]
      
      clDistr[[ i ]]$Cluster <- as.numeric(  clDistr[[ i ]]$Cluster )
      ggplot( clDistr[[ i ]] ) + 
        geom_bar( aes( x=Cluster, y=Cells), fill=fillCol,
                  stat = "identity", show.legend=F) +
        theme_minimal() +
        labs(
          x="Cluster", y="Number of cells", title=names(clDistr)[i]  ) +
        xlim(-0.5,maxX) + ylim(0,maxY)
      
    }, clDistr=clDistr, maxX=maxX, maxY=maxY
  )
  clDistr_plots <- do.call( grid.arrange, c(clDistr_plots, ncol=1) )

  if( clusterParams$uncertainty$measureUncertainty ){
    
    uncCols <- paste0( "Uncertainty_res_", clusterParams$resolution )
    clCols <-  paste0( "Clustering_res_", clusterParams$resolution )
    
    uncert_plots <- lapply(
      seq_along( uncCols ),
      function(i){
        p <-RidgePlot( inputData, features = uncCols[i], group.by = clCols[i] ) 
        p$layers[[1]]$show.legend <- F
        p$labels$x <- "Clustering uncertainty"
        p$labels$y <- "Cluster"
        p$labels$title <- paste0("Resolution = ", gsub("Clustering_res_", "", clCols[i] ) )
        return( p )
      }
    )
    uncert_plots <- do.call( grid.arrange, c(uncert_plots, ncol=1))
  } else { uncert_plots <- NULL }
  
  return(list( Cluster_Distribution=clDistr_plots,
               Uncertainty_Distribution=uncert_plots ))
}


Cluster_clUncertainty <- function( inputData, iters=100, clResolution=NULL,
                                   ret="seurat", seeds=NULL ){
  # Determine the uncertainty of  the clustering label assigned to each cell.
  # The determination if based on reapplying the clustering algorithm with a
  # set of random seeds and meassure the proportion of "cluster partners"
  # mantained by each cell along the iterations.
  ## ins:
  #  - inpuData: scRNA-seq dataset, Seurat object.
  #  - iters: number of iterations.
  #  - clResolution: resolution value passed down to the clustering algorithm.
  #  - ret: return format ("seurat" for Seurat objects, otherwise, a dataframe).
  #  - seeds: numeric vector, seeds for the clustering algorithm.
  ## outs:
  #  - Depending on 'ret', either a seurat Object with added meta data colmns,
  #    or a data frame containing reference cluster asignment and uncertainty
  #    for a given resolution.
  
  gc()
  if( is.null(seeds)){
    seeds <- runif( iters, min = 0, max=100 )
  }
  
  ref_clustering <- FindClusters(
    inputData, resolution = clResolution, 
    group.singletons = T, verbose = F 
  )@meta.data$seurat_clusters
  
  cluster_table <- do.call(cbind, lapply(
    seeds,
    function(i){
      FindClusters(inputData, resolution = clResolution, 
                   group.singletons = T, random.seed=i, verbose = F 
      )@meta.data$seurat_clusters
    }
  ) )
  
  rownames( cluster_table ) <- colnames(inputData)
  colnames( cluster_table ) <- seeds
  
  coCl <- do.call( rbind, lapply(
    1:nrow(cluster_table),
    function(i){
      neighbours <-lapply(
        1:ncol(cluster_table),
        function(j, i ){
          
          cl_i <- cluster_table[ i, j ]
          n <- rownames( cluster_table)[ cluster_table[, j] == cl_i ] 
          n <- n[ n != rownames( cluster_table)[ i ] ]
          n
          
        }, i = i 
      ) 
      neighbours <- table(unlist( neighbours ))
      mean(neighbours) / ncol(cluster_table  )
    }
  ) )
  
  names(coCl) <- rownames(cluster_table)
  inputData$cl_uncertainty <- 1 - coCl
  
  if(ret=="seurat"){
    inputData[[paste0("Clustering_res_", clResolution)]] <- ref_clustering
    inputData[[paste0("Uncertainty_res_", clResolution)]] <- 1 - coCl
    return( inputData )
  } else {
    
      ret <-  data.frame( 
        Clustering=ref_clustering,
        Uncertainty=inputData$cl_uncertainty  )
      rm(inputData)
      return( ret )
  }
}










