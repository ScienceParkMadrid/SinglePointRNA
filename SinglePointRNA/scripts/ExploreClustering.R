run_clExplore <- function( inputData=NULL,
                       range_nPCs=5, PCA_minSDdrop=0.05, jackStraw=FALSE, 
                       clResolution=NULL, final_nPC=NULL, 
                       testUncertainty=FALSE, uncIters=25,
                       plotExamples=FALSE, interactivePlots=FALSE, section=NULL
                       ){
  # Main function to explore clustering options and variability
  ## ins:
  # input dataset (Seurat Object)
  
  # test_nPCs: whether to test the effect of several numbers of principal 
  #   components in clustering.
  # range_nPCs: depending on the value type:
  #   - Single integer: test the recomended PC number +- range_nPCs (default=5)
  #   - Array of length 2: defines min and max number of PCs to test
  # PCA_minSDdrop: minimum drop in standard deviation when adding new PCs, like
  #   determining the start of the "flattening" of an elbow plot.
  # jackStraw: whether to perform jackstraw analysis on the significance of the
  #   principal components. Only for datasets *NOT* normalized with SCTransform.

  # clResolution: array or resolutions for SNN clustering (usually between 0-2)
  # final_nPC: define a number of principal components to use as input for
  #   clustering and example plots.
  # testUncertainty: whether to test clustering uncertainty for each cell.
  # uncIters: Define number of iterations of the algorithm to test uncertainty. 
  
  # interactivePlots: generate plotly interactive cluster flow plots
  # section: defines the block of parameters analized:
  #   - 1: number of principal components
  #   - 2: clustering resolution
  #   - 3: clustering uncertainty
  
  ## outs:
  # section 1: list of: 
  #   - Recommended number of PCs, standard deviation elbow plot, and 
  #     Jackstraw analysis plot (only for non-SCT normalized datasets)
  #   - If test_nPCs, clustering flow plot and co-clustering plot.
  # section 2: 
  #   - Clustering flow plot and co-clustering plot for the different
  #     values of 'resolution' tested.
  #   - Example 2D reduction plots of the dataset.
  #   - Clustering uncertainty estimation for the values of 'resolution'.


  # This script also needs zoo::rollsum()
  
  if( !'pca' %in% Reductions( inputData ) ){
    inputData <- Cluster_inputPCA( inputData )
  } else {
    if( inputData@reductions[["pca"]]@assay.used != DefaultAssay( inputData )){
      inputData <- Cluster_inputPCA( inputData )
    }
  }
  
  if( section==1 ){ # Principal Components
    
    if(is.null(jackStraw)){jackStraw <- FALSE}
    
    PCA_pcDrop <- clExplore_pcDropOff( 
      inputData, PCA_minSDdrop, elbow = T, jack = jackStraw, margin = range_nPCs )

    if( any( is.null( range_nPCs ))){ # auto PC range
      nPCs <- PCA_pcDrop$maxRecommendedPC + c( -5, 5)
    } else {
      if( length( range_nPCs ) == 1 ){ # range = recommended +- x
        nPCs <-  seq( from=PCA_pcDrop$maxRecommendedPC - range_nPCs, 
                      to=PCA_pcDrop$maxRecommendedPC + range_nPCs )
        
      } else if( length( range_nPCs ) == 2 ){ # range defined in input params
        nPCs <-  range_nPCs
      }
    } 

    nPC <- clExplore_clusterPC( inputData, nPCs, clResolution, interactivePlots )
    
    return( list( PCA_pcDrop = PCA_pcDrop,  nPC = nPC ) )
  }
  
  
  if( section==2 ){ # SNN clustering resolution
      

    if( is.null( final_nPC ) ){ # if nPC hasn't been defined
      final_nPC <- clExplore_pcDropOff( inputData, elbow = F, jack = F )$maxRecommendedPC
    }

    if( is.null( clResolution ) ){          # if resolution values aren't provided, select them
      if( ncol( inputData ) > 10000  ){     # by the size of the dataset.
        clResolution <- c( 0.2, 0.5, 1 )
      } else if( ncol(inputData) > 5000 ){           
        clResolution <- c( 0.3, 0.8, 1.2 )
      } else{
        clResolution <- c(0.4, 1, 1.5) }
    }

    inputData <- FindNeighbors( inputData, dims = 1:final_nPC )
    resolution <- clExplore_clusterRes( inputData, clResolution, interactivePlots )
    
    # make some example plots to visualize the effects of resolution variation
    if( plotExamples ){
      if( length( clResolution ) >3 ){
        exPlotRes <- c(min( clResolution ), round(mean( clResolution ), 1), max(clResolution) )
      } else { exPlotRes <- clResolution }
      
      exPlots <- clExplore_examplePlots( inputData, final_nPC, "umap", exPlotRes )
    } else { exPlots <- NULL }

    if( testUncertainty ){
      seeds <- runif( uncIters, min = 0, max=100 )
      
      ##### reducing workers to not exceed memory limits if the datasets is somewhat large:
      if(ncol(inputData)>5000){ nWork <- 2 } else { nWork = 4}
      
      clUncertainty <- bplapply( 
        clResolution,
        function(i){

          clExplore_clUncertainty( inputData, clResolution = i, seeds = seeds, ret="array" )
        }, BPPARAM = MulticoreParam( workers= nWork )
      )
      pl_clUncertainty <- ggplot(  ) +
        geom_violin( aes(
          y=unlist( lapply( clUncertainty, "[[", "Uncertainty" )),
          x= rep( clResolution, each=ncol( inputData )  )  ,
          fill=factor( rep(clResolution, each=ncol( inputData )   ))
        ), scale = "width", show.legend = F,  ) +
        geom_vline( xintercept = clResolution ) + 
        geom_jitter( aes(
          x = rep( clResolution, sapply( clUncertainty, function(i){length(unique(i$refCluster))} )  ) ,
          y = unlist( sapply( clUncertainty, function(i){ by( i$uncertainty, i$refCluster, mean ) } ) ),
          size= unlist( sapply( clUncertainty, function(i){ table(i$refCluster) } ))
        ), width=0.01, alpha=0.4, show.legend = F
        )+
        theme_minimal() +
        labs( x="SNN clustering resolution", y = "Uncertainty" ) +
        ylim( 0,1 )
      
      uncertTab <- do.call(cbind, clUncertainty)
      colnames( uncertTab ) <- paste0( "res_", rep(clResolution, each=2), "_", colnames( uncertTab ) )
    } else { 
      pl_clUncertainty  <- NULL
      uncertTab <- NULL}
    
    return( list( Resolution=resolution, examplePlots=exPlots, 
                  Uncertainty=pl_clUncertainty, uncertTab = uncertTab ) )
  }

}

##### auxiliary functions #####

clExplore_pcDropOff <- function( inputData, minDrop = 0.05, elbow=TRUE, jack=FALSE, margin=3 ){
  # Determine optimal number of principal components for a PCA:
  ## ins:
  # - inputData: Seurat object
  # - minDrop: minimum decrease in standard deviation for principal components.
  # - elbow [T|F] whether to make an elbow plot
  # - jack [T|F] whter to make and plot a jackstraw analysis of the principal components
  #     (non-SCT only)
  # - margin: number of extra PCs to test and plot
  
  ## outs: list of:
  # - maxRecommendedPC: integer, number of recomended PCs
  # - plots: list of elbow plot and jackstraw plot.
  
  maxComp <- clExplore_compLimit( inputData, minDrop )
  
  if( elbow ){ 
    pl1 <- ElbowPlot( inputData, ndims = maxComp + margin  )
    plotList <- list(ElbowPlot=pl1)
  } else {plotList <- list() }
  
  if( jack==TRUE & "SCT" %in% names( inputData@assays ) ){ # this should be blocked at the UI side
    jack=FALSE                                             # but just in case.
    showNotification("Jackstraw analysis is not available for SCT normalized datasets",
                     duration = 20, type="warning", id="exCluster_SCTjack" )    
  }
  
  if( jack ){

    inputData <- JackStraw( inputData, num.replicate = 100, dims = maxComp + margin)
    inputData <- ScoreJackStraw(inputData, dims = 1:(maxComp+margin) )
    
    pl2 <- JackStrawPlot( inputData, dims = 1:(maxComp+margin) )
    plotList <- c( plotList, list( JackStraw=pl2 )  )
  }
  
  return( list( maxRecommendedPC = maxComp, plots = plotList ) )
}

clExplore_compLimit <- function( inputData, minDrop=0.05 ){
  # Find optimal number of components of a PCA based on standard deviation
  # drop off.
  pcaStd <- Stdev( inputData[["pca"]] )
  or <- list( 1:length(pcaStd)-1 , 2:length(pcaStd) )
  pcaStd_red <- ( pcaStd[ or[[1]] ] - pcaStd[ or[[2]] ] ) / pcaStd[ or[[1]] ]
  which( zoo::rollsum( pcaStd_red < minDrop, 3, align = "right") > 2 )[1]
}

clExplore_clusterPC <- function( inputData, nPC, res=0.5, interactive=FALSE ){
  # Funtion to test clustering results according to the number of principal
  # components used as input.
  ## ins:
  # - inputData: seurat Object
  # - nPC: number (or range) of principal components to test
  # - res: provisional clustering resolution - low-ish values are recommended -
  # - interactive [T|F]: whether to make interactive (plotly) or static (ggplot)
  #     plots for the cluster flow.
  
  ## outs: list of:
  # - clusterFlowPlot: plotly or ggplot plot depicting the changes in clustering 
  #     along the PC range
  # - coClustering: ggplot box plot depicting the co-clustering proportion mantained
  #     between a nPC value and the former - i.e., how similar is a clustering using
  #     x PCs to x-1 PCs -.
  # - clusteringTab: data frame storing cluster assignments for each cell at each value
  #     of nPC.
  
  
  if( length(nPC)==2 ){
    pc <- seq( from=nPC[1], to=nPC[2])
  } else { pc <- nPC }
  pc <- pc[ pc <= ncol( inputData@reductions$pca@cell.embeddings ) ]

  cluster_flow <- do.call(cbind, lapply(
    pc,
    function(i, inputData, res){

      tmp <- FindNeighbors( inputData, dims = 1:i, verbose = F )
      tmp <- FindClusters( tmp, resolution = res , group.singletons = T )
      ret <- as.numeric( as.character( tmp$seurat_clusters ) )
      rm(tmp, inputData)
      gc()
      ret
    }, inputData=inputData, res=res
  ) )
  

  rownames( cluster_flow ) <- names(inputData$orig.ident )
  colnames( cluster_flow ) <- pc
  
  if( interactive ){
    
    clustPlot <- do.call( rbind, lapply(
      c( 2:ncol(cluster_flow) ),
      function(i){
        
        fromC <- unique( cluster_flow[, i-1 ] )
        toC <- lapply(
          fromC,
          function(j){
            table( cluster_flow[, i ][ cluster_flow[, i-1 ] == j ])
          }
        )
        fromC <-  rep( fromC, sapply(toC, length) )
        toC <- unlist(toC)
        
        return( data.frame(
          Round=as.numeric(colnames(cluster_flow)[i]),
          nPCs= colnames(cluster_flow)[i],
          fromCluster= fromC,
          toCluster= names( toC ),
          toClusterVal= as.integer( toC )
        ) )
      }
    )  )
    
    nodes <- data.frame(
      label=gtools::mixedsort( unique( c(
        paste0( clustPlot$Round -1, "_", clustPlot$fromCluster) ,
        paste0( clustPlot$Round, "_", clustPlot$toCluster) ) ) )
    )
    
    nodes$color <- scales::hue_pal()( ncol(cluster_flow)  )[ factor(gsub("_.*", "", nodes$label ) ) ]
    nodeSource <- match( paste0( clustPlot$Round -1, "_", clustPlot$fromCluster), nodes$label)
    nodeTarget <- match( paste0( clustPlot$Round , "_", clustPlot$toCluster), nodes$label)
    
    ggpl <- plot_ly(
      type = "sankey",
      arrangement = "snap",
      node = list(
        label = nodes$label,
        color = nodes$color,
        x= nodes$X,
        y= nodes$Y,
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      link = list(
        source = nodeSource-1,
        target = nodeTarget-1,
        value =  clustPlot$toClusterVal
      )
    )
    
  } else {
    
    clustPlot <- data.frame(
      table( apply( cluster_flow, 1, paste0, collapse=".") )
    )
    clustPlot <- data.frame(
      Ident= rep(clustPlot$Var1, each=ncol(cluster_flow) ),
      round= rep( 
        factor( 1:ncol(cluster_flow), levels=1:ncol(cluster_flow)),
        nrow(clustPlot) ),
      nPC =  rep( colnames(cluster_flow), nrow(clustPlot) ),
      Cluster = unlist( strsplit( as.character( clustPlot$Var1), "\\." ) ),
      Freq = rep(clustPlot$Freq, each=ncol(cluster_flow) )
    )
    
    ggpl <- ggplot( clustPlot,
                    aes(x = as.numeric(nPC), stratum = Cluster, alluvium = Ident,
                        y = Freq,
                        fill = Cluster, label = Cluster)) +
      geom_flow() +
      geom_stratum(alpha = .5) +
      geom_text(stat = "stratum", size = 3) +
      theme(legend.position = "none") +
      labs(
        x="Number of Principal Components", y=""
      )
    
  }
  
  
  coCl <- clExplore_coCluster( cluster_flow, variable = "Principal Components" )
  
  
  return( list( clusterFlowPlot=ggpl, coClustering=coCl, clusteringTab = as.data.frame(cluster_flow)  ) )
}

clExplore_clusterRes <- function( inputData, res=seq(from=0, to=2, by=0.1), interactive=FALSE ){
  # Funtion to test clustering results according to the value of clustering resolution
  ## ins:
  # - inputData: Seurat object.
  # - res: numeric vector, resolution values to test
  # - interactive [T|F]: whether to make interactive (plotly) or static (ggplot)
  #     plots for the cluster flow.
  ## outs:
  # - clusterFlowPlot: plotly or ggplot plot depicting the changes in clustering 
  #     along the PC range
  # - coClustering: ggplot bar plot depicting the co-clustering proportion mantained
  #     between a nPC value and the former - i.e., how similar is a clustering using
  #     x PCs to x-1 PCs -.
  # - clusteringTab: data frame storing cluster assignments for each cell at each value
  #     of nPC.
  
  
  cluster_flow <- do.call(cbind, lapply(
    res,
    function(i, inputData){
      as.numeric( as.character( FindClusters(
        inputData, resolution = i , group.singletons = T, verbose = F)$seurat_clusters ) )
    }, inputData=inputData
  ) )
  rownames( cluster_flow ) <- names(inputData$orig.ident)
  colnames( cluster_flow ) <- res
  
  if( interactive ){
    
    clustPlot <- do.call( rbind, lapply(
      c( 2:ncol(cluster_flow) ),
      function(i){
        
        fromC <- unique( cluster_flow[, i-1 ] )
        toC <- lapply(
          fromC,
          function(j){
            table( cluster_flow[, i ][ cluster_flow[, i-1 ] == j ])
          }
        )
        fromC <-  rep( fromC, sapply(toC, length) )
        toC <- unlist(toC)
        
        return( data.frame(
          Round=i-1,
          Granularity= colnames(cluster_flow)[i],
          fromCluster= fromC,
          toCluster= names( toC ),
          toClusterVal= as.integer( toC )
        ) )
      }
    )  )
    
    
    nodes <- data.frame(
      label=gtools::mixedsort( unique( c(
        paste0( res[clustPlot$Round ], "_", clustPlot$fromCluster) ,
        paste0( res[clustPlot$Round+1], "_", clustPlot$toCluster) ) ) )
    )
    nodes$x <- as.numeric( gsub("_.*", "", nodes$label ) ) 
    nodes$y <- as.numeric( gsub(".*_", "", nodes$label ) )
    nodes$x <- nodes$x / max( nodes$x )
    nodes$y <- nodes$y / max( nodes$y )
    nodes$color <- scales::hue_pal()( length(res) )[ rep(
      1:length(res), as.numeric( table( nodes$x) ) ) ]
    
    links <- data.frame(
      from=paste0( res[clustPlot$Round ], "_", clustPlot$fromCluster),
      to=paste0( res[clustPlot$Round +1 ], "_", clustPlot$toCluster)
    )
    links$from_N <- match( links$from, nodes$label )
    links$to_N <- match( links$to, nodes$label )
    
    ggpl <- plot_ly(
      type = "sankey",
      arrangement = "snap",
      node = list(
        label = nodes$label,
        color = nodes$color,
        #x= nodes$x,
        y= nodes$y,
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      link = list(
        source = links$from_N-1,
        target = links$to_N-1,
        value =  clustPlot$toClusterVal
      )
    )
    
  } else {
    clustPlot <- data.frame(
      table( apply( cluster_flow, 1, paste0, collapse=".") )
    )
    clustPlot$Var1 <- as.character( clustPlot$Var1 )
    
    clustPlot <- data.frame(
      Ident= rep(clustPlot$Var1, each=ncol(cluster_flow) ),
      round= rep( 
        factor( 1:ncol(cluster_flow), levels=1:ncol(cluster_flow)),
        nrow(clustPlot) ),
      Resolution =  rep( colnames(cluster_flow), nrow(clustPlot) ),
      Cluster = factor( as.numeric( unlist( strsplit(  clustPlot$Var1, "\\." ) ) ),
                        levels=0:max( cluster_flow ) ),
      Freq = rep(clustPlot$Freq, each=ncol(cluster_flow) )
    )
    
    ggpl <- ggplot( clustPlot,
      aes(x = as.numeric(Resolution)*10, stratum = Cluster, alluvium = Ident,
          y = Freq,
          fill = Cluster, label = Cluster)) +
      geom_flow() +
      geom_stratum(alpha = .5) +
      geom_text(stat = "stratum", size = 3) +
      theme(legend.position = "none") +
      labs(
        x="Clustering Resolution * 10", y=""
      )
    
  }
  
  
  coCl <- clExplore_coCluster( cluster_flow, variable="Clustering Resolution" )
  
  
  return( list(clusterFlowPlot=ggpl, coClustering=coCl, clusteringTab = as.data.frame(cluster_flow) ) )
  
}

clExplore_coCluster <- function( cluster_flow, variable=NA ){
  # Function to compare clustering results along a range of values for a given parameter
  ## ins:
  # - cluster_flow: data frame storing the clustering assignment of every cell (rows)
  #     at each value of the parameter of interest (columns)
  # - variable: parameter tested (in the app currently, 'Clustering Resolution' or 
  #     'Principal Components')
  
  ## outs: ggplot boxplot
  
  
  coCl <- do.call( cbind, lapply(
    2:ncol(cluster_flow),
    function(i){
      sapply(
        rownames(cluster_flow),
        function(j, i ){
          
          tab_r1 <- cluster_flow[,(i-1)]
          tab_r2 <- cluster_flow[,i]
          
          c_r1 <- cluster_flow[j,(i-1)]
          c_r2 <- cluster_flow[j,i]
          
          sum( tab_r1 == c_r1 & tab_r2 == c_r2 ) / sum( tab_r1 == c_r1 )
          
          
        }, i = i
      )
    }
  ) )
  
  colnames(coCl) <- colnames(cluster_flow)[2:ncol(cluster_flow)]
  
  tmp <- reshape2::melt(coCl)
  
  if( grepl("Principal",variable ) ){
    pl <- ggplot( tmp ) +
      geom_boxplot(
        aes(x=factor(Var2, levels=sort(unique(Var2))), 
            y =value, group=Var2), fill=scales::hue_pal()(1) 
      ) + 
      theme_minimal() +
      labs( x=variable, y="Co-clustering proportion" )
  } else if( grepl("Resolution",variable ) ){
    pl <- ggplot( tmp ) +
      geom_boxplot(
        aes(x=Var2, y =value, group=Var2), fill=scales::hue_pal()(1) 
      ) + 
      theme_minimal() +
      labs( x=variable, y="Co-clustering proportion" )
  }
  
  return( pl )
  
  
}

clExplore_examplePlots <- function( inputData, nPCs, reduction="tsne", res=NULL ){
  # Function to generate 2D reduction plots of the dataset ilustrating the effects of
  # different resolution values on the resulting clustering.
  
  ## ins:
  # - inputData: Seurat object
  # - nPCs: number of principal components to use.
  # - reduction ['tsne'|'umap']: dimensional reduction algorithm to use.
  # - res: clustering resolutions to plot.
  
  ## outs: gtable object with the corresponding plots
  
  if(length( nPCs ) > 1 ){
    nPCs <- seq( from = nPCs[1], to = nPCs[2] )
  }
  
  if( is.null(res) ){
    if( ncol(inputData)<5000 ){
      res <- c( 0.4, 0.8, 1.2 )
    } else{
      res <- c(0.5, 1, 1.5) }
  }
  
  pls <- lapply( 
    nPCs,
    function(i){
      lapply(
        res,
        function(j, i){
          inputData <- FindNeighbors( inputData, dims = 1:i )
          inputData <- FindClusters( inputData, resolution = j , group.singletons = T)
          if( reduction=="tsne"){
            inputData <- RunTSNE( inputData, dims = 1:i )
          } else if( reduction=="umap"){
            inputData <- RunUMAP( inputData, dims = 1:i )
          }
          
          DimPlot( inputData, reduction = reduction,label = T ) + NoLegend() + 
            labs( title=j )
        }, i=i
      )
    })
  
  pgrobs <- lapply(
    pls,
    function(i){
      do.call( "arrangeGrob", c( i, ncol=length(i) ) )
    }
  )
  
  pgrobs_resTtl <- arrangeGrob(
    textGrob( "Clustering Resolution", gp = gpar(fontsize = 13, fontface = 'bold')),
    do.call( "arrangeGrob", c( pgrobs, ncol=1 )),
    ncol=1, heights = c(0.05,0.95)
  )
  
  ttl <-  do.call( "arrangeGrob", c(
    lapply(
      nPCs,
      function(i){
        textGrob( paste0(i," Components"), gp = gpar(fontsize = 13, fontface = 'bold')) 
      }), 
    ncol=1 )  )
  
  pgrobs_resTtl_pcTtl <- arrangeGrob(
    ttl,
    pgrobs_resTtl,
    ncol=2, widths = c(0.1,0.9)
  )
  
  
  return(pgrobs_resTtl_pcTtl)
  
}

clExplore_clUncertainty <- function( inputData, iters=100, clResolution=NULL, ret="seurat", seeds=NULL ){
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

clExplore_inputPCA <- function( inputData, nFeats=2000 ){
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








