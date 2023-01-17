run_Visualization <- function( inputData, groupVar=NULL, algorithm="tsne",
                               saveMode=1, paramMode=NULL, seedSet="default",
                               tsne_perplx=30, umap_nNeighbors=30, umap_minDist=0.3
                               ){
  # Main function to generate 2D visualizations using dimensional reduction techniques.
  # Either try several random seeds with a fixed set of parameters, or try a range of
  # parameters with a fixed seed.
  ## ins:
  # - input dataset (Seurat Object)
  # - groupVar: metadata variable to color cells by (i.e.: cluster, sample of origin).
  # - algorithm [ 'tsne'|'umap' ]: dimension reduction algorithm
  # - saveMode [ 1 | 2 ]: if 1, generate and return a list of plots, if 2, calculate 
  #   dimensional reduction and return a Seurat object.
  # - seedSet ['default' | number | 'random' ]: set random seed/seeds used to calculate
  # - paramMode ['default' | 'userSet' | 'size' ] :
  #   tSNE / UMAP embeddings.
  # - tsne_perplx: tSNE parameter ,'Perplexity'
  # - umap_nNeighbors: UMAP parameter, number of neighbors.
  # - umap_minDist: UMAP parameter, minimum distance.
  
  ## outs:
  # - if saveVisuals == T : 
  #   - Seurat object
  # - if saveVisuals == F : 
  #   - list of ggplot plots
  
  # Set params
  if( seedSet=="random" ){ # generate set of random seeds.
    if( algorithm=="tsne" ){seedList <- round( runif( 6, min=0, max=100), 2 ) 
    } else if( algorithm =="umap"){seedList <- round( runif( 6, min=0, max=100) )}
  } else if ( seedSet == "default" ){ # use Seurat's default seed values.
    if( algorithm=="tsne" ){ seedList <- 1 
    } else if( algorithm =="umap"){ seedList <- 42}
  } else if (is.numeric( seedSet ) ){ # user-defined seed.
    seedList <- seedSet }
  
  #print(umap_nNeighbors)
  
  if( paramMode == "size" ){ # Set parameter ranges according to the dataset's size 
    parList <- Visu_setParaRange( ncol(inputData) )
    tsne_perplx <- parList$perplx
    umap_nNeighbors <- parList$nNeighbors
    umap_minDist <- parList$minDist
  } else if ( paramMode == "default" ){ # use Seurat's default parameter values.
    tsne_perplx <- 30
    umap_nNeighbors <- 30
    umap_minDist <- 0.3
  } 
  
  
  # Perform PCA
  if( !'pca' %in% Reductions( inputData ) ){
    inputData <- Cluster_inputPCA( inputData )
  } else {
    if( inputData@reductions[["pca"]]@assay.used != DefaultAssay( inputData )){
      inputData <- Cluster_inputPCA( inputData )
    }
  }
  if( any( grepl("FindNeighbors\\..*\\.pca", names(inputData@commands))) ){ 
    # set the same number od PCs as the clustering algorithm if applicable
    cmm <- grep("FindNeighbors\\..*\\.pca", names(inputData@commands))
    if( length( cmm ) > 1 ){ cmm <- cmm[ length(cmm) ] }
    maxDims <- max( inputData@commands[[ cmm ]]$dims  )
  } else {
    maxDims <- Visu_compLimit( inputData  )
  }
  
  
  
  
  if( saveMode==1 ){
    plotList <- Visu_bagOfPlots(
      inputData, groupVar, dims = 1:maxDims, algorithm, seedList, 
      tsne_perplx, umap_nNeighbors, umap_minDist  )
    return( plotList )
    
  } else if( saveMode == 2 ) {
    inputData <- Visu_reduce(inputData, dims = 1:maxDims, algorithm, seedList, 
                             tsne_perplx, umap_nNeighbors, umap_minDist  )
    return( inputData )
  }
  
}

##### auxiliary functions #####
Visu_inputPCA <- function( inputData, nFeats=2000 ){
  # Function to perform principal component analysis on variable features
  # of a dataset.
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

Visu_compLimit <- function( inputData, minDrop=0.05 ){
  # Find optimal number of components of a PCA based on standard deviation
  # drop off.
  pcaStd <- Stdev( inputData[["pca"]] )
  or <- list( 1:length(pcaStd)-1 , 2:length(pcaStd) )
  pcaStd_red <- ( pcaStd[ or[[1]] ] - pcaStd[ or[[2]] ] ) / pcaStd[ or[[1]] ]
  which( zoo::rollsum( pcaStd_red < minDrop, 3, align = "right") > 2 )[1]
}

Visu_setParaRange <- function( nCells ){
  # Return a list o parameters for tSNE or UMAP dimensional reduction 
  # according to the number of cells in the dataset
  ## ins:
  # - nCells: number of cells in the analyzed dataset
  ## outs: list of parameter for tSNE or UMAP:
  # - perplx: perplexity (tSNE)
  # - minDist: Minimum distante (UMAP)
  # - nNeighbors: Number of neighbors (UMAP)
  
  if( nCells > 10000  ){     
    perplx <- c(25, 40, 50)
    minDist <- c( 0.05, 0.1, 0.25 )
    nNeighbors <- c( 25, 40, 50 )
  } else if( nCells > 5000 ){ 
    perplx <- c( 20, 30, 45 )
    minDist <- c( 0.2, 0.25, 0.35 )
    nNeighbors <- c( 20, 30, 45  )
  } else if( nCells > 1000 ){
    perplx <- c( 15, 25, 35 )
    minDist <- c( 0.25, 0.35, 0.45 )
    nNeighbors <- c( 15, 25, 35  )
  } else {
    perplx <- c( 5, 15, 25  )
    minDist <- c( 0.3, 0.4, 0.5 )
    nNeighbors <- c( 5, 15, 25 )
  }
  return( list(
    perplx = perplx, minDist = minDist, nNeighbors = nNeighbors
  ))
  
}

Visu_getParList <- function( input ){
  # Function to extract the relevant input parameters to generate 2D embeddings
  # of the dataset from the app's input values.
  ## ins:
  # - input: shiny app input storage object.
  ## outs: list of 
  #  groupVar: cell grouping variable,
  # algorithm ['tsne'|'umap']: dimensional reduction algoritm,
  # paraMode ['default'|'size'|'userSet']: criteria to select reduction parameters
  # seedSet ['random'|'default'|'userSet']: criteria to set random seeds for the
  #   tSNE / UMAP algoritms
  # tsne_perplx: perplexity value
  # umap_nNeighbors: number of neighbors value
  # umap_minDist: minimum distance value
  
  
  input <- reactiveValuesToList(input)
  
  if(input$Visu_saveMode ==1 ){
    
    groupVar <- input$Visu1_groupVar
    algorithm <- input$Visu1_algorithm
    paraMode <- input$Visu1_paramMode
    tsne_perplx <- input$Visu1_perplx
    umap_nNeighbors <- input$Visu1_nNeighbors
    umap_minDist <- input$Visu1_minDist
    
    if( input$Visu1_setSeed == "userSet" ){
      seedSet <- input$Visu1_setSeedValue
    } else {
      seedSet <- input$Visu1_setSeed
    }
    
    return( list(
      groupVar=groupVar, algorithm=algorithm, paraMode=paraMode, seedSet=seedSet,
      tsne_perplx=tsne_perplx, umap_nNeighbors=umap_nNeighbors,umap_minDist=umap_minDist
    ))
    
  } else if( input$Visu_saveMode == 2 ){

    algorithm <- input$Visu2_algorithm
    paraMode <- input$Visu2_paramMode
    tsne_perplx <- input$Visu2_perplx
    umap_nNeighbors <- input$Visu2_nNeighbors
    umap_minDist <- input$Visu2_minDist
    
    if( input$Visu2_setSeed == "userSet" ){
      seedSet <- input$Visu2_setSeedValue
    } else {
      seedSet <- input$Visu2_setSeed
    }
    
    return( list(
      groupVar=NULL, algorithm=algorithm, paraMode=paraMode, seedSet=seedSet,
      tsne_perplx=tsne_perplx, umap_nNeighbors=umap_nNeighbors,umap_minDist=umap_minDist
    ))
  }
  
}


Visu_reduce <- function( inputData, dims=NULL, algorithm="tsne", seed=NULL, tsne_perplx=NULL,
                         umap_nNeighbors=NULL, umap_minDist=NULL ){ 
  # Function to run the dimensional reduction functions on the PCA of the dataset.
  ## ins:
  # - inputData: Seurat object
  # - dims: numper of principal components to use for the reduction
  # - algorithm ['tsne'|'umap']: dimensional redution algorithm
  # - tsne_perplx: perplexity value
  # - umap_nNeighbors: number of neighbors value
  # - umap_minDist: minimum distance value
  ## outs:
  # the same inputData object with an extra element in the 'reductions' slot of the Seurat object.
  
  if( algorithm == "tsne" ){
    inputData <- RunTSNE( inputData, dims=dims, seed.use = seed, perplexity = tsne_perplx )
  } else if( algorithm == "umap" ){
    inputData <- RunUMAP( inputData, dims=dims, seed.use = seed, n.neighbors= umap_nNeighbors,
                          min.dist = umap_minDist )
  }
  return(inputData)
}


Visu_bagOfPlots <- function( 
  inputData, groupVar, dims=NULL, algorithm="tsne", seedList=NULL, tsne_perplx=NULL,
  umap_nNeighbors=NULL, umap_minDist=NULL ){
  # Generate several dimensional reduction plots using a set of random seeds.
  ## ins:
  # - input dataset (Seurat Object)
  # - groupVar: metadata variable to color cells by (i.e.: cluster, sample of origin).
  # - dims: principal components to use (i.e. '1:10')
  # - algorithm: dimension
  # - seed: random seed used in tSNE / UMAP
  # - tsne_perplx: tSNE parameter ,'Perplexity'
  # - umap_nNeighbors: UMAP parameter, number of neighbors.
  # - umap_minDist: UMAP parameter, minimum distance.
  ## outs:
  # - list of ggplot plots  
  
  if( length( seedList ) == 1){ # try ranges of params
    
    if( algorithm == "umap"){
      redList <- lapply( 
        umap_nNeighbors,
        function(i){
          ret <- lapply(
            umap_minDist,
            function(j, i){
              inputData <- Visu_reduce( inputData, dims ,"umap", seedList, umap_nNeighbors = i, umap_minDist = j )
              redCoords <- data.frame( inputData@reductions$umap@cell.embeddings )
              pl <- DimPlot( inputData, reduction = "umap", label = T, group.by = groupVar ) + NoLegend() 
              pl$labels$title <- paste0( "UMAP: # neigbors=", i, " ; Min. distance=", j )
              list( plot= pl, redCoords=redCoords )
            }, i=i )
          names(ret) <- paste0("minDist_", umap_minDist )
          ret
        })
      names(redList) <- paste0( "nNeighbors_", umap_nNeighbors )
      
      pgrobs <- lapply( unlist( redList, recursive=FALSE ), "[[", "plot" )
      redCoords <- lapply( unlist( redList, recursive=FALSE ), "[[", "redCoords" )
      pgrobs <- do.call( grid.arrange, c( pgrobs, ncol=3 ) )
    } else if( algorithm == "tsne" ){
      redList <- lapply( 
        tsne_perplx,
        function(i){
          inputData <- Visu_reduce( inputData, dims, "tsne", seedList, tsne_perplx = i )
          redCoords <- data.frame(inputData@reductions$tsne@cell.embeddings)
          pl <- DimPlot( inputData, reduction = "tsne", label = T, group.by = groupVar ) + NoLegend() 
          pl$labels$title <- paste0( "tSNE: Perplexity=", i )
          list( plot= pl, redCoords=redCoords )
        })
      names(redList) <- paste0( "Perplexity_", tsne_perplx )
      pgrobs <- lapply( redList, "[[", "plot" )
      redCoords <- lapply( redList, "[[", "redCoords" )
      pgrobs <- do.call( grid.arrange, c( pgrobs, ncol=3 ) )
    }
    
  } else { # try several random seeds
    
    if( algorithm == "umap"){
      redList <- lapply( 
        seedList,
        function(i){
          inputData <- Visu_reduce( inputData, dims, "umap", i, umap_nNeighbors = umap_nNeighbors,
                                    umap_minDist = umap_minDist )
          redCoords <- data.frame(inputData@reductions$umap@cell.embeddings)
          pl <- DimPlot( inputData, reduction = "umap", label = T, group.by = groupVar ) + NoLegend() 
          pl$labels$title <- paste0( "UMAP: random seed =", i )
          list( plot= pl, redCoords=redCoords )
        })
      
    } else if( algorithm == "tsne" ){
      redList <- lapply( 
        seedList,
        function(i){
          inputData <- Visu_reduce( inputData, dims, "tsne", i, tsne_perplx = tsne_perplx )
          redCoords <- data.frame(inputData@reductions$tsne@cell.embeddings)
          pl <- DimPlot( inputData, reduction = "tsne", label = T, group.by = groupVar ) + NoLegend() 
          pl$labels$title <- paste0( "tSNE: random seed =", i  )
          list( plot= pl, redCoords=redCoords )
        })
    }
    names(redList) <- paste0( "Seed_", seedList )
    pgrobs <- lapply( redList, "[[", "plot" )
    redCoords <- lapply( redList, "[[", "redCoords" )
    pgrobs <- do.call( grid.arrange, c( pgrobs, ncol=3 ) )
    

  }
  return(list( Plots=pgrobs, redCoords=redCoords ) )
}













