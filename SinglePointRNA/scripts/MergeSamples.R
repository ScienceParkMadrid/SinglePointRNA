run_Merge <- function( inputTable=NULL, inputData=NULL, projectName="scRNA-seq",
                       phase=NULL, mode=NULL, norMode=NULL, nAnchors=2000
                       ){
  # Main function to merge/integrate scRNA-seq samples.
  ## ins:
  # inputTable: data frame containing paths and file type for scRNA-seq experiments
  # inputData: list of seurat objects 
  # projectName "Project_ABCD"
  # phase: the function is divided in two phases to allow for
  #     updates to be rendered on the UI.
  ## outs:
  # phase 1: list of Seurat objects
  # phase 2: Seurat objects result of the merging and normalization
  
  if(phase==1){ # loading datasets
    inputData <- lapply(
      1:nrow( inputTable ),
      function(i){
        Merge_loadData(
          inputTable[[ "fullPath" ]][i], inputTable[[ "inputType" ]][i],
          inputTable[[ "sampleName" ]][i])
      }    )
    names( inputData ) <- inputTable$sampleName
    inputData[ sapply( inputData, is.null ) ] <- NULL # in case any of the files are not found
    return( inputData ) # return list of seurat objects
  }
  
  if( phase==2 ){ # joining the datasets
    
    if( mode=="merge"){ # no correction!
      cellIds <- sapply( inputData, function(i){ i@project.name } )
      inputData <- merge( 
        x = inputData[[ 1 ]], y = inputData[ 2:length(inputData) ],
        project = projectName, add.cell.ids = cellIds )
      
      if( norMode == "SCT" ){
        inputData <- SCTransform(
          inputData, method = "glmGamPoi",
          vars.to.regress = NULL, 
          verbose = FALSE, return.only.var.genes = T )
        
        inputData <- RunPCA( inputData, features=rownames(inputData)  )
      } else if( norMode == "logNorm"){
        inputData <- NormalizeData( 
          inputData, normalization.method = "LogNormalize", scale.factor = 1e5  )
        inputData <- ScaleData( 
          inputData,
          vars.to.regress = NULL,
          features = rownames( inputData ) )
        inputData <- RunPCA( inputData, features=rownames(inputData)  )
      }
      
      return( inputData )
      
    } else if ( mode == "integrate" ){
      
     regVars <- Merge_checkRegressVars( inputData ) 
     
     if( any(!is.na(regVars)|is.null(regVars)) & norMode=="logNorm" ){
        
        inputData <- lapply(inputData, function(i) { # Find variable features
          i <- NormalizeData( i, normalization.method = "LogNormalize",  verbose = FALSE )
          i <- FindVariableFeatures(                 # return twice as many variable 
            i, selection.method = "vst",     # features for each dataset as required anchors
            nfeatures = 2* nAnchors )
        })

        features <- SelectIntegrationFeatures( inputData, nfeatures =  nAnchors )

        inputData <- lapply( inputData, function(i, features) {
          i <- ScaleData( i, features = features, verbose = FALSE )
          i <- RunPCA( i, features = features, verbose = FALSE )
        }, features=features )

        inputData <- FindIntegrationAnchors(
          inputData, anchor.features = features, reduction = "rpca" )
        
        inputData <- IntegrateData( inputData,   )
        inputData <- ScaleData( inputData, verbose = FALSE, vars.to.regress = regVars )
        
        inputData@project.name <- projectName
        return( inputData )
        
      } else if( any(!is.na(regVars)|is.null(regVars)) & norMode=="SCT" ){
        inputData <- lapply( 
          inputData,
          function(i, regVars){ Merge_SCTnorm( i, regVars ) }, regVars = regVars )
        gc()
        features <- SelectIntegrationFeatures( inputData, nfeatures = nAnchors)
        inputData <- PrepSCTIntegration( inputData, anchor.features = features)
        gc()
        inputData <- FindIntegrationAnchors(
          inputData, normalization.method = "SCT", anchor.features = features )
        inputData <- IntegrateData( anchorset = inputData, normalization.method = "SCT" )
        
        inputData@project.name <- projectName
        return( inputData )

      } else {
        return( NULL )
      }
     
    }
    
  }
  
}

##### auxiliary functions #####

Merge_updateInput <- function( inFile, inDir, inName, volumes, currentPaths=NULL,
                               sampleNameUpd=FALSE, SeuratList=NULL ){
  # Add new paths to scRNA-seq expertiments to the input dataframe This function
  # checks if 'inFile' or 'inDir' is already on the 'currentPaths' data frame, 
  # ands if isn't, it determines the input type and returns an updated data frame.
  # If RDS files are loaded, the sample names listed will be updated with the values
  ## ins:
  # inFile & inDir: Input values for the two shinyFiles UI buttons
  # volumes: array of reference paths.
  # currentPaths: dataframe of already added paths
  ## outs:
  # Updated data frame.

  if( ! sampleNameUpd ){ 
        # While the buttons aren't pushed, return the current DF
    if( all(is.integer( inFile )) & all(is.integer( inDir )) ){
      return( currentPaths )
    }
    
    # if any button is pushed, trigger a check and update:
    if( ! all( is.integer( inFile ) ) ){
      dataPath <- parseFilePaths(volumes, inFile)$datapath
    } else  if( !all( is.integer( inDir ) ) ) {
       dataPath <- parseDirPath(volumes, inDir )
    }
    return( Merge_addInput( dataPath, inName, currentPaths ) )
   
  } else { # in case RDS files were loaded, update the sample name with the value stored
    project_names <- sapply( SeuratList, function(i){ i@project.name })
    currentPaths$sampleName <- project_names
    return(currentPaths)
  }

}

Merge_addInput <- function( inPath, inName, currentPaths=NULL ){
  # This functions checks if 'inPath' is already on the 'currentPaths'
  # data frame, ands if isn't, it determines the input type and 
  # returns an updated data frame.
  ## ins:
  # inPath: path to a folder or file
  # currentPaths: dataframe of already added paths
  ## outs:
  # Updated data frame.
    
  if( length( inPath ) == 1 ){
    if( is.null( currentPaths ) ){
      Merge_checkInputFile( inPath, inName )
    } else {
      if( ! inPath %in% currentPaths$fullPath ){
        p <- Merge_checkInputFile( inPath, inName )
        rbind( currentPaths, p )
      } else {
        currentPaths
      }
    }
  } else if (length(inPath) > 1 ){
    do.call( rbind, lapply( inPath, Merge_addInput ))
  }
    
}


Merge_checkInputFile <- function( inFile, inName ){
  # Determine input file type (10X cell ranger output folder, h5 matrix,
  # plain text matrix...).
  # Returns data frame with 3 columns:
  #     Full path
  #     Dataset type [ 10X | H5 | plainText | RDS ]
  #     sample name
  
  if( is.null( inName )){
    sampleName <- gsub( "/$", "", gsub( ".*/", "", inFile ) )
    sampleName <- gsub( "_.*", "", sampleName )
  } else {sampleName <- inName  }
  
  if( dir.exists( inFile ) ){
    # 10X cell ranger results as input 
    inFile <- gsub( "/$", "", inFile )
    
    # if path is the top-level folder
    if("outs" %in% dir(inFile) ){
      inFile <- file.path( inFile, "outs/filtered_feature_bc_matrix/") 
      # if path is a mid-level folder
    }else if( "filtered_feature_bc_matrix" %in% dir(inFile) ){
      inFile <- file.path( inFile, "filtered_feature_bc_matrix/") 
      # if path is the bottom-level folder 
    } else if(any(grepl("barcodes", dir(inFile)) ) ){
      inFile <- file.path( inFile ) 
    }
    
    file_final <- dir(inFile)
    
    files10x <- paste0( inFile, "/" , file_final )
    
    if( ! all( file.exists( files10x ) ) ) {
      showNotification(
        paste0( "Input 10X files not found in '", inFile, "' .\n" ),
        type = "error", closeButton = TRUE, duration = NULL )
    } else {
      return( data.frame( fullPath = inFile, inputType = "10X", sampleName=sampleName ))
    }
    
  } else {
    # files as input ( h5, RDS, txt )
    if( file.exists( inFile ) ){

      sampleName <- gsub( "\\.[a-zA-Z0-9]+$", "",sampleName ) # remove extension
      
      if( grepl( "\\.RDS$", inFile, perl = T, ignore.case = T ) ){
        # RDS file as input
        return( data.frame( fullPath = file.path( inFile ), inputType = "RDS",sampleName=paste0( sampleName," (provisional)") ) )
        
      } else if( grepl( "\\.h5$", inFile, perl = T ) ){
        # H5 file as input 
        return( data.frame( fullPath = file.path( inFile ), inputType = "H5",sampleName=sampleName ))
      } else {
        # count matrix as input 
        return( data.frame( fullPath = file.path( inFile ), inputType = "plainText",sampleName=sampleName ))
      }
      
    } else {
      showNotification( paste0( "Count matrix file '", i, "' not found.\n" ),
        type = "error", duration = 10, closeButton = TRUE )
    }
  }
}

Merge_read10X <- function( inputPath, gene.column = 2 ){
  # Check if cellRanger output files are gzipped (ver>3.0) or not. Load with Seurat's
  # Read10X() or ReadMtx() respectively.
  ## Inputs:
  #   - inputPath: path to folder containing cellRanger output files
  #   - gene.column: column from features.tsv/genes.tsv to pull the gene IDs from
  #     (usually, 1 = ENSEMBL gene IDs and 2 = gene symbols )
  ## Outputs: dgCMatrix object to be used in CreateSeuratObject().
  
  files10X <- dir( inputPath )
  expected <- c( "barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz" )
  
  if( all( expected %in% files10X ) ){
    Read10X( data.dir = inputPath, gene.column = gene.column )
    
  } else {
    geneFile <- grep( "features.tsv", files10X, value = TRUE )
    cellFile <- grep( "barcodes.tsv", files10X, value = TRUE )
    matrixFile <- grep( "matrix.mtx", files10X, value = TRUE )
    
    if( length( geneFile )==0 ){ geneFile <- grep( "genes.tsv", files10X, value = TRUE )}
    
    if( any( sapply( list(geneFile, cellFile, matrixFile ), length) ==0 )  ){
      showNotification(
        HTML(paste0( "Input 10X files not found in the selected folder' .<br>",
                     "Expected 'barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz' or ",
                     "'barcodes.tsv', 'features.tsv', 'matrix.mtx'."        )),
        type = "error", closeButton = TRUE, duration = NULL )
      return(NULL)
    } else {
      ReadMtx( 
        mtx = paste0( inputPath, "/", matrixFile ), features = paste0( inputPath, "/", geneFile ), 
        feature.column = gene.column, cells =  paste0( inputPath, "/", cellFile )  )
    }
  }
}

Merge_loadData <- function( inputPath, inputType, projectName="scRNA-seq", gene.column=2 ){
  # Read single-cell RNA-seq count matrices in different formats (cellRanger output,
  # H5 or RDS )
  ## Inputs:
  # - inputPath: file/folder path
  # - inputType: [ 10X | H5 | RDS ], format of the scRNA-seq results.
  # - projectName: optional, name for the sample loaded
  # - gene.column: for cellRanger outputs, the gene ID format to select from 
  #   'features.tsv.gz' (defaults to gene symbols)
  ## Outputs: SeuratObject containing the experiment's count matrix and cell metadata.
  
  
  if( inputType == "10X" ){
    # Cell Ranger outputs
    inputData <- Merge_read10X( inputPath, gene.column )
    
  } else if( inputType == "H5" ){
    suppressWarnings( inputData <- Read10X_h5( inputPath ) )
    
  } else if( inputType == "plainText" ){
    inputData <- read.table( 
      file = inputPath, header = TRUE, as.is = TRUE, row.names = 1)
  }
  
  if( inputType == "RDS" ){
    readRDS( inputPath )
  } else {
    if( !is.null( inputData ) ){
      CreateSeuratObject( counts = inputData, project = projectName )
    }
  }
  
}


Merge_addSampleMD <- function( cellMeta, ExpMeta ){
  # Function to add sample-wide meta data to all the cells by their "orig.ident" value
  ## ins:
  # - cellMeta: current data frame of cell meta data from the seurat object
  # - ExpMeta: data frame with the sample names as row names.
  ## outs: data frame with a row per cell containing the matching sample-wide data.
  
  md_table <- Merge_readMDfile( ExpMeta )

  Samps <- unique( cellMeta$orig.ident )

  missed <- sum( !Samps %in% rownames(md_table) )

  if( missed == length(Samps) ){
    showNotification( "Row names in the table do not match sample names stored", 
      type = "error", duration = 10, closeButton = TRUE )
    return(NULL)
    
  } else if( missed %in% 1:(length(Samps)-1) ){
    missSamps <- paste0( Samps[ ! Samps %in% rownames(md_table)  ], collapse = ", ")
   
    showNotification( paste0(
      "Sample data missing for sample/s :", missSamps, ". 'NA's will be introduced"),
      type = "warning", duration = 10, closeButton = TRUE )
    matchTab <- md_table[ match( cellMeta$orig.ident, rownames( md_table  ) ), ]
  } else {

    matchTab <- md_table[ match( cellMeta$orig.ident, rownames( md_table  ) ), ]
  }
  rownames( matchTab ) <- rownames( cellMeta )
  return( matchTab )
}

Merge_readMDfile <- function( md_table ){
  # load a plain text or excel table containing sample metadata for a
  # single cell experiment.
  # Input:
  # - md_table: path to the file (xls, xlsx, txt, csv, tsv...)
  # Output: data frame.
  
  if( grepl("\\.xls", md_table ) | grepl("\\.xlsx", md_table ) ){
    inputMD <- as.data.frame( suppressMessages( read_excel( md_table )) )
    rownames( inputMD ) <- inputMD[, 1]
    inputMD <- inputMD[, -1, drop=F ]
  } else {
    if(grepl("\\.csv", md_table ) ){
      inputMD <- read.table( md_table, header = T, sep=",", row.names = 1 )
    } else if( grepl("\\.tsv", md_table ) ){
      inputMD <- read.table( md_table, header = T, sep="\t", row.names = 1 )
    } else {
      inputMD <- read.table( md_table, header = T, row.names = 1 )
    }
    
  }
  
  return( inputMD )
}

Merge_genomeSpecies <- function( x ){
  # Check the feature ids to determine the species
  # Inputs:
  # x: character vector - gene IDs
  
  if( all( grepl("^ENS", x ) ) ){ # ensembl IDs
    
    genePatts <- gen_loadGeneLists( "ENSids", "patterns" )
    Nmatch <- sapply( genePatts, function(i, x){ sum( grepl( i, x ) ) } ,x=x )
    IDtype <-"ENSEMBL Gene"

  } else {                        # gene symbols
    
    genePatts <- gen_loadGeneLists( "symbols", "patterns" )
    Nmatch <- sapply( genePatts,
      function(i, x){ sum( sum( grepl( i[1], x ) ) + sum( grepl( i[2], x ) ) )
      }, x=x
    )
    IDtype <-"Gene symbols"
  }
  
  if( any( Nmatch > 1 ) ){ Species <- names(genePatts)[ which( Nmatch == max(Nmatch) ) ]
  } else { Species <- "Other" }
  Species <- paste0( strsplit( Species,"")[[1]][1], ". ", 
    paste0(strsplit( Species,"")[[1]][-1], collapse = "" ), collapse = "" )
  
  return( c( Species = Species, IDtype = IDtype) )
}

Merge_checkRegressVars <- function( inputData ){
  # Function to check if any variables have been regressed in a SeuratObject.
  # Throws an error if a variable is not available in all objects of
  # a list.
  # Ins:
  #   - inputData: SeuratObject or list of SeuratObjects
  # Outs: array of variable names if any have been regressed, NULL otherwise

  if (is.list( inputData ) ){
    regVars <- unique( unlist( lapply( inputData, Merge_checkRegressVars ) ) )
    
    missingVars <- lapply(
      regVars,
      function(i){
        names(inputData)[sapply( inputData, function(j){ ! i %in% colnames(j@meta.data) } )]
        })
    names(missingVars) <- regVars
    
    if( any( sapply( missingVars, length )>0 ) ){
     
      missingVars <- missingVars[  sapply( missingVars, length )>0 ]
      stpMsg <- paste0( unlist( lapply( seq_along(missingVars), function(i){  paste0(
        names(missingVars)[i], " (",
        paste0( missingVars[[i]], collapse=",", recycle0 = FALSE),  ")" )
        } ) ), collapse=", " )
      
      showNotification(
        id="regVarError",
        paste0("Some regression variables are not present in all datasets: ", stpMsg ),
        type = "error", closeButton = TRUE, duration = 3600
        )
      return(NA)
      
    } else { return( regVars ) }
    
  } else {
    
    if( length( inputData@commands )>0 ){
      if( any( grepl("SCTransform", names(inputData@commands) ) |
            grepl("ScaleData", names(inputData@commands) ) )
      ){ # use whichever operation was done last
        com <- which( grepl("SCTransform", names(inputData@commands) ) |
                        grepl("ScaleData", names(inputData@commands) ) )
        if(length(com) > 1 ){com <- com[ length(com)]}
        
        return( inputData@commands[[ com ]]$vars.to.regress )
        
      } else {
        return( NULL )
      }
    } else { return( NULL ) } 
  }
  
 
    
}

Merge_checkAnchors <- function( inputData, nAnchors ){
  # Funtion to check the number of anchor features used for integration is not too
  # high compared to the number of genes available in the datasets
  ## ins:
  # InputData: list of Seurat objects
  # nAnchors: number of anchors set for integration
  ## outs:
  # If nAnchors is too high, return a vaue equal to half the number of common genes.
  
  # Find number of common genes
  avlblGenes <- table( unlist( lapply( inputData, rownames ) ) )
  avlblGenes <- length( avlblGenes[ avlblGenes == length( inputData )] )
  
  if( nAnchors >= 2 * avlblGenes ){ # too many anchors for the dataset
    nAnchors <- avlblGenes / 2
    return( nAnchors )
  } else { return( nAnchors ) }
  
}

Merge_SCTnorm <- function( inputData, regVars=NULL ){
  # Function to perform SCT normaliation on a Seurat object
  ## ins:
  # - inputData: seurat object
  # - regVars: names of variables to regress, if any.
  ## outs: Seurat object with new SCT slot in 'inputData@assays'ç
  
  if( length( inputData@commands )> 0 ){
    if( !any( grepl("SCTransform", names(inputData@commands) ) ) ){
      inputData <- SCTransform(
        inputData, assay = "RNA",  method = "glmGamPoi", vars.to.regress = regVars,
        verbose = FALSE, return.only.var.genes = T )
      
    } else {
      
      com <- which( grepl("SCTransform", names(inputData@commands) ) )
      if( length( com )>1){ com[ length(com ) ] }
      
      if( !all( regVars %in% inputData@commands[[com]]$vars.to.regress )){
        inputData <- SCTransform(
          inputData, assay = "RNA",  method = "glmGamPoi", vars.to.regress = regVars,
          verbose = FALSE, return.only.var.genes = T )
      }
    }
  } else {
    inputData <- SCTransform(
      inputData, assay = "RNA",  method = "glmGamPoi", vars.to.regress = regVars,
      verbose = FALSE, return.only.var.genes = T )
  }
  
  
    
  return( inputData )
}

Merge_getInputReport <- function( inputData, userOrg="Autodetect" ){
  # Generate a data frame containing a basic summary about the 
  # dataset loaded.
  ## Inputs:
  # - inputData: input dataset, Seurat Object
  # - userOrg: user-selected organism
  ## Outputs: data frame specifying species, gene ID type, etc.
  
  if( is.list( inputData ) ){    # Multiple samples
    inputDF <- do.call(rbind, lapply( inputData, Merge_getInputReport ) )
    inputDF$'Sample name' <- names(inputData)
    return( inputDF )
    
  } else {                       # One sample
    if( DefaultAssay( inputData ) != "RNA" ){
      DefaultAssay( inputData ) <- "RNA"
    }
    spec <- Merge_genomeSpecies( inputData@assays$RNA@data@Dimnames[[1]] )
    if( userOrg !="Autodetect"){ spec[1] <- userOrg }
    
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
    inputDF[['Sample name']] <- paste0(inputData@project.name, " (merged)" ) 
    return( inputDF )
  }
    
 
}


