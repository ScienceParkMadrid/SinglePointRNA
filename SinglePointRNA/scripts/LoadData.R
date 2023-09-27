run_Load <- function( input, metadata, projectName="scRNA-seq", cellGroups=NULL ){
  # Read single-cell RNA-seq count matrices in different formats (cellRanger output,
  # H5 or RDS )
  ## Inputs:
  #   - input: file/folder path
  #   - projectName: optional, name for the sample loaded
  #   - metadata: path to table file containing cell metadata (excel or plain text).
  #   - cellGroups: optional, variable to separate cells by.
  ## Outputs: SeuratObject containing the experiment's count matrix and cell metadata.
  
  

  #### input data and raw plots ####
  InputFile <- load_checkInputs( input, metadata )
  inputData <- load_dataset( 
    InputFile[[ "fullPath" ]], InputFile[[ "inputType" ]], projectName )
  
  if( !is.null( inputData ) ){
    if( !is.null( metadata ) ){
      inputData <- load_add_cell_md( inputData, metadata, cellGroups )
    } else {
      if( !is.null( cellGroups ) ){
        if( cellGroups %in% colnames( inputData@meta.data ) ){
          Idents( inputData ) <- cellGroups
        }
      }
    }
    return( inputData )
  }
}


#### auxiliary functions ####

load_ShinyUpldPath <- function( input, volumes ){
  # Function to determine which upload button has been clicked (choosing a file or
  # a folder) and generate the proper path to pass to the load functions.
  ## Inputs:
  # - input: shiny app's input list.
  # - volumes: char. vector, volumes or drives available used in UI file/folder selectors.
  
  file_chosen <- which(c(
    !is.integer( input$Load_inputFile ),
    !is.integer( input$Load_inputDir )
  ))
  
  if( file_chosen == 1 ){
    uplFile <- parseFilePaths(volumes, input$Load_inputFile)$datapath
  } else {
    uplFile <- parseDirPath(volumes, input$Load_inputDir)
  }
  if(!is.integer( input$Load_meta )){
    meta <- parseFilePaths(volumes, input$Load_meta)$datapath
  } else{ meta <- NULL }
  
  if(!is.integer( input$Load_DEG )){
    deg <- parseFilePaths(volumes, input$Load_DEG)$datapath
  } else{ deg <- NULL }
  
  return(list(
    dataset=uplFile,
    meta=meta,
    deg=deg
  ))
  
}

load_checkInputFile <- function( inFile ){
  # Determine input file type (10X cell ranger output folder, h5 matrix,
  # plain text matrix...).
  ## Inputs:
  # - inFile: file/folder path
  ## Output: Returns data frame with 2 columns:
  #     Full path
  #     Dataset type [ 10X | H5 | plainText | RDS ]
  
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
      return( data.frame( fullPath = inFile, inputType = "10X" ))
    }
    
  } else {
    
    if( file.exists( inFile ) ){
      
      if( grepl( "\\.RDS$", inFile, perl = T, ignore.case = T ) ){
        # RDS file as input
        return( data.frame( fullPath = file.path( inFile ), inputType = "RDS" ) )
        
      } else if( grepl( "\\.h5$", inFile, perl = T ) ){
        # H5 file as input 
        return( data.frame( fullPath = file.path( inFile ), inputType = "H5" ))
      } else {
        # count matrix as input 
        return( data.frame( fullPath = file.path( inFile ), inputType = "plainText" ))
      }
      
    } else {
      showNotification( paste0( "Count matrix file '", i, "' not found.\n" ),
        type = "error", duration = 10, closeButton = TRUE )
    }
  }
}

load_checkInputs <- function( input, metadata ){
  # Before loading data and running the analysis, check the selected files exist
  # and input options are as expected.
  ## Inputs:
  # - input: file/folder path
  # - metadata: path to table file containing cell metadata (excel or plain text).
  ## Outputs: if invalid input & path, shows notifications in the app, otherwise,
  # returns a data.frame with the full path and file type.
  
  if( ! is.null( input ) ){
    inputFile <- load_checkInputFile( input )
  } 
  # Check metadata table exists
  if( !is.null( metadata ) ){
    if(! file.exists( file.path( metadata ) ) ){
      showNotification( paste0("Metadata file '", metadata,  "' not found.\n"),
        type = "error", duration = 10, closeButton = TRUE )
    }
  }
  return(inputFile )
  
}

load_dataset <- function( inputPath, inputType, projectName="scRNA-seq", gene.column=2 ){
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
    inputData <- load_read10X( inputPath, gene.column )
    
  } else if( inputType == "H5" ){
    suppressWarnings( inputData <- Read10X_h5( inputPath ) )
    
  } else if( inputType == "plainText" ){
    inputData <- read.table( 
      file = inputPath, header = TRUE, as.is = TRUE, row.names = 1)
  }
  
  if( inputType == "RDS" ){
    inputData <- readRDS( inputPath )
  } else {
    if( !is.null( inputData ) ){
      inputData <- CreateSeuratObject( counts = inputData, project = projectName )
    }
  }
  
  return( inputData )
}

load_read10X <- function( inputPath, gene.column = 2 ){
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

load_md <- function( md_table ){
  # load a plain text or excel table containing cell or sample metadata for a
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

load_add_cell_md <- function( inputData, md_table, cellGroups=NULL ){
  # Read metadata and check rownames correspond with input cell IDs.
  # If a variable of interest is defined or the metadata table contains
  # only one variable, it will be set as the cell identity.
  ## Inputs: 
  # - inputData: input dataset, Seurat Object.
  # - md_table: path to the file (xls, xlsx, txt, csv, tsv...)
  # - cellGroups:  optional, variable to separate cells by.
  ## Output: Seurat Object including the cell metadata stored in md_table
  
  inputMD <- load_md( md_table ) 
  
  idMatch <- sum( colnames(inputData) %in% rownames( inputMD ) )
  
  if( idMatch == 0 ){
    showNotification( "Metadata and input cell IDs do not match. Discarding metadata.",
      type = "warning", duration=NULL, closeButton = TRUE, )
  } else {
    if( idMatch < ncol( inputData ) ){
      showNotification(
        paste0( 'Only ', idMatch, ' / ', ncol(inputData), ' cell IDs present in metadata.'),
        type = "warning", duration = 10, closeButton = TRUE)
    } 
    inputMD <- inputMD[ match( colnames(inputData), rownames( inputMD ) ), , drop=F ] 
    inputData <- AddMetaData( inputData, inputMD )
    
    if( !is.null( cellGroups ) ){
      Idents( inputData ) <- cellGroups
    }
    return( inputData )
  }
  
  
}

load_getInputReport <- function( inputData, userOrg="Autodetect" ){
  # Generate a data frame containing a basic summary about the 
  # dataset loaded.
  ## Inputs:
  # - inputData: input dataset, Seurat Object
  # - userOrg: user-selected organism
  ## Outputs: data frame specifying species, gene ID type, etc.
  
  if( is.list( inputData ) ){
    # Multiple samples
    inputDF <- do.call(rbind, lapply( inputData, load_getInputReport, userOrg=userOrg ) )
    inputDF$'Sample name' <- names(inputData)
    return( inputDF )
    
  } else {
    if( DefaultAssay( inputData ) != "RNA" ){
      DefaultAssay( inputData ) <- "RNA"
    }
    # One sample
    spec <- load_genomeSpecies( inputData@assays$RNA@data@Dimnames[[1]] )
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
    return( inputDF )
  }
  
}

load_genomeSpecies <- function( x ){
  # Check the feature ids to determine the species
  ## Inputs:
  # - x: character vector - gene IDs
  ## Outputs: character vector, species detected and gene ID type.
  
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

load_addDEG <- function( degPath  ){
  # Load differential gene expression results previously calculated and stored as
  # excel or RDS files.
  ## Inputs:
  # - degPath: path to the .xls / .xlsx file.
  ## Output: list containing:
  # 1st element: data frame, parameters used for differential expression
  # 2nd-last elements: data frames, differential expression results, with
  # one element per comparison.
  
  if( grepl("\\.xls", degPath ) | grepl("\\.xlsx", degPath ) ){
    sheets <- excel_sheets( degPath )
    params <- which( sheets == "Parameters" )
    deg <- which( sheets != "Parameters" )
   
    paramsTab <- as.data.frame( suppressMessages( read_excel( degPath, sheet = params )) )
    rownames( paramsTab ) <- paramsTab[,1]
    paramsTab <- paramsTab[, -1, drop=F ]

    degTab <- lapply(
      deg,
      function(i){
        tab <- as.data.frame( suppressMessages( read_excel( degPath, sheet = i ))) 
        rownames( tab ) <- tab[,1]
        tab[, -1, drop=F ]
      } )
    names(degTab) <- sheets[ deg ]
    return( list(
      Parameters=paramsTab, Tables = degTab
    ))
    
  } else {
    if( grepl("\\.RDS", degPath, ignore.case = TRUE ) ){
      readRDS( degPath  )
    }
  }
}




