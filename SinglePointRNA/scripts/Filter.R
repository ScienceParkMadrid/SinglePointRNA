run_Filter <- function( inputData, filterList, regressVars=NULL, mode="Scale" ){
  # Main function to filter cells and apply normlization and variance stabilization
  ## ins:
  # input dataset
  # basicFilters: counts, features, % mito RNA, % ribo RNA
  # extra filters: remove / keep categorical & numerical variables
  # regress metadata variable and mode
  ## outs:
  # Seurat object containing cells passing the filters and 
  # with added scaled/normalized expression data
  
  inputData <- Filter_filterCells( inputData, filterList, check=FALSE )
  
  inputData <- Filter_normRegress( inputData, regressVars, mode )

  return( inputData )
  
}

##### auxiliary functions #####


Filter_updateUIvals <- function( inputData, projectName ){
  # Function to update values displayed in the UI when a dataset is loaded
  
  currentInput <- renderText( {paste0(
    "Dataset: ", projectName, " - ", ncol(inputData), " cells." ) } )
  
  dataRange <- list(
    counts=range( inputData$nCount_RNA, na.rm = T, finite=F ), 
    feats=range( inputData$nFeature_RNA ,na.rm = T, finite=F ) 
  )
  if( "percent.mt" %in% colnames(inputData@meta.data) ){ 
    dataRange$mt=range(inputData$percent.mt, na.rm = T, finite=F ) 
    } else {dataRange$mt = c(0,100) }
  if( "percent.ribo" %in% colnames(inputData@meta.data) ){ 
    dataRange$ribo=range( inputData$percent.ribo, na.rm = T, finite=F )
    } else {dataRange$ribo = c(0,100) }
  
  # basic filters:
  rangeCount <- renderText( {paste0(
    "1: Counts per cell - Range: ", dataRange$counts[1], " - ",  dataRange$counts[2]) } )
  rangeFeat <- renderText( {paste0(
    "2: Features per cell - Range: ", dataRange$feats[1], " - ",  dataRange$feats[2]) } )
  rangeMT <- renderText( {paste0(
    "3: % of mitocondrial RNA - Range: ", dataRange$mt[1], " - ",  dataRange$mt[2]) } )
  rangeRibo <- renderText( {paste0(
    "4: % of ribosomal RNA - Range: ", dataRange$ribo[1], " - ",  dataRange$ribo[2]) } )
  
  minCount <- renderUI({ numericInput( "Filter_minCount", "Min:", min = dataRange$counts[1],
    max = dataRange$counts[2], value=dataRange$counts[1] ) }) 
  maxCount <- renderUI({ numericInput( "Filter_maxCount", "Max:", min = dataRange$counts[1],
    max = dataRange$counts[2], value=dataRange$counts[2] ) }) 
  minFeat <- renderUI({ numericInput( "Filter_minFeat", "Min:", min = dataRange$feats[1], 
    max = dataRange$feats[2], value=dataRange$feats[1] ) }) 
  maxFeat <- renderUI({ numericInput( "Filter_maxFeat", "Max:", min = dataRange$feats[1],
    max = dataRange$feats[2], value=dataRange$feats[2] ) }) 
  minMT <- renderUI({ numericInput( "Filter_minMT", "Min:", min = dataRange$mt[1],
    max = dataRange$mt[2], value=dataRange$mt[1] ) }) 
  maxMT <- renderUI({ numericInput( "Filter_maxMT", "Max:", min = dataRange$mt[1], 
    max = dataRange$mt[2], value=dataRange$mt[2] ) }) 
  minRibo <- renderUI({ numericInput( "Filter_minRibo", "Min:", min = dataRange$ribo[1],
    max = dataRange$ribo[2], value=dataRange$ribo[1] ) }) 
  maxRibo <- renderUI({ numericInput( "Filter_maxRibo", "Max:", min = dataRange$ribo[1],
    max = dataRange$ribo[2], value=dataRange$ribo[2] ) }) 
  
  # extra filters 
  Filter_numVars <- colnames( inputData@meta.data )[
    sapply( inputData@meta.data, typeof ) %in% c("numeric","double", "integer") ]
  Filter_numVars <- Filter_numVars[ !Filter_numVars  %in% c(
    "nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo" )]
  Filter_catVars <- colnames(inputData@meta.data)
  Filter_catVars <- Filter_catVars[ !Filter_catVars %in% c(
    "nCount_RNA", "nFeature_RNA", "percent.mt","percent.ribo" ) ]
  
  catVariableKeep <- renderUI({ selectInput(
    'Filter_catVariableKeep', '', c("None", Filter_catVars), selected = "None" ) })
  numVariableKeep <- renderUI({ selectInput(
    'Filter_numVariableKeep', '', c("None", Filter_numVars), selected = "None" )})
  catVariableRm <- renderUI({ selectInput(
    'Filter_catVariableRm', '', c("None", Filter_catVars), selected = "None" ) })
  numVariableRm <- renderUI({ selectInput(
    'Filter_numVariableRm', '', c("None", Filter_numVars), selected = "None" )})
  
  Filter_regressVars <- colnames(inputData@meta.data)
  if( all( c("S.Score", "G2M.Score") %in% Filter_regressVars )){
    Filter_regressVars <- c("None","cellCycle", Filter_regressVars)
    Filter_regressVars <- Filter_regressVars[ !Filter_regressVars  %in% c(
      "S.Score", "G2M.Score" ) ]
  }
  regressVariable <- renderUI({ selectInput(
    'Filter_regressVariable', '', c("None", Filter_regressVars),
    selected = "None", multiple = TRUE  )})
  
  return(list(
    currentInput = currentInput, rangeCount = rangeCount,
    rangeFeat = rangeFeat, rangeMT = rangeMT, rangeRibo = rangeRibo,
    
    minCount = minCount, maxCount = maxCount, minFeat = minFeat,
    maxFeat = maxFeat, minMT = minMT, maxMT = maxMT,
    minRibo = minRibo, maxRibo = maxRibo,
    
    catVariableKeep = catVariableKeep, catVariableRm = catVariableRm,
    numVariableKeep = numVariableKeep, numVariableRm = numVariableRm,
    
    regressVariable  = regressVariable
  ))
  
}


Filter_updateUIvals_extra <- function( 
  inputData, CatKeep = "None", CatRm = "None", NumKeep =  "None", NumRm = "None" ){
  # Function to update extra filter's available values displayed in the UI
  
  if( CatKeep != "None" ){
    catKeepVal <- renderUI( selectInput(
      "Filter_catValueKeep", '', multiple = T,
      choices=sort(unique(inputData@meta.data[[ CatKeep ]] ) ) ))
  } else { catKeepVal <- renderUI( selectInput('Filter_catValueKeep', '', "None" )) }
  
  if( CatRm != "None" ){
    catRmVal <- renderUI( selectInput(
      "Filter_catValueRm", '', multiple = T,
      choices=sort(unique(inputData@meta.data[[ CatRm ]] ) ) ))
  } else { catRmVal <- renderUI( selectInput('Filter_numVariableRm', '', "None" )) }
  
  if( NumKeep != "None" ){
    varRange <- range( inputData@meta.data[[ NumKeep ]], na.rm = T, finite=F )
    numKeepVal <- list(
      min = renderUI( numericInput( "Filter_numValue_minKeep",
        'Min:', min = varRange[1], max=varRange[2], value = varRange[1] )),
      max = renderUI( numericInput( "Filter_numValue_maxKeep",
        'Max:', min = varRange[1], max=varRange[2], value = varRange[2] ))
    )
  } else { NumKeepVal <- list(
    min=renderUI( numericInput('Filter_numValue_minKeep', '', min=NA, max=NA, value=NA, width = "150px" )),
    max=renderUI( numericInput('Filter_numValue_maxKeep', '', min=NA, max=NA, value=NA, width = "150px" ))) } 
  
  if( NumRm != "None" ){
    varRange <- range( inputData@meta.data[[ NumRm ]], na.rm = T, finite=F )
    numRmVal <- list(
      min = renderUI( numericInput( "Filter_numValue_minRm",
        'Min:', min = varRange[1], max=varRange[2], value = varRange[1] )),
      max = renderUI( numericInput( "Filter_numValue_maxRm",
        'Max:', min = varRange[1], max=varRange[2], value = varRange[2] ))
    )
  } else { numRmVal <- list(
    min=renderUI( numericInput('Filter_numValue_minRm', '', min=NA, max=NA, value=NA, width = "150px" )),
    max=renderUI( numericInput('Filter_numValue_maxRm', '', min=NA, max=NA, value=NA, width = "150px" ))) }
  
  return( list(
      catKeepVal=catKeepVal, catRmVal=catRmVal, NumKeepVal=NumKeepVal, numRmVal=numRmVal
  ) )
  
}



Filter_filterCells <- function( inputData, filterList, check=FALSE){
  # ins:
  # input dataset
  # basicFilters:
  #    - counts, features, % mito RNA, % ribo RNA
  # remove / keep categorical
  # remove / keep numerical

  if( ! "percent.mt" %in% colnames(inputData@meta.data ) ){
    filterList$Numeric$Keep <- filterList$Numeric$Keep[ - which(
      names(filterList$Numeric$Keep)=="percent.mt") ]
  }
  if( ! "percent.ribo" %in% colnames(inputData@meta.data ) ){
    filterList$Numeric$Keep <- filterList$Numeric$Keep[ - which(
      names(filterList$Numeric$Keep)=="percent.ribo") ]
  }
  
  cellSubset <- Filter_makeCellSubset( inputData, filterList )

  if(check==TRUE){
    return( list(
      Summary = Filter_makeFilteringSummary( filterList, cellSubset$failingCells ) ) )
  } else {
    inputData <- subset(
      inputData, cells = which( cellSubset[[ "subset" ]] ) )
    return( inputData )
  }
}

Filter_normRegress <- function( inputData, regressVars=NULL, mode="logNorm" ){
  # ins:
  # input dataset
  # regress metadata variables (optional)
  # normalization mode (logNormal | SCT)
  
  if( length(regressVars) == 1 ){ if(regressVars=="None"){ regressVars <- NULL }
  } else { regressVars <- regressVars[ regressVars != "None" ]  }
  if( "cellCycle" %in% regressVars ){
    regressVars <- c( regressVars[regressVars!="cellCycle"], "S.Score", "G2M.Score" )
    regressVars <- unique( regressVars[!is.na(regressVars) ] )
  }
  
  if( mode == "SCT" ){
    inputData <- SCTransform(
      inputData, method = "glmGamPoi",
      vars.to.regress = regressVars, 
      verbose = FALSE, return.only.var.genes = T )
    
    inputData <- RunPCA( inputData, features=rownames(inputData)  )
  } else if( mode == "logNorm"){
    inputData <- NormalizeData( 
      inputData, normalization.method = "LogNormalize", scale.factor = 1e5  )
    inputData <- ScaleData( 
      inputData,
      vars.to.regress = regressVars,
      features = rownames( inputData ) )
    inputData <- RunPCA( inputData, features=rownames(inputData)  )
  }
  
  gc()
  return( inputData )
}

Filter_makeCellSubset <- function( inputData, filterVars ){
  # Given a dataset and a list of variables and values, the funtion returns
  # the list of cells that meet the criteria in filterVars and a list of
  # the number of cells missing each criteria
  
  fails <- list( NA, NA, NA, NA )
  
  if( !is.null( filterVars$Numeric$Keep )){
    
    keepNum <- Filter_applyFilter( inputData, filterVars$Numeric$Keep, "num", F )
    fails[[1]] <- Filter_countFails( inputData, filterVars$Numeric$Keep, "num", 1 )

  } else { keepNum <- rep( T, ncol( inputData ) ) }
  if( !is.null( filterVars$Numeric$Remove )){
    removeNum <- Filter_applyFilter( inputData, filterVars$Numeric$Remove, "num", T )
    fails[[2]] <- Filter_countFails( inputData, filterVars$Numeric$Remove, "num", -1 )
  } else { removeNum <- rep( F, ncol( inputData ) ) }
  if( !is.null( filterVars$Categoric$Keep ) ){
    keepCat <- Filter_applyFilter( inputData, filterVars$Categoric$Keep, "cat", F )
    fails[[3]] <- Filter_countFails( inputData, filterVars$Categoric$Keep, "cat", 1 )
  } else { keepCat <- rep( T, ncol( inputData ) ) }
  if( !is.null( filterVars$Categoric$Remove ) ){
    removeCat <- Filter_applyFilter( inputData, filterVars$Categoric$Remove, "cat", T )
    fails[[4]] <- Filter_countFails( inputData, filterVars$Categoric$Remove, "cat", -1 )
  } else { removeCat <- rep( F, ncol( inputData ) ) }

  cellSubset <- ( keepNum & keepCat ) & !( removeNum | removeCat )
  failingCells <- unlist( fails, use.names = T )
  failingCells <- failingCells[ !is.na(failingCells) ]
  return( list( subset=cellSubset, failingCells=failingCells ) )
}

Filter_applyFilter <- function( inputData, varList, type="num", NA_trearment=FALSE ){
  # Given a dataset and a list of variables and values, the funtion returns
  # an array of T/F if each cell meets the conditions in varList
  
  if( type=="num"){
    CellPass <- sapply(
      seq_along( varList ),
      function(i){
        var.i <-  names( varList )[ i ]
        Filter_inRange( inputData[[ var.i ]][,1],  varList[[ i ]] )
      }
    )
  } else if( type=="cat"){
    CellPass <- sapply(
      seq_along( varList ),
      function(i){   
        var.i <-  names( varList )[ i ]
        inputData[[ var.i ]][,1] %in% varList[[ i ]]
      } )
  }
  
  if( length( varList ) > 1 ){
    CellPass <- apply(CellPass, 1, all)
  } else { CellPass <- CellPass[,1 ] }
  CellPass[ is.na(CellPass) ] <- NA_trearment
  
  return( CellPass )
}

Filter_countFails <- function( inputData, varList, type="num", value=1 ){
  # Count the number of cells that fails to meet each filtering criteria.
  # returns the list varList with an extra element in each filter, with the count of
  # cells failing to meet the standard.
  
  if( type=="num"){
    s <- sapply(
      seq_along( varList ),
      function(i, value){
        var.i <-  names( varList )[ i ]
        if( value == 1 ){
          sum( ! Filter_inRange( inputData[[ var.i ]][,1],  varList[[ i ]] ) ) 
        } else {
          sum( Filter_inRange( inputData[[ var.i ]][,1],  varList[[ i ]] ) ) 
        }
      }, value=value
    )
    
  } else if( type=="cat"){
    s <- sapply(
      seq_along( varList ),
      function(i, value ){
        var.i <-  names( varList )[ i ]
        if( value == 1 ){
          ncol(inputData) - sum( inputData[[ var.i ]][,1] %in% varList[[ i ]], na.rm = T )
        } else {
          sum( inputData[[ var.i ]][,1] %in% varList[[ i ]], na.rm = T )
        }
      }, value=value
    )
  }
  
  names( s ) <- names(varList)
  return( s )
}

Filter_makeFilteringSummary <- function( filterVars, fails ){
  
  if( !is.null( filterVars$Numeric ) ){
    
    ret <- list()
    NumTable <- data.frame(
      Variable = unlist(lapply( filterVars$Numeric, names ),use.names = F),
      Treatment = rep( names(filterVars$Numeric), sapply( filterVars$Numeric, length ) ),
      Range_Minimum = sapply( unlist(filterVars$Numeric, recursive = F ), "[[", 1 ),
      Range_Maximum = sapply( unlist(filterVars$Numeric, recursive = F ), "[[", 2 ),
      row.names = NULL
    )
    NumTable$Failing_Cells = fails[
      match( NumTable$Variable, names( fails ) )]
    ret <- c(ret, list( Numeric=NumTable ) )
  }
  
  if( !is.null( filterVars$Categoric ) ){
    CatTable <- data.frame(
      Variable = unlist(lapply( filterVars$Categoric, names ),use.names = F),
      Treatment = rep( c("Keep", "Remove"), sapply( filterVars$Categoric, length ) ),
      Value = sapply( unlist(filterVars$Categoric, recursive = F ), paste0, collapse=", " ),
      row.names = NULL
    )
    CatTable$Failing_Cells = fails[
      match( CatTable$Variable, names( fails ) )]
    ret <- c(ret, list( Categoric=CatTable ) )
  }
  
  return(ret)
}

Filter_inRange <- function( x, limits=NULL, minVal= -Inf, maxVal=Inf ){
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





