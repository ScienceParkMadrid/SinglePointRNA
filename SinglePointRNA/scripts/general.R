gen_currentData <- function( inputData ){
  if( is.null( inputData )){
    renderUI({ HTML("No dataset loaded.") })
  } else {
    renderUI( {
      HTML(paste0( "<b>Dataset: </b>", inputData@project.name, " - ",
                   ncol(inputData), " cells." ) )}
    )
  }
}

gen_save <- function( inputData=NULL, tabs = NULL, type="meta", f_name=NULL ){
  
  if(type=="meta"){
    downloadHandler(
      filename = function(){paste0( inputData@project.name, "_Cell_Meta_Data.csv")},
      content = function(file){write.csv(inputData@meta.data, file, row.names = TRUE )}
    )
  } else if ( type=="seuratObject" ){
    downloadHandler(
      filename = function(){paste0( inputData@project.name, "_SeuratObject.RDS")},
      content = function(file){ saveRDS( inputData, file )}
    )
  } else if ( type=="tables" ){
    
    if( is.data.frame( tabs ) ){ # single tab
      if( suppressWarnings( all( sort(rownames( tabs[[i]]) == 1:nrow(tabs[[i]]) ) ) ) ) {
        #keep rownames only if they are informative
        rNames=FALSE } else { rNames=TRUE }
      
      downloadHandler(
        filename = function(){ paste0(f_name, ".csv") },
        content = function(file){write.csv( tabs, file, row.names = rNames, col.names = TRUE)}
      )
    } else { # list of tables
      
      xls_file <- createWorkbook()
      lapply( 
        seq_along( tabs ),
        function(i){
          if( !is.null( tabs[[i]] ) ){ # if an element is null, skip it
            
            if( suppressWarnings( all( sort(rownames( tabs[[i]]) == 1:nrow(tabs[[i]]) ) ) ) ) {
              #keep rownames only if they are informative
              rNames=FALSE } else { rNames=TRUE }
            
            sheetName <- names(tabs)[i]
            addWorksheet( xls_file, sheetName )
            writeData( xls_file, sheetName, tabs[[i]], rowNames = rNames, colNames = TRUE )
          }
        } )
      downloadHandler(
        filename = function(){paste0( f_name, ".xlsx" )},
        content = function(file){ saveWorkbook( xls_file, file )}
      )
    }
    
    
  }
  
}

gen_noDataNotif <- function( inputData, loading=FALSE ){
  # function to generate an error pop up notification if any of the run buttons
  # are clicked on without a dataset loaded.
  
  if( is.null( inputData ) & !loading ){
    showNotification("No dataset loaded!",
                     duration = 20, type="error", id="gen_noData" )    
  } 
  
  if( is.null( inputData ) & loading){
    showNotification("Select dataset files first",
                     duration = 20, type="error", id="gen_noData" )    
  } 
}

gen_loadGeneLists <- function( IDformat, section ){
  # load Cell Cycle / mitocondrial / ribosomal gene lists from 
  # "data/CellCycleGenes/" and "data/Ribo.MT.Genes
  # Ins:
  #  - IDformat: Gene ID format (gene symbols or ENSEMBL gene IDs)
  #  - section: what lists to load:
  #    - 'patterns': regular expressions.
  #    - 'CC': cell cycle genes
  #    - 'Rb.Mt': mitocondrial and ribosomal genes
  #    - 'pathways': Pathway databases
  
  if( section=="patterns" ){
    folder <- "data/GenePatterns/"
  } else if( section=="CC" ){
    folder <- "data/CellCycleGenes/"
  } else if( section=="Rb.Mt.NRY" ){
    folder <- "data/QC.Genes/"
  } else if( section=="pathways"){
    folder="data/Pathways/"
  }
  
  flist <- dir( folder )
  
  if( IDformat=="symbols" ){
    flist <- flist[grep("^symbols", flist)]
    namesL <-   gsub("symbols_", "", gsub( "\\.txt", "", flist  ) )
  } else if(IDformat=="ENSids" ){
    flist <- flist[grep("^ENSEMBL", flist)]
    namesL <- gsub("ENSEMBLids_", "", gsub( "\\.txt", "", flist  ) )
  }  
  
  geneList <- lapply( 
    flist,
    function(i){
      
      l <- readLines( paste0( folder ,i), skipNul = TRUE )
      l <- l[ !grepl( "^#", l ) & !grepl( "^$", l ) ]
      lName <- sapply( strsplit( l, ": "  ), "[[", 1 )
      l <- strsplit( sapply( strsplit( l, ": "  ), "[[", 2 ), ", ")
      l <- lapply( l, gsub, pattern="\"", replacement="" )
      l <- lapply( l, trimws )
      names(l) <- lName
      return( l )      
    }    )
  
  if( length( flist ) > 1 ){
    names(geneList) <- namesL
  } else {
    geneList <- as.list( geneList[[1]] )
  }
  
  return( geneList )
  
}

helpButton <- function (inputId, label=HTML("<b>?</b>"), icon = NULL, 
                        width = NULL, ...)  {
  value <- restoreInput(id = inputId, default = NULL)
  tags$button(
    id = inputId, 
    display="block", 
    style= paste0("heigth: 28px; width: 28px; color: #FFF; ",
                  "font-size: 18px;  font-family: inherit; line-height: 24px ",
                  "border-color: transparent; background: #3c8dbc; ",
                  "border-radius: 50%; border: 1px solid blue",
                  "margin: 0; padding: 1px 8px;",
                  "box-shadow: 0 2px 5px 0 rgba(0,0,0,0.18), 0 1px 5px 0 rgba(0,0,0,0.15)"),
    type = "helpButton", 
    class = "btn btn-default action-button", 
    `data-val` = value, 
    list(icon(icon), label), ...
  )
}

updButtonStates <- function( input, section, btnState ){
  # Function to return the state of help buttons to display help messages.
  ## ins:
  #  - input: shiny app input list
  #  - section: type of button. Currently only used with section="_help"
  #  - btnState: previos output from updButtonStates
  ## outs:
  #  - list of:
  #    - newBSts: Times each help button has been pressed along the current session
  #    - btn: if a help button has been clicked on, the ID of said button.
  
  inL <- reactiveValuesToList(input)
  inL <- inL[ grep( section, names(inL) )]
  
  currentButtons <- sapply(inL, "[[", 1 )
  currentButtons <- currentButtons[ order(names(currentButtons) ) ]
  
  if(is.null(btnState)){ # inicial states
    return( list( bSts= currentButtons, btn=NULL ) )
    
  } else {
    # if the change is triggered by a button being clicked on
    if( all( names( currentButtons ) %in% names(btnState) ) ){ 
      btn <- names(currentButtons)[ which( currentButtons != btnState ) ]
      if( length( btn ) == 0 ){ btn <- NULL }
      
      return( list( bSts = currentButtons, btn = btn ) )
      
    } else {      # if the change is triggered by a button appearing on the layout
      return( list(
        bSts = currentButtons, btn=NULL ) )
    }
  }
}

getMssg <- function( btn ){
  if(!is.null(btn)){
    renderUI(  HTML( Hmsg[[btn]] ) )  
  }
  
}






