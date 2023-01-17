run_featPl <- function( inputData, plotType="1D", featType='gene', features=NULL,
                        groupingVar=NULL, varTypesUser=NULL  ){
  # Options:
  # input dataset
  # plotType ['1D'|'2D' ]: generate violing plots or project the feature values onto
  #   UMAP/tSNE cell projections.
  # featType ['gene'|'md']: wheter the variables to plot are gene expression or 
  #   cell metadata.
  # features: character vector of feature names
  # groupingVar: extra variable to separate values by.
  # varTypesUser ['num'|'cat']: if a variable is stored as numerica, wheter to
  #   treat it as continuous or categorical (ie: 'cluster', 'replica', etc)
 
  
  if( featType=='gene'){
    g <- sort(inputData@assays$RNA@counts@Dimnames[[1]] ) 
    if( all( features %in% g ) ){ # plotting gene expression
      
      # Switch to lognormalized RNA counts
      inputData <- featPl_logNorm( inputData )
      
      cols <- 1 + as.numeric(length(features)>1)
      
      if( plotType=="1D" ){
        if( length( groupingVar ) == 2 ){
          pl <- VlnPlot( inputData, features = features, pt.size = 0, group.by = groupingVar[1],
                   split.by = groupingVar[2], split.plot = T , ncol = cols )
        } else {
          pl <- VlnPlot( inputData, features = features, pt.size = 0, group.by = groupingVar, 
                         ncol=cols )
        }
      } else if( plotType =="2D" ){
        
        if( length( groupingVar ) == 2 ){
          pl <- FeaturePlot( inputData, features = features, ncol =  cols, order=T )
        } else {
          pl <- FeaturePlot( inputData, features = features, ncol = cols, order=T )
        }
      }
    } else { return( list( Plot=NULL ) )}
  }
  
  if( featType == "md" ){
    
    varTypes <- featPl_varTypes( inputData, features, varTypesUser )

    pl <- lapply(
      seq_along( features ),
      function(i){
        
        rast <- ncol( inputData ) > 5000  # raster plots for larger datasets
        if(rast){ ptSize <-0.5 } else{ ptSize <-1 }
        
        if( varTypes[ i ] == "cat" ){
          
          if(  plotType=="2D" ){
            DimPlot( inputData, group.by = features[i], shuffle = TRUE,
                     label = TRUE, raster = rast )
          } else { # categorical values get barplots instead of violins
            
            if( !is.null( groupingVar ) ){
              
              tab <- do.call(rbind,lapply(
                unique( inputData@meta.data[, features[i]  ] ),
                function( j ){
                  
                  tab <- inputData@meta.data
                  tab[,features[i]  ] <- factor( tab[,features[i]  ] )
                  if( is.na(j)){
                    tab <- tab[ is.na( tab[,features[i]  ] ), ]
                  } else { 
                    tab <- tab[ tab[,features[i]  ] == j & !is.na(tab[,features[i] ]), ] 
                  }
                  data.frame( 
                    Group = sort( unique( tab[, groupingVar ] ) ),
                    VariableValue = paste0(  features[i], " - ", j ),
                    count=as.numeric( table( tab[, groupingVar ] ) ),
                    percentage  = 100* as.numeric( table( tab[, groupingVar ] ) ) / nrow(tab)
                  )
                }
              ))
              ggplot( tab ) + 
                geom_bar( aes(
                  x=VariableValue, y=percentage, fill=Group), stat="identity") + #, position="dodge" ) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, hjust=1)) +
                labs( x=features[i], y="Percentage of cells", fill="Group" )
              
            
            } else {
              tab <- data.frame( 
                    VariableValue = sort( unique(  inputData@meta.data[, features[i] ] ) ),
                    percentage  = 100* as.numeric(
                      table(  inputData@meta.data[, features[i] ] ) ) / nrow( inputData@meta.data )
                  )
              ggplot( tab ) + 
                geom_bar( aes(
                  x=VariableValue, y=percentage, fill=VariableValue), stat="identity", position="dodge",
                  show.legend = FALSE ) +
                theme_minimal() +
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
                labs( x=features[i], y="Percentage of cells" )
            }
          }
         
        } else if( varTypes[ i ] == "num" ){
          if( plotType == "2D"  ){
            pl <- FeaturePlot(
              inputData, features = features[i], order=T, raster = rast, pt.size = ptSize )
            if( grepl("^Uncertainty", features[i] ) ){ # change scale to 0-1 for uncertainty scores
              pl <- pl + scale_color_gradientn( colours = c("#d3d3d3","#0000ff"), limits=c(0,1) ) 
            }
          } else {
            if( length( groupingVar ) == 2 ){
              pl <- VlnPlot( inputData, features = features[i], pt.size = 0, group.by = groupingVar[1],
                             split.by = groupingVar[2], split.plot = T  )
            } else {
              pl <- VlnPlot( inputData, features = features[i], pt.size = 0, group.by = groupingVar)
            }
          }
        }
      }
    )
    if( length( features ) > 1){
      pl <-  do.call( grid.arrange, c(pl, ncol=2 ))
    }
  }
  
 # returns
  return( pl ) 
}


#### auxiliary functions ####

featPl_logNorm <- function( inputData ){
  
  if( DefaultAssay( inputData ) != "RNA" ){
    DefaultAssay( inputData ) <- "RNA"
  }
  NormalizeData( inputData, normalization.method = "LogNormalize", scale.factor = 1e5 )
}

featPl_widgetMD <- function( widgetID, metaData ){
  # Return a selectizeInput widget with the ID 'widgetID' to select available 
  # cell meta data variables.
  
  if( grepl("2$", widgetID ) ){ 
    maxIts <- 1 
    lab <- "Select a grouping variable:"
  } else {
    maxIts <- 2
    lab <- "Select up to two grouping variables:"
    }
  
  cellVars <- sort( colnames( metaData ) )
  return(
    renderUI({ selectizeInput(
      inputId = widgetID, label=lab, choices =  cellVars,
      selected = NULL, multiple = TRUE, options = list(maxItems = maxIts) )
    })
  )
}

featPl_widgetDims <- function( widgetID ){
  # return a radioButtons widget with the ID 'widgetID' to choose the plot type
  # (violin plots or projections onto a 2D embedding )
  return(
    renderUI({
      radioButtons( widgetID, label =  "Plot type:",
        choiceNames = c("2D projection", "Violin/Bar plot"), choiceValues = c("2D", "1D")  )
    })
  )
}

featPl_widgetVarType <- function( widgetID, inputData, vars ){
  # if any of the metadata of inputData variables in "vars" is numeric, 
  # return two widgets, an explaining text and a widget asking if the 
  # variable/s should be treated as continuous or categorical,
  # otherwise return NULL.
  
  md <- inputData@meta.data
  varType <- featPl_varTypes( inputData,  vars )
  
  if( length( vars ) > 1 ){

    if( any( varType  == "num" ) ){
      
      NumVar <- vars[ which( varType  == "num" ) ]
      
      if( length(NumVar) == 1 ){
        widgetLab <- HTML( paste0( 
          "The variable '", NumVar, "' is stored as numerical. ",
          "Should it be treated as categorical or continuous?") )
      } else {
        varNames <- paste0( NumVar, collapse = "', '" )
        
        widgetLab <- HTML( paste0(
          "The variables '", varNames, "' are stored as numerical. ",
          "Should they be treated as categorical or continuous?" 
        ))
      }
      
      buttons <- radioButtons(
        widgetID, "", choiceNames = c("Categorical", "Continuous"),
        choiceValues = c("cat","num") )
      
      return( list( text = widgetLab, buttons = buttons ) )
      
    } else {
      return( list( text = NULL, buttons = NULL ) )
    }
  } else {
    if( varType  == "num" ){
      widgetLab <- HTML( paste0( 
        "The variable '", vars, "' is stored as numerical. ",
        "Should it be treated as categorical or continuous?") )
      buttons <- radioButtons(
        widgetID, "", choiceNames = c("Categorical", "Continuous"),
        choiceValues = c("cat","num") )
      
      return( list( text = widgetLab, buttons = buttons ) )
    
    } else {
      return( list( text = NULL, buttons = NULL ) )
    }
  }
  
}

featPl_varTypes <- function( inputData, vars, widgetValue=NULL ){
  # return a character vector stating the treatment of each variable
  # according to the data type and the user's demand.
  
  NumVars <- sapply(
    vars,
    function(i, tab ){
      if( is.numeric( tab[,i] ) ){ "num" } else { "cat" }
    }, tab=inputData@meta.data
  )
  
  if( !is.null( widgetValue )){
    NumVars[ NumVars == "num" ] <- widgetValue
  }
  
  return(NumVars)
}




