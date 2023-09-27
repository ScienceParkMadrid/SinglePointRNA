
# Loading packages and Sourcing scripts ----
library( shiny )
library( shinyjs )
library( shinydashboard )
library( shinyFiles )
library( dplyr )
library( ggplot2,  quietly = T ) 
library( grid,  quietly = T ) 
library( gridExtra,  quietly = T ) 
library( plotly )
library( ggalluvial )
library( openxlsx )
library( zoo )
library( readxl, quietly = T )
library( Seurat )
library( hdf5r )
library( HH )
library( stringr )

library( BiocParallel )
library( sctransform )
library( glmGamPoi )

library( scCATCH )

sapply(  dir("scripts",full.names = T), source )

# UI SIDE ----
ui <- dashboardPage( 

  title = "SinglePoint - scRNA-seq analysis GUI",
  
  dashboardHeader(
      tags$li(
        class = "dropdown",
        tags$style(".main-header {max-height: 87px}"),
        tags$style(".main-header .logo {height: 87px}")
      ),
      title = tags$img(src='nameLogo.png', width=280 ),
      titleWidth = 300,
      tags$li( class = "dropdown", 
        br(),
        p( align="left", "Contact with us through GitHub:",
        a( href="https://github.com/ScienceParkMadrid/SinglePointRNA",  
          img(src='GHicon.png', width="50px" ) ),
        HTML("&#160;&#160;&#160;")
      ) ),
      tags$li( class = "dropdown")
    ),
  

  dashboardSidebar(
    width = 300,
    br(), br(),
    sidebarMenu(
      menuItem("User's Guide", tabName = "intro", icon = icon("book") ),
      menuItem("Load data", tabName = "load", icon = icon("upload")),
      menuItem("QC analysis", tabName = "QC", icon = icon("broom")),
      menuItem("Filter Cells", tabName = "filter", icon = icon("filter")),
      menuItem("Sample Integration", tabName = "merge", icon = icon("blender")),
      menuItem("Explore Clustering Parameters", tabName = "exCluster", icon = icon("binoculars")),
      menuItem("Cluster", tabName = "cluster", icon = icon("object-group")),
      menuItem("2D visualizations", tabName = "visu", icon = icon("cookie")),
      menuItem("Differential Expression", tabName = "DiffExp", icon = icon("crosshairs")),
      menuItem("Dataset plots", tabName = "featPl", icon = icon("chart-bar")),
      menuItem("Cell Type Imputation", tabName = "cellType", icon = icon("bullseye")),
      menuItem("Altered Pathways", tabName = "paths", icon=icon("bezier-curve") )
    ),
    br(),img(src="logos.png", width=280), br()
    
  ),
  
  dashboardBody(
    useShinyjs(),
    
    tags$head(tags$link(rel="shortcut icon", href="browserIcon.ico")),
    
    tabItems(
      
      ## instructions ----
      
      tabItem(tabName = "intro",
        fluidRow(
          column(2),
          column(8,
          box(title="What's SinglePoint RNA©?", width = 12, collapsible = TRUE, collapsed = FALSE,
              includeHTML(path = "docu/Intro.html" ) ),
          box(title="User's guide", width = 12, collapsible = TRUE, collapsed = TRUE,
              includeHTML(path = "docu/guide.html" ) ),
          box(title="Case example 1: identify cell types in a sample",
              width = 12, collapsible = TRUE, collapsed = TRUE,
              includeHTML(path = "docu/Case1.html" ) ),
          box(title="Case example 2: Integrating biological replicates", 
              width = 12, collapsible = TRUE, collapsed = TRUE,
              includeHTML(path = "docu/Case2.html" ) ),
          box(title="Case example 3: Find the effects of a treatment on gene expression",
              width = 12, collapsible = TRUE, collapsed = TRUE,
              includeHTML(path = "docu/Case3.html" ) ), 
          HTML( paste0("<span style='font-size: 26px'>What's under the hood?  </span>",
            "<span style='font-size: 18px'>Learn more about computational science</span>")), 
          box(title="Dealing with large datasets: dimensional reduction",
            width = 12, collapsible = TRUE, collapsed = TRUE,
            includeHTML(path = "docu/CompBioGuides1.html" ) ),
          box(title="Classification",
            width = 12, collapsible = TRUE, collapsed = TRUE,
            includeHTML(path = "docu/CompBioGuides2.html" ) ),
          box(title="Optimization",
            width = 12, collapsible = TRUE, collapsed = TRUE, align="center",
            column(12, align="left", includeHTML(path = "docu/CompBioGuides3.html" )),
            column(12, align="center", uiOutput("intro_cb3pl" ),
              sliderInput(inputId = "intro_cb3plSel", "", 
                min = 1, max=10, value = 1, step = 1, width = "100%" )),
            column(12, align="left", includeHTML(path = "docu/CompBioGuides3b.html" ) ) )
         
          ),
          column(2))
         ),

      ## load data ----
      tabItem(tabName = "load",
        h2("Load Datasets"),
        fluidRow( column( width=6,
          fluidRow( 
            column(2,h4("Input files")),
            column(1, helpButton("Load_help1"))),
          p("Please select a file (H5 / plaint text / RDS) and click on 'Load dataset':"),
          
          shinyFilesButton(
            "Load_inputFile", "Select File", "Please select a file (H5 / plaint text / RDS)",
            multiple = FALSE, viewtype = "detail"),
          p("Or select a folder (10X output) and click on 'Load dataset':"),
          shinyDirButton(
            "Load_inputDir", "Select Folder", "Select a folder (10X output)", multiple = FALSE,
            viewtype = "detail"),
          fluidRow(
            column( 4, p("Project Name: "), 
                    textAreaInput( "Load_projectName",NULL, value = "scRNA-seq",  rows = 1)),
            column( 4, p("Organism:"), uiOutput("Load_org"))
            ),
          fluidRow(
            column(7, HTML("<h4>Additional cell metadata (optional):</h4>")),
            column(1, helpButton("Load_help2"))),
          shinyFilesButton(
            "Load_meta", "Select File", "Additional cell metadata", multiple = FALSE,
            viewtype = "detail"),
          p(),
          fluidRow(
            column(7,HTML("<h4>Load differential gene expression results:</h4>")),
            column(1, helpButton("Load_help3")) ),
          shinyFilesButton(
            "Load_DEG", "Select File", "Differential gene expression (Excel / RDS)", multiple = FALSE,
            viewtype = "detail"),
          HTML("<br><br>"),
          
          actionButton("Load_run", "Load dataset", icon = icon("upload"),
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95") ,
          actionButton("refresh_button", "Reset the app" , icon = icon("arrows-rotate"),
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          p(),p(),
          uiOutput( "Load_doneText" ),
          uiOutput( "Load_summ" ),
          uiOutput( "Load_summDE" )
        ),
        column(width=6, 
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
              uiOutput( "Load_chosenInput" ),
              uiOutput( "Load_chosenMD" ),uiOutput( "Load_chosenDEG" ), p(),
              uiOutput( "Load_help" ),p(),
              HTML(paste0("To test the application, download 10X Genomics' ",
                          "<a href=https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor",
                          "-v-3-chemistry-3-standard-3-0-0 style='color:rgb(23,59,99)'>1000 Peripheral Blood ",
                          "Mononuclear Cells</a> dataset in the three possible input formats. <br>")),
              downloadButton(
                "Load_dwlEx", "Download example datasets" )
          )
        )
        )
      ),
      
      ## QC analysis ----
      tabItem(tabName = "QC",
        h2("Quality Control"),
        fluidRow( column(width=9,
          h4("Options"),
          p("Grouping Variable (optional)"),
          fluidRow( column(3, uiOutput( "QC_groupingVar" )),column(1, helpButton("QC_help1")) ),
          fluidRow( column(3, uiOutput("QC_scoreCC")), column(1, helpButton("QC_help2") )),
          actionButton("QC_run", "Run QC analysis",
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          # outputs
          uiOutput( "QC_resText1" ), 
          tableOutput( "QC_table" ),
          uiOutput( "QC_resText2" ), 
          tableOutput( "QC_filters" ),
          uiOutput( "QC_resText3" ), 
          uiOutput( "QC_rawPlot"  ),
          uiOutput( "QC_resText4" ), 
          uiOutput( "QC_filteredPlot"),
          uiOutput( "QC_resText5" ),
          uiOutput( "QC_CCplot" )
        ),
        column(width=3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
              uiOutput( "QC_currentData" ), p(),
              uiOutput("QC_help")),
          box(title="Save & export", background = "green", width = 12, 
              downloadButton( "Save_QC_RDS", "Export dataset (RDS)" ),
              p(),
              downloadButton( "Save_QC_md", HTML("Export cell meta data<br>(plain text)") ),
              p(), HTML("To save plots, right-click on them.<br>")
          )
        )
        )
      ),
      
      ## filter/normalize cells ----
      
      tabItem(tabName = "filter",
        fluidRow( column(width=9,
          HTML("<h4>Basic Filters:</h4>"),
          fluidRow( 
            column( 3, textOutput( "Filter_rangeCount" )),
            column( 3, uiOutput( "Filter_minCount" )),
            column( 3, uiOutput( "Filter_maxCount" )) ),
          fluidRow( 
            column( 3, textOutput( "Filter_rangeFeat" )),
            column( 3, uiOutput( "Filter_minFeat" )),
            column( 3, uiOutput( "Filter_maxFeat" )) ),
          fluidRow( 
            column( 3,textOutput( "Filter_rangeMT" ) ),
            column( 3, uiOutput( "Filter_minMT" )),
            column( 3, uiOutput( "Filter_maxMT" )) ),
          fluidRow( 
            column( 3,textOutput( "Filter_rangeRibo" ) ),
            column( 3, uiOutput( "Filter_minRibo" )),
            column( 3, uiOutput( "Filter_maxRibo" )) ),
          fluidRow( 
            column(3,HTML("<h4>Additional filters:</h4><br>")),
            column(1, helpButton( "Filter_help1" ))),
          HTML("<br><b>Keep cells by a categorical variable:</b>"),
          fluidRow(
            column(4, uiOutput( "Filter_catVariableKeep" ) ),
            column(4, uiOutput( "Filter_catValueKeep" ))
          ),
          HTML("<br><b>Remove cells by a categorical variable:</b>"),
          fluidRow(
            column(4, uiOutput( "Filter_catVariableRm" ) ),
            column(4, uiOutput( "Filter_catValueRm" ))
          ),
          HTML("<br><b>Keep cells by a numerical variable:</b>"),
          fluidRow(
            column(3, uiOutput( "Filter_numVariableKeep" ) ),
            column(3, uiOutput( "Filter_numValue_minKeep" )),
            column(3, uiOutput( "Filter_numValue_maxKeep" ))
          ),
          HTML("<br><b>Remove cells by a numerical variable:</b>"),
          fluidRow(
            column(3, uiOutput( "Filter_numVariableRm" ) ),
            column(3, uiOutput( "Filter_numValue_minRm" )),
            column(3, uiOutput( "Filter_numValue_maxRm" ))
          ),
          fluidRow(
            column(4,HTML("<h4>Regress the effect of a variable on the dataset:</h4>")),
            column(1, helpButton("Filter_help2")) ),
          uiOutput( "Filter_regressVariable" ), 
          fluidRow(
            column(4,HTML("<h4>Normalization and variance stabilization:</h4>")),
            column(1, helpButton("Filter_help3")) ),
          radioButtons( "Filter_normMode", "", inline = T,
            choiceNames = c( "None", "Lognormalization", "SCTransform"),
            choiceValues = c( "None", "logNorm", "SCT")
          ), 
          fluidRow(
            column(3, actionButton("Filter_check", "Check filters",
              style="color: #fff; background-color: #3c8dbc; border-color: #015e95") ),
            column(3, actionButton("Filter_run", "Filter cells & normalize",
              style="color: #fff; background-color: #3c8dbc; border-color: #015e95") )
          ),
          #outputs
          uiOutput("Filter_doneText"),
          uiOutput("Filter_numVarsCheck"),
          uiOutput("Filter_catVarsCheck")
        ),
        column(width=3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
              uiOutput( "Filter_currentData" ),p(),
              uiOutput("Filter_help")),
          box(title="Save & export", background = "green", width = 12, 
              downloadButton( "Save_Filter_RDS", "Export dataset (RDS)" ),
              p(),
              downloadButton( "Save_Filter_md", HTML("Export cell meta data<br>(plain text)") )
          )     
        )
        )
              
      ),
      
      ## Sample Merge ----
      tabItem(tabName = "merge",
        h2("Sample integration"),
        fluidRow(column(9,
          fluidRow(
            column( 4, p("Project Name: "), 
                    textAreaInput( "Merge_projectName",NULL, value = "scRNA-seq",  rows = 1)),
            column( 4, p("Organism:"), uiOutput("Merge_org"))
          ),
          
          fluidRow(
            column(4, p( "Select input files or folders:"),
              column(6, shinyFilesButton(
                "Merge_addFile", "Add files","Please select a file (H5 / plaint text / RDS)",
                multiple = FALSE, viewtype = "detail") ),
              column(6, shinyDirButton(
                "Merge_addDir", "Add folders","Please select a folder (10X output)",
                multiple = FALSE, viewtype = "detail") ),
              uiOutput("Merge_FileMsg")
            ),
            column(4, p("Sample Name:" ), uiOutput("Merge_sName") ),
            column(2, br(),
              actionButton( "Merge_add", "Add sample" ,
                style="color: #fff; background-color: #3c8dbc; border-color: #015e95")
              )
          ),
          HTML("<br><h4>(Optional) Select a sample-wide metadata table: </h4>"),
          shinyFilesButton(
            "Merge_meta", "Select File", "Additional cell metadata", multiple = FALSE,
            viewtype = "detail"),
          h4("Dataset list:"),
          uiOutput("Merge_displayInput"),
          fluidRow( 
            column(2, radioButtons("Merge_mergeType", "Type of merge", choiceValues = c("1","2"),
              choiceNames = c("Simple merge", "Sample integration") ) ),
            column(1, helpButton("Merge_help1"))),
          fluidRow( 
            column(4,HTML("<br><b>Normalization and variance stabilization:</b>"),
              uiOutput( "Merge_mergeNormMode") ),
            column(1, p(HTML("<br>")), helpButton("Merge_help2")) ),
          uiOutput("Merge_mergeOpts_2"),
          actionButton( "Merge_run", "Merge datasets" ,
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          actionButton("refresh_button2", "Reset the app" , icon = icon("arrows-rotate"),
                       style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          uiOutput("Merge_ProgressText"),
          uiOutput( "Merge_inputTable" ) 
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
              uiOutput( "Merge_currentData" ), p(),
              uiOutput("Merge_help")),
          box(title="Save & export", background = "green", width = 12, 
              downloadButton( "Save_Merge_RDS", "Export dataset (RDS)" ),
              p(),
              downloadButton( "Save_Merge_md", HTML("Export cell meta data<br>(plain text)") )
          )    
        )
        )
      ),
      
      ## Explore clustering options ----
      tabItem(tabName = "exCluster",
        h2("Explore clustering options"),
        fluidRow( column(9, 
          fluidRow( 
            column(4, h3("Principal Component Analysis")),
            column(1, p(),helpButton("ExCl_help1"))
          ),
          checkboxInput( "exCluster_testNPC", "Test different numbers of PCs", value = FALSE ),
          conditionalPanel(
            condition = "input.exCluster_testNPC",
            fluidRow( 
              column( 4,numericInput( "exCluster_minSDdrop", "Minimum drop in standard deviation",
                value = 0.05, min = 0.01, max=0.1   )), 
              column(1, helpButton( "ExCl_help2" )  )),
            fluidRow( 
              column(4, numericInput( "exCluster_marginNPC", "Margin", value = 5, step = 1  )),
              column(1, helpButton( "ExCl_help3" ) ) ),
            uiOutput( "exCluster_jackStraw" ),
            actionButton("exCluster_run1", "Run test",
              style="color: #fff; background-color: #3c8dbc; border-color: #015e95")
          ),
          fluidRow( 
            column(4, h3("Clustering")),
            column(1, p(),helpButton("ExCl_help4"))
          ),
          checkboxInput( "exCluster_testRes", "Test different clustering resolutions", value = FALSE ),
          conditionalPanel( 
            condition = "input.exCluster_testRes",
            fluidRow( 
              column( 4, checkboxInput( 
                "exCluster_clUncert", "Check clustering uncertainty", value = FALSE )),
              conditionalPanel( 
                condition = "input.exCluster_clUncert",
                column( 4, numericInput( "exCluster_clUncertIters", "Uncertainty - iterations:", value=25 ))
              ),
              column( 1, helpButton( "ExCl_help5" )) 
            ),
            sliderInput(
              "exCluster_resolution", "Resolution range to test", 
              value = c(0,2), min = 0, max=2, step = 0.1, ticks = TRUE ),
            
            numericInput( "exCluster_finalnPC", "number of PCs (optional)", value = NA  ),
            fluidRow( column(2, checkboxInput( 
              "exCluster_plotEx", "Plot examples", value = FALSE )),
              column(1, helpButton( "ExCl_help6" ))
            ),
            actionButton("exCluster_run2", "Run test",
              style="color: #fff; background-color: #3c8dbc; border-color: #015e95")
            
          ),
          checkboxInput( "exCluster_plotly", "Generate interactive plots", value = TRUE ),
          
          uiOutput("exCluster_PCAresText1"),
          uiOutput("exCluster_nPCText"),
          fluidRow(
            column( 6, uiOutput("exCluster_PCAresText2"), uiOutput("exCluster_PCAcomps")),
            column( 6, uiOutput("exCluster_PCAresText3"), uiOutput("exCluster_PCCoCl"))
          ),
          uiOutput("exCluster_PCAresText4"),
          uiOutput("exCluster_PCflow"), 
          uiOutput("exCluster_ResresText1"),
          uiOutput("exCluster_ResresText2"),
          uiOutput("exCluster_ResCoCl"),
          uiOutput("exCluster_ResresText3"),
          uiOutput("exCluster_Resflow"),
          uiOutput("exCluster_ResresText4"),
          uiOutput("exCluster_ResUncertainty"),
          uiOutput("exCluster_ResresText5"),
          uiOutput("exCluster_ResExamples")
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
            uiOutput( "ExCl_currentData" ), p(),
            uiOutput( "ExCl_help" )),
          box(title="Save & export", background = "green", width = 12, 
            downloadButton( "Save_ExCl_nPCs", HTML("Export clustering at each<br>number of PCs (plain text)")),
            p(),
            downloadButton( "Save_ExCl_res", HTML("Export clustering at each<br>resolution (plain text)")),
            p(),
            downloadButton( "Save_ExCl_unCert", HTML("Export uncertainty at each<br>resolution (plain text)")),
            p(), HTML("To save plots, right-click on them.<br>")
         )
        )
        )
      ),
      
      ## Cluster cells ----
      tabItem(tabName = "cluster",
        h2("Cluster cells in the dataset"),
        fluidRow( column(9,    
          fluidRow(
            column(3, radioButtons( "cluster_autoParams", "Clustering parameters",
              choiceNames = c("Automatic selection", "User-defined parameters"), 
              choiceValues = c(1,2), selected = 1 )),
            column(1, helpButton("Cluster_help1"))),
          fluidRow(
            column( width=4, uiOutput( "cluster_resolution" )),
            column( width=4,uiOutput( "cluster_resValue" ))),
          fluidRow(
            column( width=4, uiOutput( "cluster_nPCs" )),
            column( width=4, uiOutput( "cluster_PCValue" ))),
          fluidRow(
            column( width=4, uiOutput( "cluster_makePlot" )),
            column( width=4, uiOutput( "cluster_embedding" ))),
          fluidRow(
            column( width=4, uiOutput( "cluster_clUncert" )),
            column( width=4, uiOutput( "cluster_clUncertIters" )),
            column(1, uiOutput("Cluster_help2")) ),
          actionButton( "Cluster_run", "Run clustering",
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          
          uiOutput("Cluster_resText1"),
          uiOutput("Cluster_resParams"),
          uiOutput("Cluster_resText2"),
          fluidRow(
            column( width = 2,uiOutput("cluster_summary1")),
            column( width = 2,uiOutput("cluster_summary2")),
            column( width = 2,uiOutput("cluster_summary3")),
          ),
          uiOutput("Cluster_resText3"),
          uiOutput("cluster_summPlots"),
          uiOutput("Cluster_resText4"),
          uiOutput("cluster_2Dplots") 
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
            uiOutput( "Cluster_currentData" ), p(),
            uiOutput( "Cluster_help" )),
          box(title="Save & export", background = "green", width = 12, 
            downloadButton( "Save_cluster_RDS", "Export dataset (RDS)" ),
            p(),
            downloadButton( "Save_cluster_md", HTML("Export cell meta data<br>(plain text)") ),
            p(), HTML("To save plots, right-click on them.<br>") )
        )
        )
      ),
      
      ## 2D visualizations ----
      tabItem(tabName = "visu",
        h2("Generate 2D visualizations"),
        fluidRow( column(9,
          fluidRow(
            column( 5, radioButtons( "Visu_saveMode", "", choiceNames = c(
              "Explore options", "Generate an embedding for the dataset"),
              choiceValues = c(1,2), selected = character(0) ) ),
            column( 1, helpButton("Visu_help1"))
          ),
           # mode 1: example plots
          fluidRow(
            column( 4, uiOutput( "Visu1_groupVar" ), uiOutput( "Visu1_paramMode" ) ),
            column( 2, uiOutput( "Visu1_algorithm" )),
            column( 1, uiOutput( "Visu_help2" ))
           ),
           fluidRow( # user-defined parameters
             column( 2, uiOutput( "Visu1_perplx" ) ),
             column( 2, uiOutput( "Visu1_nNeighbors" ) ),
             column( 2, uiOutput( "Visu1_minDist" ) ),
             column( 1, uiOutput( "Visu_help3" ) )
           ),
           fluidRow(
             column( 3, uiOutput( "Visu1_setSeed" )),
             column( 3, uiOutput( "Visu1_setSeedValue" ) ),
             column( 1, uiOutput( "Visu_help4" ) )
           ), 
           # mode 2: store reduction in dataset
           fluidRow(
             column( 4, uiOutput( "Visu2_paramMode" ) ),
             column( 2, uiOutput( "Visu2_algorithm" ) ),
             column( 1, uiOutput( "Visu_help5" ) )
           ),
           fluidRow( # user-defined parameters
             column( 2, uiOutput( "Visu2_perplx" ) ),
             column( 2, uiOutput( "Visu2_nNeighbors" ) ),
             column( 2, uiOutput( "Visu2_minDist" ) ),
             column( 1, uiOutput( "Visu_help6" ) )
           ),
           fluidRow(
             column( 3, uiOutput( "Visu2_setSeed" )),
             column( 3, uiOutput( "Visu2_setSeedValue" ) ),
             column( 1, uiOutput( "Visu_help7" ) )
           ),
           
           actionButton( "Visu_run", "Generate visualizations",
             style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
           
           uiOutput("Visu_resText1"),
           uiOutput("Visu_resPlots")
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
            uiOutput( "Visu_currentData" ),p(),
            uiOutput( "Visu_help")),
          box(title="Save & export", background = "green", width = 12, 
            downloadButton( "Save_Visu_RDS", "Export dataset (RDS)" ),
            p(),
            downloadButton( "Save_Visu_reds", HTML("Export cell coordinates<br>(plain text)") ),
            p(), HTML("To save plots, right-click on them.<br>") 
          )
        ) 
        )
      ),
      
      ## Differential expression ----
      tabItem(tabName = "DiffExp",
        h2("Differential expression"),
        fluidRow( column(9,
          fluidRow(
            column( 4, h4("Select the grouping variable"),
                    uiOutput("DiffExp_Grouping1" ) ),
            column( 4, uiOutput( "DiffExp_uiText1" ),
                    uiOutput( "DiffExp_Grouping2" )),
            column(1, helpButton("DiffExpr_help1"))
          ),
          h4("Choose an analysis mode:"),
          fluidRow( 
            column( 4,uiOutput( "DiffExp_mode" ) ),
            column( 1, helpButton("DiffExpr_help2") ) ),
            uiOutput( "DiffExp_uiText2" ),
          fluidRow(
            column( width = 4, 
                    uiOutput("DiffExp_1v1fg" ) ),
            column( width=4, uiOutput( "DiffExp_1v1bg" ))
          ),
          h4("DE Parameters"),
          fluidRow(
            column(width=4, HTML("<br>Discard genes expressed by a small percentage of the cells:")),
            column( width=4, numericInput(
              "DiffExp_minpct", "Minimum % of cells:", min=0, max=100, value=25 ) )
          ),
          fluidRow(
            column(width=4, HTML("<br>Set a Log<sub>2</sub>( Fold Change ) threshold:")),
            column( width=4, numericInput(
              "DiffExp_minLFC", "Minimum ±LogFC:", value=0.5 ) )
          ),
          fluidRow(
            column(width=4, HTML("<br>Define the mininum number of cells in a group:")),
            column( width=4, numericInput(
              "DiffExp_minCells", "Minimum # of cells:", value=20 ) )
          ),
          checkboxInput( "DiffExpr_rmMitoG", "Remove mitocondrial genes", value = TRUE ),
          checkboxInput( "DiffExpr_rmRiboG", "Remove ribosomal genes", value=TRUE ),
          actionButton( "DiffExp_run", "Calculate Differential Expression",
                style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          # outputs
          uiOutput( "DiffExpr_resText1" ),
          uiOutput( "DiffExpr_resParams" ),
          uiOutput( "DiffExpr_resText2" ),
          uiOutput( "DiffExpr_resSummary" ),
          uiOutput( "DiffExpr_resText3" ),
          uiOutput( "DiffExpr_resGroupSel" ),
          uiOutput( "DiffExpr_tabTopDE" ),
          uiOutput( "DiffExpr_resText4" ),
          uiOutput( "DiffExpr_plotHeatmap" )
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
              uiOutput( "DiffExpr_currentData" ),p(),
              uiOutput( "DiffExpr_help")),
          box(title="Save & export", background = "green", width = 12, 
              downloadButton( "Save_DEG_tabs", HTML("Export differential<br>expression (Excel)" )),
              p(), HTML("To save plots, right-click on them.<br>") )
        ) 
        )
      ),
      
      
      ## Dataset plots ----
      tabItem(tabName = "featPl",
        fluidRow( column(9,
          h2("Generate feature and meta-data plots"),
          HTML("<h3>Plot gene expression values</h3>"),
          fluidRow(
            column(5,
              selectizeInput(
                "featPl_feats", label="Select up to four genes:", choices = NULL,
                selected = NULL , multiple = TRUE, options = list(maxItems = 4) ) ),
              column(3, uiOutput("featPl_2D")),
              column(1, helpButton("featPl_help1"))
            ),
            fluidRow(
              column( 4, uiOutput("featPl_Grouping1" ) ),
              column(1, uiOutput("featPl_help2") )
            ),
            actionButton("featPl_run1", label = "Plot Gene Expression",
              style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
            uiOutput("featPl_plot1"),
            
            HTML("<h3>Plot cell metadata values</h3>"),
            fluidRow(
              column(5,
                selectizeInput(
                  "featPl_mdFeats", label="Select up to four variables:", choices = NULL,
                  selected = NULL , multiple = TRUE, options = list(maxItems = 4) ), 
                uiOutput( "featPl_varTypeText" ),
                uiOutput( "featPl_varType" ),
                fluidRow(
                  column(4, uiOutput( "featPl_Grouping2" )),
                  column(1, uiOutput("featPl_help4")) )
              ),
              column(3, uiOutput("featPl_2D_md")),
              column(1, helpButton("featPl_help3"))
            ),
            actionButton("featPl_run2", "Plot Cell Metadata",
              style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
            uiOutput("featPl_plot2"),
          ),
          column(3,
            box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE,
              uiOutput( "featPl_currentData" ), p(),
              uiOutput("featPl_help")),
            box(title="Save & export", background = "green", width = 12 ,
              downloadButton( "Save_featPl_RDS", "Export dataset (RDS)" ),
              p(),
              downloadButton( "Save_featPl_md", HTML("Export cell meta data<br>(plain text)") ),
            p(), HTML("To save plots, right-click on them.<br>") )
          )
        )  
      ),
      
      ## cell type imputation ----
      tabItem(tabName = "cellType",
        h2("Cell type Imputation"),
        fluidRow( column(9,
          fluidRow(
            column(4, uiOutput( "cType_mode" )), column(1, helpButton( "cType_help1" ))
          ),
          uiOutput( "cType_GroupVar" ),
          fluidRow( column(3, uiOutput( "cType_mtext" )), column(1, HTML("<br>"),uiOutput("cType_help2.3")) ),
          
           
          uiOutput( "cType_m1cancer" ),
          fluidRow(
             column( 3, uiOutput( "cType_m1tissues" ) ),
             column( 3, uiOutput( "cType_m1subtissues" ) )
          ),
           
          fluidRow(
             column( 4, uiOutput( "cType_m2ctName" ) ),
             column( 4, uiOutput( "cType_m2ctMarkers" ) ),
             column( 1, uiOutput("cType_m2add"),  uiOutput("cType_m2rm"),
                     style = "margin-top: 25px;")
          ),
          uiOutput( "cType_m2text2" ),
          uiOutput("cType_m2DisplayMark"),
          uiOutput( "cType_plot" ),
          actionButton( "cType_run", "Annotate possible cell types",
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
           
          uiOutput("cType_resText"),
           # scCATCH
          uiOutput( "cType_autoResTab" ),
          uiOutput( "cType_resText2"),
          uiOutput( "cType_dimplot" ),
           # manual
          uiOutput( "cType_manualTabC" ),
          uiOutput( "cType_resText3"),
          uiOutput( "cType_manualTabG" ),
          uiOutput( "cType_selectGroup"),
          uiOutput( "cType_ringPlot")
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
            uiOutput( "cType_currentData" ), p(),
            uiOutput("cType_help")),
          box(title="Save & export", background = "green", width = 12, 
            downloadButton( "Save_cType_RDS", "Export dataset (RDS)" ),
            p(),
            downloadButton( "Save_cType_md", HTML("Export cell meta data<br>(plain text)") ),
            p(),
            downloadButton( "Save_cType_cellType", HTML("Export cell type results<br>(plain text)") ),
            p(),
            downloadButton( "Save_cType_markerInfo", HTML("Export marker information<br>(plain text)") ),
            p(), HTML("To save plots, right-click on them.<br>")
          )
        )
        )
      ),
      
      ## altered pathways ----
      tabItem(tabName = "paths",
        h2("Altered pathways"),
        
        fluidRow( column(9, 
          fluidRow(
            column(9, HTML(
              paste0("<h3>Find metabolical pathways enriched in ",
                     "differentially expressed genes</h3>"))),
            column(1, helpButton("paths_help1"))
          ),
          uiOutput("paths_GroupVar"),
          HTML("<h4>Select pathway collections: </h4>"),
          uiOutput( "paths_CPselect" ), 
          HTML("<h4>Select differential expression results to test: </h4>"),
          uiOutput( "paths_DEGselect" ),
          
          actionButton("run_paths", "Run Analysis",
            style="color: #fff; background-color: #3c8dbc; border-color: #015e95"),
          
          # outputs
          uiOutput( "paths_resText1" ),
          uiOutput( "paths_summary" ),
          fluidRow( column(8, uiOutput( "paths_resText2" )), 
                    column(1, uiOutput("paths_help2") ) ),
          fluidRow(
            column(4, uiOutput("paths_resplot_text1"), uiOutput( "paths_resPlot_sel1" ) ),
            column(6, uiOutput("paths_resplot_text2"),uiOutput( "paths_resPlot_sel2" ))
          ),
          uiOutput( "run_pathsPlot"),
          uiOutput( "paths_resPlot" )
        ),
        column(3,
          box(title="Help & Info", background = "light-blue", width = 12, collapsible = TRUE, 
            uiOutput( "paths_currentData" ), p(),
            uiOutput("paths_help")),
          box(title="Save & export", background = "green", width = 12, 
            downloadButton( "Save_path", HTML("Export enrichment<br>results (excel)" )),
            p(),
            p(), HTML("To save plots, right-click on them.<br>")
          )
        )
        )
              
      )
      
    ) ) ) 

# SERVER SIDE ----
server <- function(input, output, session) {
  
  ## initialize storage objects ----
  
  gen_Msg <- reactiveValues( currentData=NULL  ) # store general messages
  inputDataset <- reactiveValues( Seurat=NULL, SeuratList=NULL )
  iD_summary <- reactiveValues( Table=NULL )
  QC_result <- reactiveValues( QC = NULL )
  Filter_check <- reactiveValues( Filters=NULL )
  Merge_inputDF <- reactiveValues( DF = data.frame( fullPath=c(), inputType=c(), sampleName=c())  )
  exF_Vars <- reactiveValues(CatKeep = "None", CatRm = "None", NumKeep =  "None", NumRm = "None" )
  Filter_fList <- reactiveValues( Filters = NULL )
  expCluster <- reactiveValues( nPCs=NULL, clResolution=NULL, clUncertainty = NULL )
  diffExpr <- reactiveValues( Parameters=NULL, Tables = NULL, Summary = NULL, Plots=NULL )
  Visu_data <- reactiveValues( Results=NULL )
  cType_markers <- reactiveValues( Ref=NULL, results=NULL )
  helpBtnStates <- reactiveValues( bSts=NULL, btn=NULL )
  DatasetPlots <- reactiveValues( geneExpr=NULL, datasetFeats=NULL )
  alteredPaths <- reactiveValues( Tables=NULL, Summary=NULL, backGrounds=NULL )
  
  observe({
    gen_Msg$currentData <- gen_currentData( inputDataset$Seurat )
  })
  
  observeEvent( input$refresh_button, { shinyjs::refresh() })
  observeEvent( input$refresh_button2, { shinyjs::refresh() })
  # Set user paths
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  
  ## Intro ----
  
  output$intro_cb3pl <- renderUI({
    plF <- paste0("docu/cb3_pl/round_", input$intro_cb3plSel, ".png")
    renderImage(list(src = plF, contentType = 'image/png', width="80%", heigth="80%" ), deleteFile = FALSE )
  })
  
  ## Load data ----
  
  shinyFileChoose(input, 'Load_inputFile', roots=volumes, filetypes=c('RDS','rds','Rds', 'h5', "txt"))
  shinyDirChoose(input, 'Load_inputDir', roots=volumes, allowDirCreate = FALSE)
  shinyFileChoose(input, 'Load_meta', roots=volumes, filetypes=c('xls',"xlsx","txt","tsv"))
  shinyFileChoose(input, 'Load_DEG', roots=volumes, filetypes=c('xls',"xlsx","RDS"))
  
  output$Load_summ <- renderUI({ HTML( "No datasets loaded" ) })
  
  output$Load_org <- renderUI({
    avblOrgs <- unique( gsub( "\\.txt", "", gsub( "^[A-Za-z]*_", "", dir("data/QC.Genes/") ) ) )
    avblOrgs <- sapply( strsplit( avblOrgs,""), function(i) {
      paste0( i[1], ". ", paste0(i[-1], collapse = "" ), collapse = "" )
    })
    selectInput("Load_org", NULL, multiple = F, choices = c("Autodetect", avblOrgs, "Other")) 
    })
  
  output$Load_chosenInput <- renderUI({
    if( !is.integer( input$Load_inputDir ) ){
      HTML(paste0( "<b>Input dataset folder selected:</b><br>&nbsp&nbsp&nbsp",
                   parseDirPath(volumes, input$Load_inputDir) ))
    } else if(!is.integer( input$Load_inputFile  )){
      HTML(paste0( "Input dataset file selected:<br>&nbsp&nbsp&nbsp",
                   parseFilePaths(volumes, input$Load_inputFile)$datapath ) )
    } else {""}
  })
  
  output$Load_chosenMD <- renderUI({
    if( !is.integer( input$Load_meta ) ){
      HTML(paste0( "Metadata file selected:<br>&nbsp&nbsp&nbsp",
                   parseFilePaths(volumes, input$Load_meta)$datapath ) )
    } else {""}
  })

  observeEvent( input$Load_DEG, { # differential expression is loaded independently
    if(  !is.integer( input$Load_DEG ) ){
      output$Load_chosenDEG <- renderUI({
        HTML(paste0( "Differential Expression Results selected:<br>&nbsp&nbsp&nbsp",
                     parseFilePaths(volumes, input$Load_DEG)$datapath ) )
      })
      deg <- parseFilePaths(volumes, input$Load_DEG)$datapath
    } else {
      output$Load_chosenDEG <- NULL
      deg <- NULL
    }
    if( !is.null( deg ) ){
      deg <- load_addDEG( deg )
      diffExpr$Parameters <- deg$Parameters
      diffExpr$Tables <- deg$Tables
    }
    output$Load_summDE <- renderUI({ renderTable({  diffExpr$Parameters  }, rownames = TRUE) })
  })
  
  observeEvent( input$Load_run, {
    
    showNotification( 
      HTML("Loading data<br>Please wait..."),
      duration = 20, type="message", id="Load_busy" )    
    anyFile <- any(c(
      !is.integer( input$Load_inputFile ),
      !is.integer( input$Load_inputDir ) ) )
    
    if( anyFile ){
      uplPaths <- load_ShinyUpldPath( input, volumes )
      
      inputDataset$Seurat <- run_Load(
        input=uplPaths$dataset, metadata=uplPaths$meta,
        projectName = input$Load_projectName  )
      
      if( !is.null( inputDataset$Seurat )){
        iD_summary$Table <- load_getInputReport( inputDataset$Seurat, input$Load_org )
        
        colnames(iD_summary$Table) <- gsub("\\.", " ", colnames(iD_summary$Table))

        output$Load_summ <- renderUI({ renderTable({ iD_summary$Table }) })
        output$Load_doneText <- renderUI({ HTML( "Done!" ) })
        
        if(!is.null( uplPaths$meta )) { load_missMD( uplPaths$meta, colnames( inputDataset$Seurat ) ) }
      }
    } else { gen_noDataNotif( inputDataset$Seurat, loading=T )}
  }, ignoreInit = TRUE)
  
  
  
  
  ## QC analysis ----
  
  #initial values:
  output$QC_groupingVar <- renderUI({selectInput('QC_groupingVar', '', "No dataset loaded" )})
  output$QC_currentData <- renderUI({ HTML("No dataset loaded.") })
  output$QC_scoreCC <- renderUI({ radioButtons(
    "QC_scoreCC", "Score Cell Cycle", choices = c("No"), selected = "No")  })
  
  observeEvent( inputDataset$Seurat, {
    output$QC_currentData <- gen_Msg$currentData
    output$QC_groupingVar <- renderUI({
      selectInput('QC_groupingVar', '',
        c("None", colnames( inputDataset$Seurat@meta.data ) )) })
    
    avblOrgs <- unique(gsub( "\\.txt", "",gsub( "^[A-Za-z]*_", "", dir("data/QC.Genes/"))))

    if( gsub("\\. ", "", iD_summary$Table$Species[1]) %in% avblOrgs ){
      output$QC_scoreCC <- renderUI({ radioButtons(
        "QC_scoreCC", "Score Cell Cycle", choices = c("Yes","No"), selected = "No")
      })
    } else {
      output$QC_scoreCC <- renderUI({ radioButtons(
        "QC_scoreCC", "Score Cell Cycle", choices = c("No"), selected = "No")
      })
    }

  })
  
  observeEvent( input$QC_run, {

    if( input$QC_groupingVar %in% c( "No dataset loaded", "None") ){
      cellGroups <- NULL
    } else { cellGroups <- input$QC_groupingVar }
    
    gen_noDataNotif( inputDataset$Seurat )
    if( ! is.null( inputDataset$Seurat ) ){
      showNotification( id = "QC_busy",
        HTML( paste0( "Testing cell quality<br>",
             "This might take a while...")),type = "message",duration = NULL) 
      QC_result$QC <- run_QCanalysis(
        inputDataset$Seurat, cellGroups, input$QC_scoreCC, iD_summary$Table )
    } 
    
  })
  # update input data and plot results
  observeEvent( QC_result$QC, {
    
    inputDataset$Seurat@meta.data <- QC_result$QC[["cellMD"]] 
    
    output$QC_resText1 <- renderUI({ HTML({"<h3>Summary: </h3>" }) })
    output$QC_table <- renderUI({ renderTable(
      { QC_result$QC[["reportTable"]] }, rownames = TRUE, colnames = TRUE) })
    
    output$QC_resText2 <- renderUI({ HTML({ "<h3>Default filters: </h3>" }) })
    output$QC_filters <- renderTable({
      QC_result$QC[["defaultFilters"]] },  rownames = TRUE, colnames = TRUE)
    
    output$QC_resText3 <- renderUI({ HTML({ "<h3>Raw data: </h3>" }) })
    output$QC_rawPlot <- renderUI({ renderPlot({
      gridExtra::grid.arrange( QC_result$QC[["rawPlot"]] ) }, 
      height = 800) 
    })
    
    output$QC_resText4 <- renderUI({ HTML({ "<h3>Filtered data: </h3>" }) })
    output$QC_filteredPlot <- renderUI({ renderPlot({
      gridExtra::grid.arrange( QC_result$QC[["filteredPlot"]] )  },
      height = 800 ) 
    })
    
    if( input$QC_scoreCC != "No"){
      output$QC_resText5 <- renderUI({ HTML({ "<h3>Cell Cycle scoring: </h3>" }) })
      output$QC_CCplot <- renderUI({ renderPlot({ 
        QC_result$QC[["CCscoringPlot"]]  }, height = 600, width = 800 )
      })
    }
    
    
  })
  
  ## Filter & Normalize ----
  
    ### f1: inicial values ----
  
  output$Filter_currentData <-  renderUI({ HTML("No dataset loaded.") })
  output$Filter_doneText <- NULL
  
  output$Filter_rangeCount <- renderText("1: Counts per cell - Range: No dataset loaded.")
  output$Filter_rangeFeat <- renderText("2: Features per cell - Range: No dataset loaded.")
  output$Filter_rangeMT <- renderText("3: % of mitocondrial RNA - Range: No dataset loaded.")
  output$Filter_rangeRibo <- renderText("4: % of ribosomal RNA - Range: No dataset loaded.")
  
  output$Filter_minCount <- renderUI({ numericInput( "Filter_minCount", "Min:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_maxCount <- renderUI({ numericInput( "Filter_minCount", "Max:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_minFeat <- renderUI({ numericInput( "Filter_minFeat", "Min:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_maxFeat <- renderUI({ numericInput( "Filter_minFeat", "Max:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_minMT <- renderUI({ numericInput( "Filter_minMT", "Min:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_maxMT <- renderUI({ numericInput( "Filter_minMT", "Max:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_minRibo <- renderUI({ numericInput( "Filter_minRibo", "Min:", min=NA, max=NA, value=NA, width = "150px" ) })
  output$Filter_maxRibo <- renderUI({ numericInput( "Filter_minRibo", "Max:", min=NA, max=NA, value=NA, width = "150px" ) })
  
  output$Filter_catVariableKeep <- renderUI( selectInput('Filter_catVariableKeep', '', "None" ))
  output$Filter_catValueKeep <- renderUI( selectInput('Filter_catValueKeep', '', "None" ))
  output$Filter_numVariableKeep <- renderUI( selectInput('Filter_numVariableKeep', '', "None" ))
  output$Filter_numValue_maxKeep <- renderUI( numericInput('Filter_numValue_maxKeep', '', min=NA, max=NA, value=NA, width = "150px" ))
  output$Filter_numValue_minKeep <- renderUI( numericInput('Filter_numValue_minKeep', '', min=NA, max=NA, value=NA, width = "150px" ))
  
  output$Filter_catVariableRm <- renderUI( selectInput('Filter_catVariableRm', '', "None" ))
  output$Filter_catValueRm <- renderUI( selectInput('Filter_catValueRm', '', "None" ))
  output$Filter_numVariableRm <- renderUI( selectInput('Filter_numVariableRm', '', "None" ))
  output$Filter_numValue_maxRm <- renderUI( numericInput('Filter_numValue_maxRm', '', min=NA, max=NA, value=NA, width = "150px" ))
  output$Filter_numValue_minRm <- renderUI( numericInput('Filter_numValue_minRm', '', min=NA, max=NA, value=NA, width = "150px" ))
  
  output$Filter_regressVariable <- renderUI( selectInput('Filter_regressVariable', '', "None" ))
  
  
    ### f2: updating on load ----
  observeEvent( inputDataset$Seurat, {
    
    FilterUpd <- Filter_updateUIvals( inputDataset$Seurat, input$projectName )
    
    output$Filter_currentData <-  gen_Msg$currentData
    output$Filter_rangeCount <-  FilterUpd[["rangeCount"]]
    output$Filter_rangeFeat <-  FilterUpd[["rangeFeat"]]
    output$Filter_rangeMT <-  FilterUpd[["rangeMT"]]
    output$Filter_rangeRibo <-  FilterUpd[["rangeRibo"]]
    
    output$Filter_minCount <-  FilterUpd[["minCount"]]
    output$Filter_maxCount <-  FilterUpd[["maxCount"]]
    output$Filter_minFeat <-  FilterUpd[["minFeat"]]
    output$Filter_maxFeat <-  FilterUpd[["maxFeat"]]
    output$Filter_minMT <-  FilterUpd[["minMT"]]
    output$Filter_maxMT <-  FilterUpd[["maxMT"]]
    output$Filter_minRibo <-  FilterUpd[["minRibo"]]
    output$Filter_maxRibo <-  FilterUpd[["maxRibo"]]
    
    output$Filter_catVariableKeep <-  FilterUpd[["catVariableKeep"]]
    output$Filter_catVariableRm <-  FilterUpd[["catVariableRm"]]
    output$Filter_numVariableKeep <-  FilterUpd[["numVariableKeep"]]
    output$Filter_numVariableRm <-  FilterUpd[["numVariableRm"]]
    
    output$Filter_regressVariable <-  FilterUpd[["regressVariable"]]
  })
  
  observeEvent( c(input$Filter_catVariableKeep, input$Filter_catVariableRm,
                  input$Filter_numVariableKeep, input$Filter_numVariableRm), {
    exF_Vars$CatKeep <- input$Filter_catVariableKeep
    exF_Vars$CatRm <- input$Filter_catVariableRm
    exF_Vars$NumKeep <- input$Filter_numVariableKeep
    exF_Vars$NumRm <- input$Filter_numVariableRm
    
    exF_Value <- Filter_updateUIvals_extra( inputDataset$Seurat,
      exF_Vars$CatKeep, exF_Vars$CatRm, exF_Vars$NumKeep, exF_Vars$NumRm )
    
    output$Filter_catValueKeep <- exF_Value$catKeepVal
    output$Filter_catValueRm <- exF_Value$catRmVal
    output$Filter_numValue_minKeep <- exF_Value$NumKeepVal$min
    output$Filter_numValue_maxKeep <- exF_Value$NumKeepVal$max
    output$Filter_numValue_minRm <- exF_Value$numRmVal$min
    output$Filter_numValue_maxRm <- exF_Value$numRmVal$max
  } )
    ### f3: run filter ----
  
  # prepare filter parameter list
  observeEvent(  c(input$Filter_run , input$Filter_check ),{
    
    gen_noDataNotif( inputDataset$Seurat )
    #output$Filter_doneText <- NULL
    
    if(!is.null( inputDataset$Seurat ) ){
      fl <- list(
        Numeric=list(
          Keep=list(
            nCount_RNA=c( input$Filter_minCount, input$Filter_maxCount ),
            nFeature_RNA=c( input$Filter_minFeat, input$Filter_maxFeat ),
            percent.mt=c( input$Filter_minMT, input$Filter_maxMT ),
            percent.ribo=c( input$Filter_minRibo, input$Filter_maxRibo ) ),
          Remove = NULL
        ),
        Categorical = list( Keep = NULL, Remove = NULL ))
      
      if( exF_Vars$NumKeep !="None" ){
        fl$Numeric$Keep <- c( fl$Numeric$Keep, list(
          c( input$Filter_numValue_minKeep, input$Filter_numValue_maxKeep  ) ) )
        names( fl$Numeric$Keep )[ 5 ] <- exF_Vars$NumKeep
      }
      if( exF_Vars$CatKeep !="None" ){
        fl$Categorical$Keep <- list( input$Filter_catValueKeep ) 
        names( fl$Categorical$Keep ) <- exF_Vars$CatKeep
      } 
      if( exF_Vars$NumRm !="None" ){
        fl$Numeric$Remove <- list( 
          c( input$Filter_numValue_minRm, input$Filter_numValue_maxRm  ) ) 
        names( fl$Numeric$Remove ) <- exF_Vars$NumRm
      } 
      if( exF_Vars$CatRm !="None" ){
        f <- c( input$Filter_catValueRm  )
        fl$Categorical$Remove <- list( input$Filter_catValueRm )
        names( fl$Categorical$Remove ) <- exF_Vars$CatRm
      }
      
      Filter_fList$Filters <- fl
    }
    
  }, ignoreInit = TRUE )
  
  # Run filters
  observeEvent( input$Filter_check,{
    if( !is.null( inputDataset$Seurat )){
      Filter_check$Filters <- Filter_filterCells( 
        inputData=inputDataset$Seurat, filterList = Filter_fList$Filters, check=TRUE)
      output$Filter_doneText <- renderUI({ HTML("<h3>Check completed.</h3>") }) 
    }
  } )
  
  observeEvent( input$Filter_run, {
    if( !is.null( inputDataset$Seurat )){
      showNotification( id = "filter_busy",
        HTML( paste0( "Filtering cells and normalizing expression values.<br>",
              "This might take a while...")),type = "message",duration = NULL) 
      
      Filter_check$Filters <- Filter_filterCells(
        inputData=inputDataset$Seurat,
        filterList = Filter_fList$Filters, check=TRUE)
    
      inputDataset$Seurat <- run_Filter(
        inputData=inputDataset$Seurat, filterList = Filter_fList$Filters, 
        regressVars = input$Filter_regressVariable,
        mode=input$Filter_normMode )
      output$Filter_doneText <- renderUI({ HTML("<h3>Filter/normalization completed.</h3>") }) 
    }
  })
  
  # filter outputs
  observeEvent( Filter_check$Filters,{
    output$Filter_numVarsCheck <- renderUI({ renderTable({
      Filter_check$Filters[["Summary"]][["Numeric"]]
    }, caption="Numeric filters" ) })
    output$Filter_catVarsCheck <- renderUI({ renderTable({
      Filter_check$Filters[["Summary"]][["Categoric"]]
    }, caption="Categorical filters") })
  } )
  
  
  ## Sample Merge ----
  
    ### m1: input display update ----
  
  shinyFileChoose(input, 'Merge_addFile', roots=volumes,filetypes=c('RDS','rds','Rds', 'h5', "txt"))
  shinyDirChoose(input, 'Merge_addDir', roots=volumes, allowDirCreate = FALSE)
  shinyFileChoose(input, 'Merge_meta', roots=volumes, filetypes=c('xls',"xlsx","txt","tsv"))
  
  output$Merge_currentData <- renderUI({ HTML("No dataset loaded.") })
  output$Merge_displayInput <- renderUI({HTML("No datasets loaded" ) })
  output$Merge_sName <- renderUI({
    textInput( "Merge_sName", NULL, width = "100%", placeholder = "Write a sample name") })  
  

  output$Merge_org <- renderUI({
    avblOrgs <- unique( gsub( "\\.txt", "", gsub( "^[A-Za-z]*_", "", dir("data/QC.Genes/") ) ) )
    avblOrgs <- sapply( strsplit( avblOrgs,""), function(i) {
      paste0( i[1], ". ", paste0(i[-1], collapse = "" ), collapse = "" )
    })
    selectInput("Merge_org", NULL, multiple = F, choices = c("Autodetect", avblOrgs, "Other")) 
  })
  
  observeEvent( c(input$Merge_addFile, input$Merge_addDir), {
    dp <- NULL
    if( ! all( is.integer( input$Merge_addFile ) ) ){
      dp <- parseFilePaths(volumes, input$Merge_addFile)$datapath
    } else  if( !all( is.integer( input$Merge_addDir ) ) ) {
      dp <- parseDirPath(volumes, input$Merge_addDir )
    }
    if( !is.null( dp ) ){
      output$Merge_FileMsg <- renderUI({ HTML( paste0( "selected:<br>", dp )) })
    } else {  output$Merge_FileMsg <- NULL }
   
  })
  
  observeEvent( c(input$Merge_add ), {
    Merge_inputDF$DF <- Merge_updateInput(
      input$Merge_addFile, input$Merge_addDir, inName = input$Merge_sName, volumes, Merge_inputDF$DF )
    output$Merge_FileMsg <- NULL
    output$Merge_sName <- renderUI({
      textInput( "Merge_sName", NULL, width = "100%", placeholder = "Write a sample name") })  
  })
  
  observeEvent( input$Merge_mergeType, {
    
    if( input$Merge_mergeType == "2" ){ # integration parameters
      
      output$Merge_mergeNormMode <- renderUI({ # remove "none" as an option for normalization
        radioButtons( "Merge_mergeNormMode", "", inline = T,
        choiceNames = c( "Lognormalization", "SCTransform"), choiceValues = c( "logNorm", "SCT") ) })
      
      output$Merge_mergeOpts_2 <- renderUI({ # change number of anchors to 3K if SCT
        if( input$Merge_mergeNormMode == "SCT"){
          numericInput( "Merge_nAnchors", "Number of integration features",  value = 3000 ) 
        } else {
          numericInput( "Merge_nAnchors", "Number of integration features",  value = 2000 ) }
      })
    } else {
      output$Merge_mergeNormMode <- renderUI({ radioButtons( "Merge_mergeNormMode", "", inline = T, 
        choiceNames = c( "None", "Lognormalization", "SCTransform"), 
        choiceValues = c( "None", "logNorm", "SCT") )
        })
    }
  })
  observeEvent( Merge_inputDF$DF, {
    output$Merge_displayInput <- renderUI({ renderTable(
      { Merge_inputDF$DF },colnames = T, rownames = F
    )})
  })
  
    ### m2: load and merge  ----
  
  # loading datasets
  observeEvent( input$Merge_run,{
    
    if( nrow( Merge_inputDF$DF )> 0 ){
      inputDataset$SeuratList <- run_Merge( Merge_inputDF$DF, phase=1)
      
      if( any( Merge_inputDF$DF$inputType == "RDS" ) ){ #update sample names
        Merge_inputDF$DF <- Merge_updateInput( currentPaths = Merge_inputDF$DF, 
          sampleNameUpd=TRUE, SeuratList=inputDataset$SeuratList )
      }
      
      # IU update:
      inReport <- Merge_getInputReport( inputDataset$SeuratList, input$Merge_org )
      showNotification( id = "merge_busy1",
        HTML( paste0( "Loading datasets.<br>",
              "This might take a while...")),type = "message",duration = NULL) 
      output$Merge_inputTable <- renderUI({ renderTable(
        inReport, rownames = TRUE, colnames = TRUE) })
      
    } else {
      gen_noDataNotif( inputDataset$SeuratList, loading=T  )
    }
  } )
  
  # Merging datasets
  observeEvent( inputDataset$SeuratList, {
    inReport <- Merge_getInputReport( inputDataset$SeuratList, input$Merge_org  )
    
    if( input$Merge_mergeType == "1" ){ # simple merge, no batch correction
      showNotification( id = "merge_busy2",
        HTML( paste0( "Merging datasets.<br>",
              "This might take a while...")),type = "message",duration = NULL) 
      inputDataset$Seurat <- run_Merge(
        inputData = inputDataset$SeuratList, phase=2, mode="merge",
        projectName = input$Merge_projectName )
    
    } else if( input$Merge_mergeType == "2"  ){ # integration 
      Merge_nAnchors <- Merge_checkAnchors( inputDataset$SeuratList, input$Merge_nAnchors )
      showNotification( id = "merge_busy2",
       HTML( paste0( "Integrating datasets.<br>",
             "This might take a while...")),type = "message",duration = NULL) 
      inputDataset$Seurat <- run_Merge(
        inputData = inputDataset$SeuratList, phase=2, mode="integrate",
        norMode=input$Merge_mergeNormMode, nAnchors=Merge_nAnchors,
        projectName = input$Merge_projectName )
    }
    
    if(!is.null(inputDataset$Seurat)){
      if(  !is.integer( input$Merge_meta ) ){ # add sample-wide meta
        mdPath <- parseFilePaths(volumes, input$Merge_meta)$datapath
        mdTab <-  Merge_addSampleMD( inputDataset$Seurat@meta.data, mdPath )
        if(!is.null( mdTab )){
          inputDataset$Seurat <- AddMetaData( inputDataset$Seurat, mdTab )    
        }
      }
      
      iD_summary$Table <- Merge_getInputReport( inputDataset$Seurat, input$Merge_org  )
  
      fullReport <- rbind( inReport, iD_summary$Table )
      colnames(iD_summary$Table) <- gsub("\\.", " ", colnames(iD_summary$Table))
      
      output$Merge_inputTable <- renderUI({ renderTable(
        fullReport, rownames = FALSE, colnames = TRUE) })
      
      inputDataset$SeuratList <- NULL # to save memory
      gc()
      output$Merge_ProgressText <- renderUI({renderText("Datasets merged and normalized.")})
    }
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent( inputDataset$Seurat,{
    if(!is.null(inputDataset$Seurat))
    output$Merge_currentData <- gen_Msg$currentData
  } )
  
  ## Explore clustering options ----
  
    ### exC1: Initial values & UI elements ----
  
  output$ExCl_currentData <- renderUI({ HTML("No dataset loaded.") })
  observeEvent( inputDataset$Seurat,{
    output$ExCl_currentData <- gen_Msg$currentData
  })
  observeEvent( input$exCluster_testNPC,{
    
    if( input$exCluster_testNPC & !is.null( inputDataset$Seurat ) ){
        if( DefaultAssay( inputDataset$Seurat ) != "SCT" ){
          output$exCluster_jackStraw <- renderUI( checkboxInput( 
            "exCluster_jackStraw", "Use JackStraw analysis on principal components", value = FALSE )
          )
        } else { output$exCluster_jackStraw  <- NULL }
    } else {
      output$exCluster_jackStraw  <- NULL
    }
  } )
  
    ### exC2: run cluster exploration  ----
  
  observeEvent(input$exCluster_run1,{
    
    gen_noDataNotif( inputDataset$Seurat )
    
    if(!is.null( inputDataset$Seurat )){
    
      if(ncol( inputDataset$Seurat )>10000 ){ defRes <- 0.2 
      } else if( ncol( inputDataset$Seurat )>5000 ){ defRes <- 0.3
      }else{ defRes <- 0.4 }
      
      showNotification( id = "exCl_busy1",
        HTML( paste0( "Testing clustering parameters.<br>",
              "This might take a while...")),type = "message",duration = NULL) 
      
      expCluster$nPCs <- run_clExplore( 
        inputDataset$Seurat, input$exCluster_marginNPC, 
        input$exCluster_minSDdrop, input$exCluster_jackStraw, 
        clResolution = defRes, section = "1", interactivePlots =  input$exCluster_plotly
      )
      
      gc()
    }
  })
  
  observeEvent(input$exCluster_run2,{
    gen_noDataNotif( inputDataset$Seurat )
    if(!is.null( inputDataset$Seurat )){
     
      res <- seq( from=input$exCluster_resolution[1], to=input$exCluster_resolution[2], by=0.1 )
      
      if( is.na(input$exCluster_finalnPC)){
        exCluster_finalnPC <- NULL } else {exCluster_finalnPC <- input$exCluster_finalnPC }
      showNotification( id = "exCl_busy2",
        HTML( paste0( "Testing clustering parameters.<br>",
              "This might take a while...")),type = "message",duration = NULL) 
      expCluster$clResolution <- run_clExplore(
        inputDataset$Seurat, final_nPC = exCluster_finalnPC,
        clResolution = res, plotExamples = input$exCluster_plotEx, 
        testUncertainty = input$exCluster_clUncert,
        uncIters = input$exCluster_clUncertIters,
        section = "2", interactivePlots =  input$exCluster_plotly
      )
      gc()
    }
  }  )
  
    ### exC3: outputs ----
  
  observeEvent( expCluster$nPCs, {
    output$exCluster_PCAresText1 <- renderUI({ HTML("<h3><b>Influence of the number of Principal Components</h></b>") })
    
    output$exCluster_nPCText <- renderUI({  HTML(
      paste0( "<h4>Recommended number of principal components: ",
              expCluster$nPCs$PCA_pcDrop$maxRecommendedPC, "<h4><br>" )
    ) })
    
    if( "JackStraw" %in% names(expCluster$nPCs$PCA_pcDrop$plots) ){
      output$exCluster_PCAresText2  <- renderUI( {renderText("JackStraw analysis of princiapl components:")})
      output$exCluster_PCAcomps <- renderUI( {renderPlot(
        expCluster$nPCs$PCA_pcDrop$plots$JackStraw
      ) })
    } else {
      output$exCluster_PCAresText2  <- renderUI( {renderText("Elbow plot: Decrease in Standard Deviation")})
      output$exCluster_PCAcomps <- renderUI( {renderPlot(
        expCluster$nPCs$PCA_pcDrop$plots$ElbowPlot
      ) })
    }
    
    output$exCluster_PCAresText3  <- renderUI( { HTML("Co-clustering conservation:")})
    output$exCluster_PCCoCl <- renderUI({ renderPlot(
      expCluster$nPCs$nPC$coClustering
    )})
    
    output$exCluster_PCAresText4  <- renderUI( {HTML("<br>Clustering flow plot:")})
    if( input$exCluster_plotly ){

      output$exCluster_PCflow <- renderUI({ plotlyOutput("plotly_1", height = "650px" ) })
      output$plotly_1 <- renderPlotly({ expCluster$nPCs$nPC$clusterFlowPlot })
    } else {
      output$exCluster_PCflow <- renderUI( {renderPlot(
        expCluster$nPCs$nPC$clusterFlowPlot
      ) } )
    }
  })
  
  observeEvent( expCluster$clResolution, {
    output$exCluster_ResresText1 <- renderUI({ HTML("<h3><b>Influence of SNN clustering resolution</h></b>") })
    
    output$exCluster_ResresText2  <- renderUI( {renderText("Co-clustering conservation:")})
    output$exCluster_ResCoCl <- renderUI( {renderPlot(
      expCluster$clResolution$Resolution$coClustering
    ) })
    
    output$exCluster_ResresText3  <- renderUI( { HTML("<br>Clustering flow plot")})
    if( input$exCluster_plotly ){
      output$exCluster_Resflow <- renderUI({ plotlyOutput("plotly_2", height = "650px" ) })
      output$plotly_2 <- renderPlotly({ expCluster$clResolution$Resolution$clusterFlowPlot })
      
    } else {
      output$exCluster_Resflow <- renderUI( {renderPlot(
        expCluster$clResolution$Resolution$clusterFlowPlot
      ) } )
    }
    
    if( !is.null( expCluster$clResolution$Uncertainty)  ){
      output$exCluster_ResresText4  <- renderUI( {renderText("Cluster assignment uncertainty:")})
      output$exCluster_ResUncertainty <- renderUI( {renderPlot(
        expCluster$clResolution$Uncertainty
      ) })
    }
    
    if( input$exCluster_plotEx ){
      output$exCluster_ResresText5  <- renderUI( {renderText("Example Plots")})
      output$exCluster_ResExamples <- renderUI( {renderPlot(
        grid.arrange( expCluster$clResolution$examplePlots )
      ) })
    }
    
  })
  ## Cluster cells ----
  
    ### cl1: Initial values and UI ----
  
  output$Cluster_currentData <- renderUI({ HTML("No dataset loaded.") })
  observeEvent( inputDataset$Seurat,{
    output$Cluster_currentData <- gen_Msg$currentData  })
  
  observeEvent( input$cluster_autoParams , {
    if( input$cluster_autoParams == 2 & !is.null(inputDataset$Seurat ) ){ # user-defined parameters
      output$cluster_resolution <- renderUI( radioButtons( "cluster_resolution",
        "Custering Resolution", choiceNames=c("Automatic",
        paste0( c("Low","Medium", "High"), " resolution" ), "Define value" ),
        choiceValues = c("Automatic", "low","mid","high","Define value")
      ))
      output$cluster_nPCs <- renderUI(radioButtons( "cluster_nPCs",
        "Number of principal components", choices=c("Automatic", "Define value" )
      ))
      output$cluster_makePlot <- renderUI( checkboxInput(
        "cluster_makePlot", "Generate 2D plot examples", value = FALSE ) )
      output$cluster_clUncert <- renderUI( checkboxInput(
        "cluster_clUncert", "Measure clustering uncertainty", value = FALSE ) )
      output$Cluster_help2 <- renderUI({ helpButton("Cluster_help2")})
      
    } else {
      output$cluster_resolution <- NULL
      output$cluster_resValue <- NULL
      output$cluster_nPCs <- NULL
      output$cluster_PCValue <- NULL
      output$cluster_makePlot <- NULL
      output$cluster_embedding <- NULL
      output$cluster_clUncert <- NULL
      output$cluster_clUncertIters <- NULL
      output$Cluster_help2 <- NULL
    }
  })
  observeEvent( input$cluster_resolution, {
    if( input$cluster_resolution =="Define value" ){
      output$cluster_resValue <- renderUI( numericInput("cluster_resValue", 
        label = "Resolution value:", min = 0, max=2, value=NA, step = 0.05  ))
    } else { output$cluster_resValue <- NULL }
  } )
  observeEvent( input$cluster_nPCs, {
    if( input$cluster_nPCs =="Define value" ){
      output$cluster_PCValue <- renderUI( numericInput("cluster_PCValue", 
        label = "Number of PCs:", value=NA  ))
    } else { output$cluster_PCValue <- NULL }
  } )
  observeEvent( input$cluster_makePlot, {
    if( input$cluster_makePlot==TRUE ){
      output$cluster_embedding <- renderUI( selectInput(
        "cluster_embedding", "Select embedding algorithm:", choices=c("tSNE", "UMAP", "Both") ))
    } else { output$cluster_embedding <- NULL }
  } )
  observeEvent( input$cluster_clUncert, {
    if( input$cluster_clUncert ){
      output$cluster_clUncertIters <- renderUI( numericInput(
        "cluster_clUncertIters", label = "Number of Iterations:", value=100  ))
    } else { output$cluster_clUncertIters <- NULL }
  } )
  
    ### cl2: run clustering and outputs ----
  
  observeEvent( input$Cluster_run, {
    # prepare input parameters
    
    gen_noDataNotif( inputDataset$Seurat )
    if( !is.null( inputDataset$Seurat) ){
      if( input$cluster_autoParams == 1 ){ # automatic parameters
        resolution <- "Automatic"
        nPCs <- "Automatic"
        makePlot <- FALSE
        testUncert <- FALSE
      } else {
        
        if( input$cluster_resolution == "Define value" ){
          resolution <- input$cluster_resValue
        } else { resolution <- input$cluster_resolution }
        if( input$cluster_nPCs == "Define value" ){
          nPCs <- input$cluster_PCValue
        } else { nPCs <- input$cluster_nPCs }
        if( is.null(input$cluster_makePlot) ){makePlot <- F
        } else if( input$cluster_makePlot==T ){ makePlot <- T } else { makePlot <- F } 
        if( is.null(input$cluster_clUncert) ){testUncert <- F
        } else if(input$cluster_clUncert==T){ testUncert <- T } else { testUncert <- F} 
      }
      
      parList <- Cluster_checkParams( inputDataset$Seurat, nPCs, resolution,
        makePlot, input$cluster_embedding, testUncert, input$cluster_clUncertIters )
      
      if( length(parList$repeats)>0 ){                       # check if clustering resolutions have already been used. 
        colMatch <- colnames(inputDataset$Seurat@meta.data)  # If they have, results will be overwritten
        colMatch <- colMatch[ ! colMatch %in% c( 
          paste0( "Clustering_res_", parList$repeats ),paste0( "Uncertainty_res_", parList$repeats )  )]
        inputDataset$Seurat@meta.data <-  inputDataset$Seurat@meta.data[ colMatch ]
      }
      showNotification( id = "Cluster_busy",
        HTML( paste0( "Clustering cells.<br>",
              "This might take a while...")),type = "message",duration = NULL) 
      inputDataset$Seurat <- run_Cluster( inputDataset$Seurat, parList )
      clSummary <- Cluster_summary( inputDataset$Seurat, parList )
      
      #### cl2.2: outputs ----
      output$Cluster_resText1 <- renderUI({ HTML( "<h3><b>Clustering results</h></b><br><br><h4>Parameters used:</h4>" ) })
      output$Cluster_resParams <- renderUI({ renderTable({
        clSummary$clusterParams
      })})
      output$Cluster_resText2 <- renderUI({HTML( "<br><h4>Distribution of cells among clusters:</h4>" )})
      
      if( length( clSummary$clusterDistribution ) == 1){ # 1 resolution value
        output$cluster_summary1 <- renderUI( renderTable(
          clSummary$clusterDistribution[[1]], rownames = F, colnames = T, caption= "",
          caption.placement = getOption("xtable.caption.placement", "top") ))
        output$cluster_summary2 <- NULL
        output$cluster_summary3 <- NULL
      } else { # 3 resolution values
        output$cluster_summary1 <- renderUI( renderTable(
          clSummary$clusterDistribution[[1]], rownames = F, colnames = T, caption= names(clSummary[[1]])[1],
          caption.placement = getOption("xtable.caption.placement", "top") ))
        output$cluster_summary2 <- renderUI( renderTable(
          clSummary$clusterDistribution[[2]], rownames = F, colnames = T, caption= names(clSummary[[1]])[2],
          caption.placement = getOption("xtable.caption.placement", "top") ))
        output$cluster_summary3 <- renderUI( renderTable(
          clSummary$clusterDistribution[[3]], rownames = F, colnames = T, caption= names(clSummary[[1]])[3],
          caption.placement = getOption("xtable.caption.placement", "top") ))
      }
      
      if( is.null( clSummary$summaryPlots$Uncertainty_Distribution ) ){
        output$Cluster_resText3 <- renderUI({ HTML( "<br>Number of cells in each cluster:" ) })
        output$cluster_summPlots <- renderUI( 
          renderPlot( { grid.arrange(clSummary$summaryPlots[[1]]) },
                      height =  400*length(parList$resolution), width = 600 ) )
      } else {
        output$Cluster_resText3 <- renderUI({ HTML(
          "<br>Number of cells and clustering uncertainty in each cluster:" ) })
        output$cluster_summPlots <- renderUI( renderPlot( {
          do.call(grid.arrange, c( clSummary$summaryPlots, ncol=2 ) )  },
          height = 500*length(parList$resolution), width = 1200 ) )
      }
      
      if( makePlot == T){
        
        clPlots <- Cluster_plot( inputDataset$Seurat, parList )
        if( length(clPlots) == 1  ){ # 1 resolution, 1 embedding
          emb <- c(tsne="tSNE", umap="UMAP")[ parList$make2Dplots$embedding ]
          output$Cluster_resText4 <- renderUI({ HTML( paste0( "<br><h4>",emb," reduction plot:</h4>" )) })
          output$cluster_2Dplots <- renderUI( renderPlot( {  clPlots[[1]][[1]] }, height = 600, width=800  ) )
        } else if( length(clPlots) == 3 ){  # 3 resolutions, 1 embedding
          emb <- c(tsne="tSNE", umap="UMAP")[ parList$make2Dplots$embedding ]
          output$Cluster_resText4 <- renderUI({ HTML( paste0( "<br><h4>",emb," reduction plot:</h4>" )) })
          output$cluster_2Dplots <- renderUI( renderPlot( {
            do.call( grid.arrange, c(clPlots, ncol=3 ) ) }, height = 450, width=1500  ) )
        } else { # 3 resolutions, 2 embeddings
          output$Cluster_resText4 <- renderUI({ HTML( "<br><h4>tSNE and UMAP reduction plots:</h4>" ) })
          if( length(parList$resolution)==1){ cols=2 } else {cols=3}
          output$cluster_2Dplots <- renderUI( renderPlot( {
            do.call( grid.arrange, c(clPlots, ncol=cols ) ) }, height = 1000, width=1500 ) )
        }
        
      }
    }
  } )
  
  ## Cell Visualizations ----
    ### cv1: initial values and UI updates ----
  output$Visu_currentData <- renderUI({ HTML("No dataset loaded.") })
  observeEvent( inputDataset$Seurat,{
    output$Visu_currentData <- gen_Msg$currentData  })
  
  observeEvent( input$Visu_saveMode, { 
    if(!is.null( inputDataset$Seurat )){
      cellVars <- sort( colnames( inputDataset$Seurat@meta.data ) )
    } else { cellVars <- c("No dataset loaded")}
    
    if( input$Visu_saveMode == 1){
      output$Visu_help2 <- renderUI({ helpButton("Visu_help2") })
      output$Visu_help4 <- renderUI({ helpButton("Visu_help4") })
      
      output$Visu1_groupVar <- renderUI({ selectInput( 'Visu1_groupVar', 
        "Grouping variable (color):", choices =  cellVars, selected = NULL , multiple = FALSE )
      })
      output$Visu1_algorithm <- renderUI({ radioButtons(
        "Visu1_algorithm", "Algorithm", choiceNames = c("tSNE", "UMAP" ), choiceValues = c("tsne","umap")
      ) })
      output$Visu1_paramMode <- renderUI({ radioButtons( "Visu1_paramMode", "General parameter definition",
        choiceNames = c("Default values", "Set based on dataset size", "Set manually" ),
        choiceValues = c( "default", "size", "userSet" )
      ) })
      output$Visu_help5 <- NULL
      output$Visu_help7 <- NULL
      output$Visu2_paramMode <- NULL
      output$Visu2_algorithm <- NULL
      output$Visu2_perplx <- NULL
      output$Visu2_nNeighbors <- NULL
      output$Visu2_minDist <- NULL
      output$Visu2_setSeed <- NULL
      output$Visu2_setSeedValue <- NULL
      
    } else if( input$Visu_saveMode==2 ){
      output$Visu_help5 <- renderUI({ helpButton("Visu_help5") })
      output$Visu_help7 <- renderUI({ helpButton("Visu_help7") })
      
      output$Visu2_algorithm <- renderUI({ radioButtons(
        "Visu2_algorithm", "Algorithm", choiceNames = c("tSNE", "UMAP" ), choiceValues = c("tsne","umap")
      ) })
      output$Visu2_paramMode <- renderUI({ radioButtons( "Visu2_paramMode", "General parameter definition",
        choiceNames = c("Default values", "Set manually" ),
        choiceValues = c( "default", "userSet" )
      ) })
      output$Visu_help2 <- NULL
      output$Visu_help4 <- NULL
      output$Visu1_groupVar <- NULL
      output$Visu1_paramMode <- NULL
      output$Visu1_algorithm <- NULL
      output$Visu1_perplx <- NULL
      output$Visu1_nNeighbors <- NULL
      output$Visu1_minDist <- NULL
      output$Visu1_setSeed <- NULL
      output$Visu1_setSeedValue <- NULL
    }
    
  })
  observeEvent( c(input$Visu1_paramMode, input$Visu1_algorithm),{ # Parameter setting options
    if( input$Visu1_paramMode  == "default" ){
      output$Visu1_setSeed <- renderUI({ radioButtons( "Visu1_setSeed", "Set seed:",
        choiceNames = c("Random sample", "Default value", "Set manually" ),
        choiceValues = c("random", "default", "userSet" ) )
      })
      output$Visu_help3 <- NULL
    } else if (input$Visu1_paramMode  == "size" ){
      output$Visu1_setSeed <- renderUI({ radioButtons( "Visu1_setSeed", "Set seed:",
        choiceNames = c( "Default value", "Set manually" ),
        choiceValues = c("default", "userSet" ) ) 
      })
      output$Visu_help3 <- NULL
    } else if( input$Visu1_paramMode == "userSet"){
      output$Visu_help3 <- renderUI({ helpButton("Visu_help3") })
      if( input$Visu1_algorithm =="tsne" ){
        output$Visu1_perplx <- renderUI({ numericInput(
          "Visu1_perplx", "Perplexity:", min=5, max=50, value=30 )})
        output$Visu1_nNeighbors <- NULL
        output$Visu1_minDist <- NULL
      } else if( input$Visu1_algorithm =="umap" ){
        output$Visu1_nNeighbors <- renderUI({ numericInput(
          "Visu1_perplx", "Number of neighbors:", min=5, max=50, value=30 )})
        output$Visu1_minDist <- renderUI({ numericInput(
          "Visu1_perplx", "Minimum distance:", min=0.001, max=0.5, value=0.3 )})
        output$Visu1_perplx <- NULL
      }
    }
  })
  observeEvent( c(input$Visu1_setSeed, input$Visu1_algorithm ), {
    if(!is.null(input$Visu1_setSeed) ){
      if( input$Visu1_setSeed=="userSet" ){
        if( input$Visu1_algorithm=="tsne"){val <- 1} else {val <- 42}
        output$Visu1_setSeedValue <- renderUI({ numericInput( "Visu1_setSeedValue",
          "Select seed value:", min=0, max=100, step = 1, value=val
        )})
      } else {
        output$Visu1_setSeedValue <- NULL
      }
    }
    
  } )
  observeEvent( c(input$Visu2_paramMode, input$Visu2_algorithm),{
    if( input$Visu2_paramMode  == "default" ){
      output$Visu_help3 <- NULL
      output$Visu2_setSeed <- renderUI({ radioButtons( "Visu2_setSeed", "Set seed:",
        choiceNames = c( "Default value", "Set manually" ),
        choiceValues = c( "default", "userSet" ) )
      })
      output$Visu2_nNeighbors <- NULL
      output$Visu2_minDist <- NULL
      output$Visu2_perplx <- NULL
    } else if( input$Visu2_paramMode == "userSet"){
      output$Visu_help6 <- renderUI({ helpButton("Visu_help6") })
      if( input$Visu2_algorithm =="tsne" ){
        output$Visu2_perplx <- renderUI({ numericInput(
          "Visu2_perplx", "Perplexity:", min=5, max=50, value=30 )})
        output$Visu2_nNeighbors <- NULL
        output$Visu2_minDist <- NULL
      } else if( input$Visu2_algorithm =="umap" ){
        output$Visu2_nNeighbors <- renderUI({ numericInput(
          "Visu2_nNeighbors", "Number of neighbors:", min=5, max=50, value=30 )})
        output$Visu2_minDist <- renderUI({ numericInput(
          "Visu2_minDist", "Minimum distance:", min=0.001, max=0.5, value=0.3 )})
        output$Visu2_perplx <- NULL
      }
    }
  })
  observeEvent( c( input$Visu2_setSeed, input$Visu2_algorithm ),{
    if( !is.null( input$Visu2_setSeed )){
      if( input$Visu2_setSeed=="userSet" ){
        if( input$Visu2_algorithm=="tsne"){val <- 1} else {val <- 42}
        output$Visu2_setSeedValue <- renderUI({ numericInput( "Visu2_setSeedValue",
          "Select seed value:", min=0, max=100, step = 1, value=val
        )})
      } else {
        output$Visu2_setSeedValue <- NULL
      }
    }
  } )
    ### cv2: execution and outputs ----
  
  observeEvent( input$Visu_run,{
    gen_noDataNotif( inputDataset$Seurat )
    if(!is.null(inputDataset$Seurat)){
      if( input$Visu_saveMode == 1 ){ # generate plots
        parList <- Visu_getParList( input )
        
        showNotification( id = "Visu_busy",
          HTML( paste0( "Generating plots.<br>",
                "This might take a while...")),type = "message",duration = NULL) 
        Visu_data$Results <- run_Visualization(
          inputDataset$Seurat, parList$groupVar, parList$algorithm,
          saveMode=1, parList$paraMode, parList$seedSet, parList$tsne_perplx,
          parList$umap_nNeighbors, parList$umap_minDist
        )
        
        output$Visu_resText1 <- renderUI({ HTML(
          "<h3>2D visualization examples:</h3>" )})
        output$Visu_resPlots <- renderUI({ renderPlot({
          grid.arrange( Visu_data$Results$Plots )
        }, height = 500*dim(Visu_data$Results$Plots)[1], width = 500*dim(Visu_data$Results$Plots)[2] ) })
        
      } else if( input$Visu_saveMode == 2 ){ # save one reduction in input dataset
        
        parList <- Visu_getParList( input )
  
        inputDataset$Seurat <- run_Visualization(
          inputDataset$Seurat, parList$groupVar, parList$algorithm,
          saveMode=2, parList$paraMode, parList$seedSet, parList$tsne_perplx,
          parList$umap_nNeighbors, parList$umap_minDist
        )
        alg <- c(tsne="tSNE", umap="UMAP")
        alg <- alg[ parList$algorithm ]
        output$Visu_resText1 <- renderUI({ HTML(
          paste0("<h3>",alg, " reduction saved in Seurat object.</h3>" ) )})
        
        if( any( grepl( "Clustering_", colnames(inputDataset$Seurat@meta.data ) ) ) ){
  
          clColumn <- max( which( grepl( # select the last clustering performed
            "Clustering_", colnames(inputDataset$Seurat@meta.data) ) ) )
          clColumn <- colnames(inputDataset$Seurat@meta.data)[ clColumn  ] 
          
          output$Visu_resPlots <- renderUI({ renderPlot({
            DimPlot( inputDataset$Seurat, reduction = input$algorithm, label=TRUE, group.by = clColumn )
          }, height = 600, width = 800 ) })
        } else { output$Visu_resPlots  <- NULL}
        
      }
    }
  })
  
  ## Differential expression ----
  
    ### d1: initial values and UI updates ----
  output$DiffExpr_currentData <- renderUI({ HTML("No dataset loaded.") })

  output$DiffExp_Grouping1 <- renderUI({ selectInput( 'DiffExp_Grouping1', '', choices =  c("None") ) })
  output$DiffExp_mode <- renderUI({
    radioButtons(
      "DiffExp_mode", label = '', selected = "1 VS rest",
      choices =  c("1 VS rest", "1 VS 1" ) ) 
  })
  
  observeEvent( inputDataset$Seurat, {
    output$DiffExpr_currentData <- gen_Msg$currentData
    
    cellVars <- sort( colnames( inputDataset$Seurat@meta.data ) )
    output$DiffExp_Grouping1 <- renderUI({ selectInput( 'DiffExp_Grouping1', 
      '', choices =  c("None", cellVars), selected = "None", multiple = FALSE )
    })
  })
  observeEvent( input$DiffExp_Grouping1, { # secondary grouping variable
    if(!is.null(input$DiffExp_Grouping1) & input$DiffExp_Grouping1 != "None" ){
      cellVars <- sort( colnames( inputDataset$Seurat@meta.data ) )
      cellVars <- cellVars[ !cellVars %in% input$DiffExp_Grouping1 ]
      output$DiffExp_Grouping2 <- renderUI({selectInput(
        'DiffExp_Grouping2', '', c("None", cellVars ), multiple = FALSE ) })
      output$DiffExp_uiText1 <- renderUI( HTML( "<br>(Optional) Select a secondary grouping variable:" ) )
    } else {
      output$DiffExp_Grouping2 <- NULL
    }
  })
  observeEvent( input$DiffExp_Grouping2,{
    if(!is.null(input$DiffExp_Grouping2) & input$DiffExp_Grouping2 != "None" &
       input$DiffExp_Grouping1 != "None" ){
      output$DiffExp_mode <- renderUI({
        radioButtons(
          "DiffExp_mode", label = '', selected = "1 VS rest",
          choices =  c("1 VS rest", "1 VS 1", "Conditional" ) ) 
      })
    } else {
      radioButtons(
        "DiffExp_mode", label = '', selected = "1 VS rest",
        choices =  c("1 VS rest", "1 VS 1" ) )
    }
    
  })
  
  observeEvent( c(input$DiffExp_mode, input$DiffExp_Grouping1 ),{
    
    if(!is.null(inputDataset$Seurat) & 
       all(!is.null( input$DiffExp_Grouping1 ) & input$DiffExp_Grouping1!= "None") ){
      if( input$DiffExp_mode =="1 VS 1" ){
      output$DiffExp_uiText2 <- renderUI( HTML( "<br>(Optional) Choose two groups to compare:" ) )
      
      if( is.null(input$DiffExp_Grouping2) | input$DiffExp_Grouping2 == "None" ){
        subgroup_choices <- unique( inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping1 ]] )
      } else {
        subgroup_choices <- unique( paste0( 
          inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping1 ]], " - ",
          inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping2 ]]        ) )
      }
      output$DiffExp_1v1fg <- renderUI({selectInput( 'DiffExp_1v1fg', '', 
        c( "None", subgroup_choices ), selected="None", multiple = FALSE ) })
      output$DiffExp_1v1bg <- renderUI({selectInput( 'DiffExp_1v1bg', '', 
        "None", selected="None", multiple = FALSE ) })
      } else if( input$DiffExp_mode =="Conditional" ){
        output$DiffExp_uiText2 <- renderUI( HTML( paste0(
          "<br>(Optional) Select a foreground value for '", input$DiffExp_Grouping1, "': ")))
        subgroup_choices <- unique( inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping1 ]] )
        output$DiffExp_1v1fg <- renderUI({selectInput( 'DiffExp_1v1fg', '', 
           c( "None", subgroup_choices ), selected="None", multiple = FALSE ) })
        output$DiffExp_1v1bg <- NULL 
      }else {
        output$DiffExp_uiText2 <- NULL
        output$DiffExp_1v1fg <- NULL
        output$DiffExp_1v1bg <- NULL 
      }
    }  
      
      
  }, ignoreInit = TRUE )
  
  observeEvent( input$DiffExp_1v1fg, {
    if( is.null(input$DiffExp_Grouping2) | input$DiffExp_Grouping2 == "None" ){
      subgroup_choices <- unique( inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping1 ]] )
    } else {
      subgroup_choices <- unique( paste0( 
        inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping1 ]], " - ",
        inputDataset$Seurat@meta.data[[ input$DiffExp_Grouping2 ]]        ) )
    }
    subgroup_choices <- subgroup_choices[ !subgroup_choices %in% input$DiffExpr_1v1fg ]
    
    if( input$DiffExp_mode != "Conditional"){
      output$DiffExp_1v1bg <- renderUI({selectInput( 'DiffExp_1v1bg', '', 
        c( "None", subgroup_choices ), selected="None", multiple = FALSE ) })
    }
  }, ignoreInit = TRUE )
  
    ### d2: run DiffExpr. ----
  

  
  observeEvent( input$DiffExp_run, {
    gen_noDataNotif( inputDataset$Seurat )
    
    if( !is.null( inputDataset$Seurat ) ){

      if( input$DiffExp_Grouping1 !="None" ){
        showNotification( id = "DiffExpr_busy",
          HTML( paste0( "Calculating Differential Expression<br>",
            "This might take a while...")),type = "message",duration = NULL) 
        
        de <- run_DiffExp( inputDataset$Seurat, input$DiffExp_Grouping1, input$DiffExp_Grouping2,
          input$DiffExp_mode, input$DiffExp_1v1fg, input$DiffExp_1v1bg, input$DiffExp_minpct,
          input$DiffExp_minLFC, input$DiffExp_minCells, input$DiffExpr_rmMitoG, input$DiffExpr_rmRiboG
        )
        removeNotification( id = "DiffExpr_busy")
        
        diffExpr$Tables <- de$Tables
        diffExpr$Parameters <- de$Parameters
        diffExpr$Summary <- de$Summary
        diffExpr$Plots <- de$Plots
      }
    }
  }  )
    ### d3: outputs ----
  observeEvent( c(diffExpr$Tables,inputDataset$Seurat), {
    
    if(!is.null( diffExpr$Tables ) & !is.null(inputDataset$Seurat) ){
      
      if( is.null( diffExpr$Summary ) & !is.null(inputDataset$Seurat)){ # DE loaded from file
        diffExpr$Summary <- DiffExpr_summary(
          inputDataset$Seurat,  diffExpr$Tables,
          Var1 = diffExpr$Parameters["Grouping Variable",],
          Var2 = diffExpr$Parameters["Secondary grouping variable",],
          mode = diffExpr$Parameters["Mode of comparison",],
          minLFC = diffExpr$Parameters["Minimum % of cells expressing a gene",],
          m1v1fg = diffExpr$Parameters["Foreground group",],
          m1v1bg = diffExpr$Parameters["Background group",]
        )
      }
      
      colnames(diffExpr$Summary ) <- gsub( "\\.", " ", colnames(diffExpr$Summary) )
      
      output$DiffExpr_resText1 <- renderUI(
        HTML("<h3>Differential expression results</h3><br><h4>Parameters used:</h4>") )
      output$DiffExpr_resParams <- renderUI( renderTable( diffExpr$Parameters, rownames = TRUE  ))
      output$DiffExpr_resText2 <- renderUI(
        HTML("<h4>Summary:</h4>") )
      output$DiffExpr_resSummary <- renderUI( renderTable( diffExpr$Summary  ))
      output$DiffExpr_resText3 <- renderUI( HTML(
        "<h4>Top differentially expressed genes</h4><br>Select a comparison:" ))
      output$DiffExpr_resGroupSel <- renderUI( selectInput( "DiffExpr_resGroupSel",
        "", names( diffExpr$Tables ), multiple = FALSE
      ) )
      if(!is.null( diffExpr$Plots$Heatmap )){
        output$DiffExpr_resText4 <- renderUI( HTML(
          "Heatmap of top differentially expressed genes"
        ))
         output$DiffExpr_plotHeatmap <- renderUI({ renderPlot({
        diffExpr$Plots$Heatmap }, height = 100+length( diffExpr$Tables)*100 )})
      }
     
    }
  } )
  
  observeEvent( input$DiffExpr_resGroupSel, {
    output$DiffExpr_tabTopDE <- renderUI( renderTable(
      diffExpr$Tables[[ input$DiffExpr_resGroupSel ]][1:20,], rownames = TRUE ) )
  } )
  
  ## Dataset plots ----
    ### pl1: initial values and UI updates ----
  output$featPl_currentData <- renderUI({ HTML("No dataset loaded.") })
  output$featPl_plot <- NULL
  
  observeEvent( iD_summary$Table, { # this part is triggered only on loading a new dataset.
    g <- sort(inputDataset$Seurat@assays$RNA@counts@Dimnames[[1]] )
    updateSelectizeInput(session=session, inputId = "featPl_feats", choices = g, server=T )
  } )
  
  observeEvent( inputDataset$Seurat,{ # Check if there's a tSNE/UMAP reduction stored and upd params
    output$featPl_currentData <- gen_Msg$currentData
    
    cellVars <- sort( colnames( inputDataset$Seurat@meta.data ) )
    dimReds <-  names(inputDataset$Seurat@reductions) 
    
    updateSelectizeInput( session=session, inputId = "featPl_mdFeats", choices = cellVars, server=T )
    
    
    if( any(c("tsne", "umap") %in% dimReds ) ){
      output$featPl_2D <- featPl_widgetDims("featPl_2D")
      output$featPl_2D_md <- featPl_widgetDims("featPl_2D_md")
    } else {
      output$featPl_Grouping1 <- featPl_widgetMD( "featPl_Grouping1", inputDataset$Seurat@meta.data )
      output$featPl_Grouping2 <- featPl_widgetMD( "featPl_Grouping2", inputDataset$Seurat@meta.data )

    }
  })
  
  observeEvent( c( input$featPl_2D, input$featPl_2D_md ),{ # Remove grouping vars for 2D projections
    if(!is.null(input$featPl_2D )){
    if( input$featPl_2D =="2D" ){ output$featPl_Grouping1 <- NULL 
    } else {
      output$featPl_Grouping1 <- featPl_widgetMD( 'featPl_Grouping1', inputDataset$Seurat@meta.data )
      output$featPl_help2 <- renderUI({ helpButton("featPl_help2") })
    }}
    if(!is.null(input$featPl_2D_md ))
    if( input$featPl_2D_md =="2D" ){ output$featPl_Grouping2 <- NULL 
    } else {
      output$featPl_Grouping2 <- featPl_widgetMD( 'featPl_Grouping2', inputDataset$Seurat@meta.data )
      output$featPl_help4 <- renderUI({ helpButton("featPl_help4") })
    }
  } )
  
  observeEvent( input$featPl_mdFeats,{
    if( !is.null( input$featPl_mdFeats ) ){
      
      wids <- featPl_widgetVarType( "featPl_varType", inputDataset$Seurat, input$featPl_mdFeats )
      output$featPl_varTypeText <- renderUI({ wids[[ "text" ]] })
      output$featPl_varType <- renderUI({ wids[[ "buttons" ]] })
    } else {
      output$featPl_varTypeText <- NULL
      output$featPl_varType <- NULL
    }
  } )
  
  
    ### pl2: generate plots ----
  observeEvent( input$featPl_run1,{ # Plot gene expression values
    
    gen_noDataNotif( inputDataset$Seurat )
    if(!is.null( inputDataset$Seurat )){
      if( !is.null(input$featPl_feats ) ){ 
        if( any(c("tsne", "umap") %in% names(inputDataset$Seurat@reductions) ) ){
          ptype <- input$featPl_2D
        } else { ptype="1D" }
        
        plH <- 500 * (1+ round(length( input$featPl_feats )/4 ))
        plW <- 600 * (1 + as.numeric( length( input$featPl_feats )>1) ) 
        
        DatasetPlots$geneExpr <- run_featPl( 
          inputDataset$Seurat, plotType=ptype, featType='gene', features=input$featPl_feats, 
          groupingVar=input$featPl_Grouping1  )
        gc()
        
        output$featPl_plot1 <- renderUI({ renderPlot({
          DatasetPlots$geneExpr }, height = plH, width = plW )
          })
      } 
    }
  })
  observeEvent( input$featPl_run2,{ # Plot cell metadata
    
    gen_noDataNotif( inputDataset$Seurat )
    if(!is.null( inputDataset$Seurat )){
      if( !is.null(input$featPl_mdFeats ) ){ 
        if( any(c("tsne", "umap") %in% names(inputDataset$Seurat@reductions) ) ){
          ptype <- input$featPl_2D_md
        } else { ptype="1D" }
        
        if(!is.null( input$featPl_varType )){
          varTypesUser <- input$featPl_varType
        }else{ varTypesUser <- "cat" }
        
        plH <- 500 * (1+ round(length( input$featPl_mdFeats )/4 ))
        plW <- 600 * (1 + as.numeric( length( input$featPl_mdFeats )>1) ) 
  
        DatasetPlots$datasetFeats <- run_featPl( 
          inputDataset$Seurat, plotType=ptype, featType='md', features=input$featPl_mdFeats, 
          groupingVar=input$featPl_Grouping2, varTypesUser =  varTypesUser )
        gc()
        
        output$featPl_plot2 <- renderUI({ renderPlot({
          if( is.grob( DatasetPlots$datasetFeats  ) ){
            grid.arrange( DatasetPlots$datasetFeats ) 
          } else { DatasetPlots$datasetFeats }
          }, height = plH, width = plW )
        })
      } 
    }
    
    
    
    
    
  })
  
  
  ## Cell type imputation ----
    ### ct1: initial values & UI updates ----
  
  output$cType_currentData <- renderUI({ HTML("No dataset loaded.") })
  observeEvent( inputDataset$Seurat,{
    output$cType_currentData <- gen_Msg$currentData  })
  
  # Add the option to use scCATCH when the dataset is human/mouse & IDs are symbols
  output$cType_mode <- renderUI({ radioButtons(
    "cType_mode", "Annotation mode", choiceNames = c("Manual" ),choiceValues = c(2))})
  
  observeEvent( iD_summary$Table, {
    if( iD_summary$Table$Species[1] %in% c("H. sapiens","M. musculus") &
        !grepl("^ENS", iD_summary$Table[["Feature ID type"]][1] ) ){
      output$cType_mode <- renderUI({ radioButtons(
        "cType_mode", "Annotation mode", choiceNames = c("Automatic (scCATCH)", "Manual" ),
        choiceValues = c(1,2))})
    } else {
      output$cType_mode <- renderUI({ radioButtons(
        "cType_mode", "Annotation mode", choiceNames = c("Manual"),choiceValues = c(2))})
    }
    output$cType_plot <- renderUI({ checkboxInput("cType_plot","Make plots", value = TRUE) })
  } )
  
  observeEvent( input$cType_mode,{
    if( input$cType_mode == 1 ){
      output$cType_mtext <- renderUI({ HTML("<h3>scCATCH options</h3><br>Select a reference:")})
      output$cType_help2.3 <- renderUI({ helpButton( "cType_help2" )})
      output$cType_m1cancer <- renderUI({ radioButtons(
        "cType_m1cancer", "State:", choiceNames = c("Healthy","Cancer"),
        choiceValues = c("healthy","cancer"), selected="healthy"
      )})
      
      cellVars <- sort( colnames( inputDataset$Seurat@meta.data ) )
      
      output$cType_GroupVar <- renderUI( {selectInput( "cType_GroupVar",
        "Grouping variable:", choices =  cellVars, multiple = FALSE )
        
      })
      output$cType_m2ctName <- NULL
      output$cType_m2ctMarkers <- NULL
      output$cType_m2add <- NULL
      output$cType_m2rm <- NULL
      
    } else if(  input$cType_mode == 2 ){
      output$cType_mtext <- renderUI({ HTML("<h3>Manual Annotation</h3><br>Add cell types and makers:")})
      output$cType_help2.3 <- renderUI({ helpButton( "cType_help3" )})
      output$cType_m2ctName <- renderUI({ textInput(
        "cType_m2ctName", "Cell Type Name:", placeholder = "Adipocyte")})
      output$cType_m2ctMarkers <- renderUI({ textInput(
        "cType_m2ctMarkers", "Markers (separated by commas):", placeholder = "ADIPOQ,PPARG,LIPE")})
      
      output$cType_m2add <- renderUI({ actionButton("cType_m2add", HTML("Add")) })
      output$cType_m2rm <- renderUI({ actionButton("cType_m2rm", HTML("Clear list")) })

      if(is.null( inputDataset$Seurat ) & is.null(diffExpr$Parameters) ){
        output$cType_GroupVar <- renderUI({ HTML("No dataset loaded") })
      } else if( !is.null( inputDataset$Seurat ) & is.null(diffExpr$Parameters)  ){
        output$cType_GroupVar <- renderUI({ HTML("Differential expression not tested") })
      } else {
        output$cType_GroupVar <- renderUI({ HTML( paste0(
          "Grouping variable: ", diffExpr$Parameters["Grouping Variable",] ) )})
      }
      output$cType_m1cancer <-NULL
      output$cType_m1tissues <- NULL
      output$cType_m1subtissues <- NULL
    }
  } )
  observeEvent( input$cType_m1cancer,{
    load("data/scCATCH_tissues.RData")
    sctype <- scCATCH_tissues[[ gsub( "\\. ","",iD_summary$Table$Species) ]]
    output$cType_m1tissues <- renderUI({ selectInput(
      "cType_m1tissues", "Tissue type:",  multiple = FALSE,
      choices = names( sctype[[ input$cType_m1cancer  ]]  )  ) })
  } )
  observeEvent( input$cType_m1tissues,{
    load("data/scCATCH_tissues.RData")
    sctype <- scCATCH_tissues[[ gsub( "\\. ","",iD_summary$Table$Species) ]]
    output$cType_m1subtissues <- renderUI({ selectInput(
      "cType_m1subtissues", "Tissue subtype:",  multiple = FALSE,
      choices = sctype[[ input$cType_m1cancer  ]][[ input$cType_m1tissues ]]  ) })
  })
  
  # throw warnings if there are not enough scCATCH markers in the dataset
  observeEvent( input$cType_m1subtissues, {
    cType_checkGenes( inputDataset$Seurat, auto_markers=input$cType_m1subtissues, iDsum=iD_summary$Table )
  }, ignoreNULL = TRUE)
  
  observeEvent( input$cType_m2add,{ # (manual) update marker list  
    cType_markers$Ref <- cType_updInput( input$cType_m2ctName, input$cType_m2ctMarkers, cType_markers$Ref )
  } )
  observeEvent( input$cType_m2rm,{ # (manual) clean marker list
    cType_markers$Ref <- NULL
  } )
  
  observeEvent( cType_markers$Ref, {
    #  throw warnings if the dataset is missing any of the manual markers.
    if( !is.null( inputDataset$Seurat )){
      cType_checkGenes( inputDataset$Seurat, manual_markers = cType_markers$Ref )
    }
    # update marker display
    output$cType_m2DisplayMark <- renderUI({ renderTable({
      cType_displayMarkers( cType_markers$Ref )
    })})
  }, ignoreNULL = TRUE )
  
    ### ct2: run cell type imputation ----
  
  observeEvent( input$cType_run,{ 
    gen_noDataNotif( inputDataset$Seurat )
    if( !is.null( inputDataset$Seurat ) ){
      showNotification( id = "cType_busy",
        HTML( paste0( "Running cell type imputation<br>",
        "This might take a while...")),type = "message",duration = NULL) 
      
      if( input$cType_mode == 1 ){ # scCATCH
        errorCheck <- cType_checkGenes( inputDataset$Seurat, auto_markers=input$cType_m1subtissues, iDsum=iD_summary$Table )
        
        if(errorCheck != "scCATCH_lowMarker"){
          cType_markers$results <- run_cellType( 
          inputDataset$Seurat, 1, input$cType_plot, input$cType_GroupVar, 
          m1_species = iD_summary$Table$Species[1], m1_tissue = input$cType_m1subtissues )
        inputDataset$Seurat <- AddMetaData( inputDataset$Seurat, cType_markers$results$cell.metaData)
        }
        
        
      } else if( input$cType_mode == 2  ){ # manual annotation
        
        if( !is.null(diffExpr$Tables)){
          cType_markers$results <- run_cellType(
            inputDataset$Seurat, 2, input$cType_plot, groupVar = diffExpr$Parameters["Grouping Variable",], 
            m2_markers = cType_markers$Ref, m2_DEparams = diffExpr$Parameters, 
            m2_diffExpr = diffExpr$Tables )
          inputDataset$Seurat <- AddMetaData( inputDataset$Seurat, cType_markers$results$cell.metaData)
        }
      }
    
    }
    
  })
  
    ### ct3: outputs ----
  observeEvent( cType_markers$results, {
    
    if( input$cType_mode == 1){
      output$cType_resText <-renderUI({ HTML( paste0(
        "<h3>Cell type imputation results</h3><br>Cell type summary:"
      )) })
      output$cType_autoResTab <- renderUI({ renderTable(
        {cType_markers$results$imputation}, rownames = FALSE, colnames=TRUE
      )})
      
      output$cType_resText2 <-renderUI({ HTML( paste0(
        "<br>There are ", nrow(cType_markers$results$markerInfo),
        " potential marker genes in CellMatch database for ", 
        iD_summary$Table$Species," on ", input$cType_m1subtissues
      )) })
      
      if(!is.null(cType_markers$results$plots)){
        output$cType_dimplot <- renderUI({ renderPlot({
          grid.arrange( cType_markers$results$plots )
        }, width = 800, height = 600 ) })
      } else {  output$cType_dimplot <- NULL }
      
      output$cType_manualTabC <- NULL
      output$cType_resText3 <- NULL
      output$cType_manualTabG <- NULL
      output$cType_selectGroup <- NULL
      output$cType_ringPlot <- NULL
      
    } else if( input$cType_mode == 2 ){
      output$cType_resText <-renderUI({ HTML( paste0(
        "<h3>Cell type imputation results</h3><br>Cell type summary:"
      )) })
      output$cType_manualTabC <- renderUI({ renderTable(
        {cType_markers$results$imputation }, rownames = FALSE, colnames=TRUE ) })
      
      output$cType_resText3 <- renderUI({ HTML( "<br>Marker summary:" ) })
      
      output$cType_manualTabG <- renderUI({ renderTable(
        {cType_markers$results$markerInfo }, rownames = FALSE, colnames=TRUE ) })
      
      output$cType_selectGroup <- renderUI({ selectInput(
        "cType_selectGroup", "Group:", names( cType_markers$results$plots ), multiple = F )  
      }) 
      output$cType_resText2 <- NULL
      output$cType_autoTab <- NULL
      output$cType_dimplot <- NULL
      
    }
  })
  
  observeEvent( input$cType_selectGroup, {
    if( input$cType_mode == 2 ){
      output$cType_ringPlot <- renderUI({ renderPlot({
        grid.arrange( cType_markers$results$plots[[ input$cType_selectGroup ]] )
      }, width = 1200, height = 600 ) })
    } else { output$cType_ringPlot <- NULL }
    
    
  })
  
  ## Altered Pathways ----
  
  ### p1: UI values ----
  output$paths_currentData <- renderUI({ HTML("No dataset loaded.") })
  output$paths_CPselect <- renderUI({ selectInput(
    "paths_CPselect", label = "",    choices=c("No dataset loaded") 
  )})
  output$paths_DEGselect <- renderUI({ selectInput(
    "paths_DEGselect", label = "",    choices=c("No dataset loaded") 
  )})
    
  observeEvent( c(inputDataset$Seurat, diffExpr$Tables), {
    output$paths_currentData <- gen_Msg$currentData  
    
    if(is.null(diffExpr$Tables)){
      output$paths_CPselect <- renderUI({ selectInput(
        "paths_CPselect", label = "",    choices=c("Differential expression not tested") 
      )})
      output$paths_DEGselect <- renderUI({ selectInput(
        "paths_DEGselect", label = "",    choices=c("Differential expression not tested") 
      )})
    } else {
      output$paths_GroupVar <- renderUI({ HTML( paste0(
        "Differential expression results loaded. <br>Grouping variable: ",
        diffExpr$Parameters["Grouping Variable",] ) )})
      
      # limit by species & ID format
      opts <- grep( gsub( "\\. ","", iD_summary$Table$Species[1]),
                    dir("data/Pathways/" ), value = T )
      if( grepl( "^ENS", iD_summary$Table[1, "Feature ID type"] )){
        opts <- opts[ grep("^ENSEMBL", opts ) ]
      } else { opts <- opts[ grep("^symbols", opts ) ] }
      
      opts <- unique( 
        gsub("[a-zA-Z]*_", "", gsub( "_[a-zA-Z]*\\.txt", "", opts ) ) )
      
      output$paths_CPselect <- renderUI({ selectInput( 
        "paths_CPselect", label = "",  multiple = TRUE, selected = "All",
        choices=c("All", opts ), 
        ) })
      output$paths_DEGselect <- renderUI({ selectInput( 
        "paths_DEGselect", label = "",  multiple = TRUE, selected = "All",
        choices=c("All", names(diffExpr$Tables)), 
      ) })
      
    }
  })
  
  ### p2: run pathway enrichment ----
   observeEvent( input$run_paths,{
     
     gen_noDataNotif( inputDataset$Seurat )  # Pop up notifs if missing data
     if(is.null(diffExpr$Tables)){
       showNotification( "Differential expression not tested",
         duration = 20, type="error", id="gen_noDEGs" )
     }
     
    if( !is.null(inputDataset$Seurat) & !is.null(diffExpr$Tables)){
      
      showNotification( id = "paths_busy",
        HTML( paste0( "Calculating Enrichment in metabolical pathways.<br>",
          "This might take a while...")),type = "message",duration = NULL) 
      
      mSig_Paths <- paths_loadDBs( iD_summary$Table, input$paths_CPselect ) # load gene sets and subset DE results
      DEGtabs <- paths_selectDEresult( diffExpr$Tables, input$paths_DEGselect )
     
      ap <- run_paths( inputDataset$Seurat, DEtabs = DEGtabs, pathDB = mSig_Paths,
                                 DEparams = diffExpr$Parameters )
      
      removeNotification( id="paths_busy" )
      alteredPaths$Tables <- ap$Tables
      alteredPaths$Summary <- ap$Summary
      alteredPaths$backGrounds <- ap$backGrounds
      
      
     }
   })
  
  ### p3: outputs & plot options ----
  
  observeEvent( alteredPaths$Tables,{

    if( !all( is.null( alteredPaths$Tables ))){
      
      output$paths_resText1 <- renderUI(HTML(paste0( 
        "<h3>Metabolical pathways altered by differential expression</h3>",
        "<br>Results summary:<br>")) )
      output$paths_summary <- renderUI({renderTable({
        tab <- alteredPaths$Summary
        colnames(tab) <- gsub("\\."," ", colnames(tab))
        tab
        })  })
      output$paths_resText2 <- renderUI({ HTML(
        paste0( "Plot the DE state of genes (up/down regulation,",
        "no change, etc) expressed in the analized groups of cells ",
        "that belong to particular pathways.<br><br><br>"))})
      output$paths_help2 <- renderUI( helpButton("paths_help2") )
      output$paths_resplot_text1 <- renderUI({ HTML("Select DE analysis results:")})
      
      # initial state of the plot selectors
      output$paths_resPlot_sel1 <- renderUI({ selectInput(
        "paths_resPlot_sel1", NULL, c( "All", names(alteredPaths$Tables) ), 
        selected = "All", multiple=T ) })

      output$run_pathsPlot <- renderUI({ actionButton(
        "run_pathsPlot", "Generate plots",
        style="color: #fff; background-color: #3c8dbc; border-color: #015e95")
      })
    } else {
      output$paths_resText1 <- NULL
      output$paths_summary <- NULL
      output$paths_resPlot_sel1 <- NULL
      output$paths_resPlot_sel2 <- NULL
      output$paths_resPlot <- NULL
    }
  } )
  
  observeEvent( c( input$paths_resPlot_sel1 ),{
    output$paths_resplot_text2 <- renderUI({ HTML("Select pathway gene sets:")})
    
    if( all( input$paths_resPlot_sel1 =="All" ) | length(input$paths_resPlot_sel1)>1 ){
      pathList <- unique(unlist( sapply( alteredPaths$Tables, function(i){i$list} ) ) )
      output$paths_resPlot_sel2 <- renderUI({ selectizeInput(
        "paths_resPlot_sel2", NULL, NULL, multiple = F ) })
      updateSelectizeInput(
        session=session, inputId = "paths_resPlot_sel2",
        choices = pathList, server=T )
    } else {
      pathList <- unique(unlist( sapply( alteredPaths$Tables, function(i){i$list} ) ) )
      output$paths_resPlot_sel2 <- renderUI({ selectizeInput(
        "paths_resPlot_sel2", NULL, NULL, multiple = T ) })
      updateSelectizeInput(
        session=session, inputId = "paths_resPlot_sel2",
        choices = pathList, server=T )
      
    }
  } )
  observeEvent( input$run_pathsPlot,{
    if(!is.null(input$paths_resPlot_sel1) | !is.null(input$paths_resPlot_sel2)){
      pathPlot <- paths_likertPlots( 
        diffExpr$Tables, diffExpr$Parameters,
        alteredPaths$backGrounds, alteredPaths$Tables,
        plotComp = input$paths_resPlot_sel1, plotPath = input$paths_resPlot_sel2 )
      if( input$paths_resPlot_sel1 =="All" ){
        plW <- 800  
        plH <- 250 + 30*length( alteredPaths$Tables )
      } else {
        plW <- 1200  
        if( !is.null( input$paths_resPlot_sel2)){
          plH <- 250 + 30*length( input$paths_resPlot_sel2 )
        } else {
          plH <- 550
        }
      }
      output$paths_resPlot <- renderUI({ renderPlot(
        pathPlot, width = plW, height = plH
      )})
    }
    
  } )
  
  ## SAVE AND EXPORT FILES ----
  
  observeEvent( inputDataset$Seurat,{
    # cell meta data (CSV)
    output$Save_QC_md <- gen_save( inputData = inputDataset$Seurat, type="meta" )
    output$Save_Filter_md <- gen_save( inputData =  inputDataset$Seurat, type="meta" )
    output$Save_Merge_md <- gen_save( inputData = inputDataset$Seurat, type="meta" )
    output$Save_cluster_md <- gen_save( inputData = inputDataset$Seurat, type="meta" )
    output$Save_featPl_md <- gen_save( inputData = inputDataset$Seurat, type="meta" )
    output$Save_cType_md <- gen_save( inputData = inputDataset$Seurat, type="meta" )
    # Seurat object (RDS)
    output$Save_QC_RDS <- gen_save( inputData = inputDataset$Seurat, type="seuratObject" )
    output$Save_Filter_RDS <- gen_save( inputData =  inputDataset$Seurat, type="seuratObject" )
    output$Save_Merge_RDS <- gen_save( inputData = inputDataset$Seurat, type="seuratObject" )
    output$Save_cluster_RDS <- gen_save( inputData = inputDataset$Seurat, type="seuratObject" )
    output$Save_Visu_RDS <- gen_save( inputData = inputDataset$Seurat, type="seuratObject" )
    output$Save_featPl_RDS <- gen_save( inputData = inputDataset$Seurat, type="seuratObject" )
    output$Save_cType_RDS <- gen_save( inputData = inputDataset$Seurat, type="seuratObject" )
  } )
  
  # Export miscellaneous tables (csv & excel/openxsls)
  observeEvent( expCluster$nPCs, {
    output$Save_ExCl_nPCs <- gen_save(
      tabs= expCluster$nPCs$nPC$clusteringTab, 
      type="tables",
      f_name = paste0( inputDataset$Seurat@project.name,"_clustering_by_number_of_PCs") )
  }, ignoreInit = TRUE )
  
  observeEvent( expCluster$clResolution, {
    output$Save_ExCl_res <- gen_save( 
      tabs =  expCluster$clResolution$Resolution$clusteringTab, 
      type="tables",
      f_name = paste0( inputDataset$Seurat@project.name,"_clustering_by_resolution"))
    
    output$Save_ExCl_unCert <- gen_save( 
      tabs = expCluster$clResolution$uncertTab,
      type="tables",
      f_name = paste0( inputDataset$Seurat@project.name,"_clustering_uncertainty") )
  }, ignoreInit = TRUE )
  
  
  observeEvent( Visu_data$Results,{
    
    output$Save_Visu_reds <-  gen_save( tabs=Visu_data$Results$redCoords, type="tables",
      f_name = paste0( inputDataset$Seurat@project.name,"_Reduction_coordinates" ) )
  }, ignoreInit = TRUE )
  
  observeEvent( diffExpr$Tables,{
    degTabs <- c(                     # save with param table
      list(as.data.frame(diffExpr$Parameters, row.names = rownames( diffExpr$Parameters) )), 
      diffExpr$Tables ) 
    names(degTabs) <- c("Parameters", names(diffExpr$Tables) )
    output$Save_DEG_tabs <-  gen_save( tabs=degTabs , type="tables",
      f_name =  paste0( inputDataset$Seurat@project.name,"_differential_expression") )
  }, ignoreInit = TRUE )
  
  observeEvent( alteredPaths$Tables,{
    pathwayTabs <- c(                     # save with param table
      list(as.data.frame(diffExpr$Parameters, row.names = rownames( diffExpr$Parameters) )), 
      alteredPaths$Tables ) 
    names(pathwayTabs) <- c("DE_Parameters", names(alteredPaths$Tables) )
    output$Save_path <-  gen_save( tabs=pathwayTabs , type="tables",
      f_name =  paste0( inputDataset$Seurat@project.name,"_altered_pathways") )
  }, ignoreInit = TRUE )
  
  
  observeEvent( cType_markers$results, {
    output$Save_cType_cellType <- gen_save( tabs=cType_markers$results$imputation, type="tables",
      f_name = paste0( inputDataset$Seurat@project.name, "_cell_type_imputation") )
    output$Save_cType_markerInfo <- gen_save( tabs=cType_markers$results$markerInfo, type="tables",
      f_name = paste0( inputDataset$Seurat@project.name, "_cell_type_markerInformation") )
  }, ignoreInit = TRUE )
  
  
  ## Help Messages ----
  
  observe({
    
    btnState <- updButtonStates( input, "_help", helpBtnStates$bSts )
    
    helpBtnStates$bSts <- btnState[["bSts"]]
    helpBtnStates$btn <- btnState[["btn"]]
    
  })
  
  observeEvent( helpBtnStates$bSts, {
    if(!is.null( helpBtnStates$btn )){
      section <- gsub( "_.*$", "", helpBtnStates$btn )
      output[[paste0( section, "_help")]] <- getMssg( btn=helpBtnStates$btn )
    }
  })
  
  output$Load_dwlEx <- downloadHandler(
    filename = "Example_datasets.zip",
    content = function(file) {
      file.copy("data/Example_Datasets.zip", file) },
    contentType = "application/zip"
  )
}





shinyApp(ui, server)















