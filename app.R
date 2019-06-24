#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library( shiny)
library( CTDquerier)
library( httr)
library( GEOmetadb)
library( plyr)
library( tidyverse)
library( DT)
library( xml2)
library( stringr)
library( reshape2)
library( shinycssloaders)
library( ggpubr)

## source some code
source( "shiny_metabolomics.R", local = TRUE)
source( "shiny_proteomics.R", local = TRUE)
source( "shiny_transcriptomics.R", local = TRUE)
source( "shiny_general.R", local = TRUE)
source( "shiny_compound.R", local = TRUE)
source( "shiny_output_texts.R", local = TRUE)
source( "plot_CTD_information.R", local = TRUE)

######
# CONSTANTS
pubchem_compound <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
pubmed_url <- "https://www.ncbi.nlm.nih.gov/pubmed/"
geo_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="


######
# load the compTox database
######
load( file = "data/comptox.RData")
load( file = "data/CTD_chemicals.RData")

# to query proteomexchange database
#library( rpx)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   #titlePanel("Search for multi-omics data sets in various public sources."),
  
  titlePanel(
     fluidRow( column( 6, h1( "MOD-Finder"), HTML( "<h3>Search for <b>M</b>ulti-<b>O</b>mics <b>D</b>ata sets in various public sources.</h3>")),
               column( 2, 
                       tags$style(".topimg {
                            padding-left:10px;
                            padding-right:10px;
                          }"),
                         div( class="topimg", img( height = 110, src="cefic_logo.jpg", class="pull-right"))),
               column( 2, 
                       tags$style(".topimgtwo {
                            margin-left:50px;
                          }"),
               div( class="topimgtwo", img( height = 110, src="LRI_logo.jpg"))),
               column( 2, 
                       tags$style(".topimgthree {
                            padding-left:10px;
                            padding-right:10px;
                          }"),
                div( class="topimgthree", img( height = 120, src="UFZ_logo.jpg", class="pull-left")))
               ),
     windowTitle = "MOD-Finder"
   ) ,
   
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     
      sidebarPanel(
        width = 3,
        
        h4( "Step I.a - Enter Compound"),
        
        textInput( inputId = "compoundName", label = "Compound Name or ID", placeholder = "Enter Compound", width = '400px'),
        
        #br(),
        #hr(),
        br(),
        
        h4( "Step I.b - Select Outputs" ),

        checkboxGroupInput( inputId = "addInformation", choices = c( "Integrated chemical information (CTDbase)" = "CTDbase"), selected = c("CTDbase"), label = "Additional Information" ),
        
        checkboxGroupInput( inputId = "omicsLayer", choices = c( "Transcriptome" = "trans", "Proteome" = "prot", "Metabolome" = "meta"), label = "Omics Layer" ),
        
        actionButton( inputId = "refineCompound", label = "Refine"),
      
        br(),
        br(),
        hr(),
        br(),
        
        h4( "Step II - Filter the exact Compound"),
        
        selectInput( inputId = "compoundList", choices = NULL,
                     label = "Specify the compound to search for:"),
        
        radioButtons( inputId = "splitCompound", 
                    label = "Name splitting", choices = c( "Search with refined compound")),
        
        actionButton( inputId = "startSearch", label = "Search")
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        
        
        tabsetPanel( id = "tabset",
                     tabPanel( value = "about",
                             title = "About",
                             h1("Welcome to MOD-Finder"),
                             p( "MOD-Finder is an easy-to-use R Shiny tool to search for compound-related omics datasets in different layer: transcriptome, proteome, metabolome."),
                             br(),
                             h3( "MOD-Finder Workflow"),
                             p( "Given a certain chemical name or ID, MOD-Finder uses the CompTox Dashboard and Pubchem databases to digitally identify the compound. In a next step, omics data sets, that are somehow linked with this compound are searched in public databases. Therefore, the user can choose which specific omics layer are targeted, e.g., transcriptomics, proteomics, and/or metabolomics. Furthermore, compound-specific information can be retrieved from additional sources (CTDbase) and will be visualized to provide insights into the perturbations that are known to be triggered by a specific chemical."),
                             br(),
                             h3( "Public omics databases that can be queried by MOD-Finder:"),
                             p( fluidRow( column( 3, a(img( src = "geo_logo.jpg", height = 80), href = "https://www.ncbi.nlm.nih.gov/geo/" )),
                                       column( 3, a( img( src = "AE_logo.png", height = 80), href = "https://www.ebi.ac.uk/arrayexpress/")),
                                       column( 3, a( img( src = "pride_transparent.png", height = 80), href ="https://www.ebi.ac.uk/pride/archive/"))
                             ),
                             br(), 
                             fluidRow( column( 3, a( img( src = "metabolights_logo.jpg", height = 80), href = "https://www.ebi.ac.uk/metabolights/")),
                                       column( 3, a( img( src = "met_workbench_logo.jpg", height = 80), href = "https://www.metabolomicsworkbench.org/")),
                                       column( 3, a( img( src = "MeRy-B2_logo.png", height = 80), href = "http://services.cbib.u-bordeaux.fr/MERYB/")),
                                       column( 3, a( img( src = "metabolonote_logo135px.png", height = 80), href = "http://metabolonote.kazusa.or.jp/Main_Page"))
                             )),
                             br(),
                             br(),
                             h3( "How to use MOD-Finder"),
                             p( h4( "1) Specify the compound of interest"),
                                p( "- In ", strong("Step I.a"), " insert either the compound name, a synonym, or the ID (Pubchem, CAS, or Comptox Dashboard)"),
                                h4( "2) Select the ouputs"),
                                p("- Specify if additional compound information should be retrieved, integrated, and visualized"),
                                p("- Specify in which omics layer data sets should be searched for"),
                                p("- Click the 'Refine' buttion"),
                                h4( "3) Refine the exact search string(s)"),
                                p( "- When a compound of interest results in multiple potential search strings, use the provided options in ", strong( "Step II"), " to specify either the exact search string that should be used for the hunt for omics data sets or search with all combinations." ),
                                p( "- Click the 'Search' button")
                             ),
                             br(),
                             h3( "Further Reading"),
                             p( "For more detailed information, please have a look at the Application Note: Canzler ", em("et al."), ", 2019, ", em("MOD-Finder: Identify mulit-omics data sets related to defined chemical exposure"), " (submitted)"),
                             br(),
                             h3( "Funding"),
                             p( "This work was funded by the Cefic Long-Range Research Initiative Programme (Project ", a( "C5-XomeTox", href = "http://cefic-lri.org/projects/c5-xometox-evaluating-multi-omics-integration-for-assessing-rodent-thyroid-toxicity/"), ").")
                            
                    )
        )
  )
)
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  chem <- reactiveValues( searchCompound = FALSE)

  observeEvent( input$refineCompound, {
    
    chem$searchCompound <- FALSE

    ## remove tabs displaying search results in case 'search' was hit again
    ## otherwise, the old tabs will also be displayed
    if( !is.null(chem$displayedTabs)){
      chem$displayedTabs <- removeTabs( session, chem$displayedTabs)
    }
    
    ## get the compound to search for
    chem$refineCompound <- input$refineCompound
    chem$compound <- isolate( input$compoundName)

    withProgress( message = "Search for compounds.", value = 0, {

      ## first of all, search comptox dashboard for compounds
      incProgress( amount = 0.5, detail = "Query comptox dashboard.")
      compound_list <- get_comptox_by_input( chem$compound)
      
      ## if an emtpy result is returned, search pubchem
      if( compound_list$name[1] == "EMPTY"){
       
        incProgress( amount = 0.2, detail = "Query pubchem database.")
        compound_list <- get_cids_by_name( chem$compound)
        
      }
 
      ## remove "NA" from compound_list$name
      nas <- which( is.na(compound_list$name)) 

      if( length( nas) > 0){
        compound_list <- compound_list[-nas,]
        if( nrow( compound_list) == 0){
          compound_list <- data.frame( cid=c(0000), name=c("EMPTY"))
        }
      }

      ## display the pubchemIDs with their associated names, in case there are entries left in compound_list
      ## otherwise display "No compound found!"
      if( compound_list$name[1] == "EMPTY"){

        showModal( modalDialog(
          title = "Compound ERROR",
          "The entered search string lead to NO RESULTS!",
          easyClose = TRUE
        ))
        
      }
      
      if( compound_list$name[1] != "CAS" && compound_list$name[1] != "EMPTY"){
        
        list_input <- paste( compound_list$cid, compound_list$name, sep = " - ")
          updateSelectInput( session, inputId = "compoundList",
                             choices = list_input
                             )
          
      } else{

        updateSelectInput( session = session,
                           inputId = "compoundList",
                           choices = "No compound found!"
                           )
      }
    })
  })
  
  
  observeEvent( input$compoundList, {
    
    if( input$compoundList != ""){
    
      if( input$compoundList == "No compound found!"){
      
        updateRadioButtons( inputId = "splitCompound", session = session,
                            label = "Search options: ",
                            choices = c( "Don't search", "Search with original input compound"))
      
      } else {
      
        inputCompound <- tolower(isolate( input$compoundName))
        cmp <- paste( unlist( strsplit( x = tolower( input$compoundList), split = " - "))[-1], collapse = "-")
      
        # in case refined compound is the same as the initial input

        if( cmp == inputCompound ){
          updateRadioButtons( inputId = "splitCompound", session = session,
                              label = "Name splitting - Search with:",
                              choiceNames = c( paste0( cmp, " (refined input)")),
                              choiceValues = c( cmp))
        } else {
    
          cN <- c( paste0( inputCompound , " (initial input)"), paste0( cmp, " (refined input)"))
          cV <- c( inputCompound, cmp)
          
          cmp_list <- unlist( strsplit( x = cmp, split = "-"))
          
          if( length( cmp_list) > 1){
            
            for( i in seq( 2, length( cmp_list))){ 
              n <- paste( cmp_list[i:length( cmp_list)], collapse = "-")
              if( !n %in% cV){
                cV <- c( cV, n)
                cN <- c( cN, n)
              }
            }
          }
          
          cN <- c( cN, "all choices")
          cV <- c( cV, "all")
          
          updateRadioButtons( inputId = "splitCompound", session = session,
                              label = "Name splitting - Search with:",
                              choiceNames = cN,
                              choiceValues = cV)
        }
      }
    }
    
  })
  
  observeEvent( input$startSearch, {
    
    ## first, remove existing tabs in main window
    if( !is.null( chem$displayedTabs)){
      chem$displayedTabs <- removeTabs( session, chem$displayedTabs)
    }
    
    chem$searchCompound <- TRUE
    chem$selectedCompound <- input$compoundList

    chem$exactCompoundList <- NULL
    chem$results <- NULL
    chem$error <- NULL
    chem$ids <- NULL
    
    if( chem$selectedCompound == "No compound found!" & input$splitCompound == "Don't search") { 
    
      chem$error <- "ERROR: No compound found!"
      chem$exactCompound <- NULL
      chem$exactCompoundList <- NULL
      chem$cid <- NULL

    } else{
      
      ## in case there is no pubchem ID but the user wants to search with the initial compound name
      ## be aware that there is no specified cid yet!
      if( chem$selectedCompound == "No compound found!"){
      
        chem$exactCompound <- input$compoundName
        chem$exactCompoundList <- c( input$compoundName)
        chem$cid <- NULL
              
      } else{
    
        ## the selected compound from the selectList has the format: cid - name
        ## hence, we have to split the string at the first '-'
        string_elements <- unlist( strsplit( chem$selectedCompound, " - ", fixed = TRUE))
        chem$cid <- string_elements[1]
        chemical_name <- paste( string_elements[ -1], collapse = " - ")
        chem$exactCompound <- paste( string_elements[ -1], collapse = " - ")
        
        if( input$splitCompound == "all"){
          
          cmp <- unlist( strsplit( x = tolower( input$compoundList), split = " - "))[-1]
          cmp_list <- unlist( strsplit( x = cmp, split = "-"))
          
          chem$exactCompoundList <- c( input$compoundName, cmp)
          
          if( length( cmp_list) > 1){
            
            for( i in seq( 2, length( cmp_list))){ 
              n <- paste( cmp_list[i:length( cmp_list)], collapse = "-")
              if( !n %in% chem$exactCompoundList){
                chem$exactCompoundList <- c( chem$exactCompoundList, n)
              }
            }
          }
          
        } else {
          
          chem$exactCompoundList <- c( input$splitCompound)

        }
      }
      
      
      ## in case the compound List contains parentheses
      ## add search strings with exchanged parentheses, e.g., ( -> [
      m <- grep( "\\(", chem$exactCompoundList)
      changed1 <- gsub( "\\((.*)\\)", "[\\1]", chem$exactCompoundList[m])
      m <- grep( "\\[", chem$exactCompoundList)
      changed2 <- gsub( "\\[(.*)\\]", "(\\1)", chem$exactCompoundList[m])
      chem$exactCompoundList <- c( chem$exactCompoundList, changed1, changed2)
      
      
      
      ## initialize a progress bar showing information about the current step
      withProgress( message="Collecting compound information", value=0, {
      
        if( !is.null( chem$cid)){ 
          
        ## collect information from pubchem using PUG REST
        ## first: CID and synonyms
        incProgress( amount = 0.2, detail = "Collecting Pubchem Information")
          
          
        ## Retrieve a list of synonyms through a RESTful HTTP request
        ## by means of the CID
        pubchem <- tryCatch( 
          
          get_synonyms_by_cid( chem$cid),
            
          warning = function(cond){
            message( "Querying Pubchem resulted in a warning.")
            return( list())
          },
            error = function(cond){
            message( "Querying CTDbase resulted in an error.")
            return( list())
          })

        
        chem$synonyms <- unique( tolower( pubchem$synonyms))

        
        ## second: general information, such as SMILE and descriptions
        chem$results <- tryCatch(
          
          get_compound_information( chem$cid),

          warning = function(cond){
            message( "Querying Pubchem resulted in a warning.")
            return( data.frame())
          },
          error = function(cond){
            message( "Querying CTDbase resulted in an error.")
            return( data.frame())
        })
          
        ## get pubchem ID mappings
        pubchem <- tryCatch(
          
          get_description_by_cid( chem$cid),
          
          warning = function(cond){
            message( "Querying Pubchem resulted in a warning.")
            return( list())
          },
          error = function(cond){
            message( "Querying CTDbase resulted in an error.")
            return( list())
        })
        

        ## create a data frame which stores all IDs collected during the runtime
        chem$ids <- data.frame( Database = c( "Pubchem"), IDs = c( sprintf( '<a href=\"https://pubchem.ncbi.nlm.nih.gov/compound/%s\" target=\"_blank\">%s</a>', chem$cid, chem$cid)), stringsAsFactors = FALSE)
       
        ## if the synonyms list has more than 10 entries, show only the first ten
        ## keep the complete list in chem$synonyms
        ## and save the output list as a single string in chem$synonyms_list
        if( length( chem$synonyms > 10)){
          chem$synonyms_list <- paste0( paste( chem$synonyms[1:10], collapse = " | "), " ... (10 of ", length( chem$synonyms), ")" )
        } else{
          chem$synonyms_list <- paste( chem$synonyms, collapse = " | ")
        }

        
        incProgress( amount = 0.2, detail = "Collecting CompTox Information")
        comptox <- get_comptox_ids_by_cid( chem$cid)
        
        
        ## combine the information into two tables
        ## one focusing on the chemical properties and descriptions
        ## and the other focusing on the IDs conversions and mappings
        chem$results <- rbind( "Chemical Name" = chemical_name, chem$results)
        chem$results <- rbind( chem$results, "Synonyms" = chem$synonyms_list)
        chem$ids <- rbind( chem$ids, comptox)
        chem$ids <- rbind( chem$ids, data.frame( Database = pubchem$sources, IDs = pubchem$functional_link))

        
        } else{
          
          chem$error <- "No compound information found for initial input name!"
          
        }
      

        
        
        
        ## check which output layer is checked and add the respective tab in the panel
        ## either search for additional information at CTDbase or
        ## search for data sets only in case the respective omics layer is selected
        progress_time = 0.6 / ( length( input$omicsLayer) + length( input$addInformation))
        
        
        if( "meta" %in% input$omicsLayer){
          
          ## collect potential metabolome data sets
          ## use ?? to search Metabolights and Metabolomics Workbench
          incProgress( amount = progress_time, detail = "Collecting Metabolome Data Sets")
          chem$metabolome <- get_metabolomeXchange_by_list( chem$exactCompoundList)
          
          ## chem$metabolome is a list that contains [1]$msg a character vector of errors that occured
          ## and a [2]$df containing the data set table
          
          ## in case there was an error
          amount <- as.character( nrow( chem$metabolome[2]$df))
          if( length( chem$metabolome[1]$msg) > 0 && nrow( chem$metabolome[2]$df) == 0){
            amount <- "ERROR"
          }

          ## create a dynamic tab to display the previously gathered information about metabolomics
          prependTab( inputId = "tabset",
                      session = session,
                      tabPanel( value = "Metabolome",
                                title = paste0( "Metabolome (", amount, ")"),
                                h3(textOutput( outputId = "headerMeta")),
                                htmlOutput( outputId = "headerMetaError"),
                                br(),
                                DT::dataTableOutput( outputId = "tableMeta"))
          )
          
        }
        

        if( "prot" %in% input$omicsLayer){
          
          ## collect potential proteome data sets
          ## use the rpx package to search amongst proteinXchange-uploaded data sets
          incProgress( amount = progress_time, detail = "Collecting Proteome Data Sets")
          chem$proteome <- get_proteome_by_list( chem$exactCompoundList)

          ## chem$proteome is a list that contains [1]$msg a character vector of errors that occured
          ## and a [2]$df containing the data set table

          amount <- as.character( nrow( chem$proteome[2]$df))
          if( length( chem$proteome[1]$msg) > 0 && nrow( chem$proteome[2]$df) == 0){
            amount <- "ERROR"
          }
                              
          ## create a dynamic tab to display the previously gathered information about proteomics
          prependTab( inputId = "tabset",
                     session = session,
                     tabPanel( value = "Proteome",
                               title = paste0( "Proteome (", amount, ")"),
                               h3(textOutput( outputId = "headerProt")),
                               htmlOutput( outputId = "headerProtError"),
                               br(),
                               DT::dataTableOutput( outputId = "tableProt"))
          )
          
        }
        

        
        if( "trans" %in% input$omicsLayer ){
          
          ## collect potential transcriptome data sets
          ## use GeoMetadb to search for uploaded data sets were the title contains the compound name
          incProgress( amount = progress_time, detail = "Collecting Transcriptome Data Sets")
          chem$transcriptome <- get_transcriptome_by_list( chem$exactCompoundList)
          chem$transcriptomeAE <- get_arrayExpress_by_list( chem$exactCompoundList)
          
          
          ## chem$transcriptome and chem$transcriptomeAE are lists that contains [1]$msg a character vector of errors that occured
          ## and a [2]$df containing the data set table
          
          ## in case there was an error
          amount <- nrow( chem$transcriptome[2]$df) + nrow( chem$transcriptomeAE[2]$df)
          if( amount == 0){
            if( nrow( chem$transcriptome[2]$df) + nrow( chem$transcriptomeAE[2]$df) == 0)
              amount <- "ERROR"
          }
          amount <- as.character( amount)
          
          ## create a dynamic tab to display the previously gathered information about transcriptomics
          prependTab( inputId = "tabset",
                     session = session,
                     tabPanel( value = "Transcriptome",
                               title = paste0( "Transcriptome (", amount, ")"),
                               h3( textOutput( outputId = "headerTrans")),
                               htmlOutput( outputId = "headerTransError"),
                               DT::dataTableOutput( outputId = "tableTrans"),
                               br(),
                               h3( textOutput( outputId = "headerTransAE")),
                               htmlOutput( outputId = "headerTransAEError"),
                               DT::dataTableOutput( outputId = "tableTransAE"))
          )
          
        }
        
        
        if( "CTDbase" %in% input$addInformation ){
          
          ## use the CTDbase to gather additional information of the selected compound
          incProgress( amount = progress_time, detail = "Collecting CTDbase information")
          chem$plot <- FALSE
          
          ## if CAS number or MESH ID is present, query CTD base to generate some informational plots
          ## gather CTDbase information on chemical-gene interactions and plot the top40 genes affected by the analyzed compound
          ## query CTDbase with chemical name
          ##
          ## FIRST: try to use the CAS-RN to extract the compound name from the CTD file
          if( "CAS-RN" %in% chem$ids$Database ){
            
            match <- which( chem$ids$Database == "CAS-RN")
            chem$ctd_cas <- chem$ids$IDs[ match]
            chem$ctd_chemical <- ctd$Chemical.Name[ which( ctd$CasRN == chem$ctd_cas)]

          }
          ## 
          ## SECOND: if that fails, try to search with the Mesh ID, if present
          if( length(chem$ctd_chemical) == 0 && "MeSH" %in% chem$ids$Database ){
            
            chem$ctd_cas <- NULL
            match <- which( chem$ids$Database == "MeSH")
            chem$ctd_mesh <- chem$ids$IDs[ match]

            chem$ctd_chemical <- ctd$Chemical.Name[ which( ctd$CasRN == chem$ctd_mesh)]
            
          }
          ##
          ## THIRD: if that also fails, check if the exact compound name is already given
          ## therefore, apply a simple string search
          if( length(chem$ctd_chemical) == 0){
            
            chem$ctd_cas <- NULL
            chem$ctd_chemical <- ctd$Chemical.Name[ which( ctd$Chemical.Name == chem$exactCompound)]
            
          }
          

          ## when a compound name was retrieved from the CTD base file
          ## run the query_ctd_chem command, to gather all necessary information
          ## that is needed for plotting
          if( length(chem$ctd_chemical) != 0 ){            

            chem$ctd_chem <- tryCatch( 
              query_ctd_chem( chem$ctd_chemical, filename = "data/CTD_chemicals.tsv.gz"),
              warning = function(cond){
                message( "Querying CTDbase resulted in a warning.")
                return( NULL)
              },
              error = function(cond){
                message( "Querying CTDbase resulted in an error.")
                return( NULL)
              })

            chem$plot <- TRUE

          }


          if( chem$plot == TRUE){
            ## create a dynamic tab to display the previously CTDbase-gathered information about chemical-gene interactions, pathways, and diseases
            prependTab( inputId = "tabset",
                       session = session,
                       tabPanel( value = "CTDbase",
                                 title = paste0( "Additional Information"),
                                 h3( paste0( "Information on gene interactions through exposure with ", chem$ctd_chemical)),
                                 p( get_generalInformationDescription( chem$ctd_chemical)),
                                 withSpinner( plotOutput( outputId = "plotGeneralInformationGene")),
                                 br(),
                                 br(),
                                 h3( paste0( "Gene expression affected by ", chem$ctd_chemical)), 
                                 p( get_chemicalGeneDescription( chem$ctd_chemical)), 
                                 withSpinner( plotOutput( outputId = "plotChemicalGene", height = "600px")),
                                 br(),
                                 br(),
                                 h3( paste0( "Diseases that are associated with ", chem$ctd_chemical)), 
                                 p( get_diseasesDescription( chem$ctd_chemical)),
                                 withSpinner( plotOutput( outputId = "plotDiseases", height = "700px")),
                                 br(),
                                 br(),
                                 h3( paste0( "Pathway enrichment caused by ", chem$ctd_chemical)), 
                                 p( get_pathwaysDescription( chem$ctd_chemical)),
                                 withSpinner(plotOutput( outputId = "plotPathways", height = "600px"))
                                 
                       )
                     
            )
          }else{
            
            ## in case there are no suitable plots,
            ## the tab will still be displayed, but a message for the user will be given
            ## furthermore, this will circumvent issues, when tabs are removed later on,
            ## e.g., when an additional compound will be searched
            prependTab( inputId = "tabset",
                       session = session,
                       tabPanel( value = "CTDbase",
                                 title = paste0( "Additional Information"),
                                 h3( "No CTD-based information could be retrieved...")
                       ))
          }
        }
        
        
        
        ## add the compound information panel
        prependTab( inputId = "tabset",
                    session = session,
                    select = TRUE,
                    tabPanel( value = "Compound",
                              title = "Compound Information",
                              h3(textOutput( outputId = "headerCompound")),
                              tableOutput( outputId = "tableCompound"),
                              br(),
                              h3( textOutput( outputId = "headerCompoundIDs")),
                              DT::dataTableOutput( outputId = "tableCompoundIDs"))
        ) 
        
        

      })

    }
    
    chem$displayedTabs <- c( input$addInformation, input$omicsLayer)
    
  })
  
  
  ## Render the header of the compound information tab
  ## if errors occurred, display an information 
  output$headerCompound <- renderPrint( {
    
    if( chem$searchCompound == FALSE) return()

    if( "error" %in% names(chem) & !is.null( chem$error)){
      cat( "An ERROR occurred: \n", chem$error, "\n\nPlease reload this page and try again!\n")
    } else if( "exactCompound" %in% names(chem) & !is.null( chem$exactCompound)){
      cat( "Compound Information: ", chem$exactCompound)
    }else{
      cat( "Something unpredictable happend.\n\nPlease reload this page and try again!\n")
    }

  })
  
  
  ## Render the table that collects compound-related information
  ## such as SMILES, InChI codes, etc...
  output$tableCompound <- renderTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    if( !is.null( chem$results)){
      head( chem$results, n = nrow( chem$results))}},
  
      rownames = TRUE,
      colnames = FALSE
  )
  
  
  ## Render the header of the ID conversion table
  output$headerCompoundIDs <- renderPrint( {
    
    if( chem$searchCompound == FALSE) return()
    if( !is.null( chem$ids)){
      cat( "ID conversion and mapping")
    }
    
  })
  
  
  ## Render the data table of the ID conversion
  output$tableCompoundIDs <- DT::renderDataTable({
    
   if( chem$searchCompound == FALSE ) return()
    
    if( !is.null( chem$ids)){
      head( chem$ids, n = nrow( chem$ids))}
    },

    rownames = FALSE,
    escape = FALSE
  
  )

  
  ## Render the plot displaying InteractionActions of the queried compound
  output$plotGeneralInformationGene <- renderPlot({
    
    if( chem$searchCompound == FALSE) return()    
    if( chem$plot){
      plot_interaction_actions( chem$ctd_chem, chem$ctd_chemical, chem$ctd_cas)
    }
    
  })
  
  
  ## Render the bar plot displaying the affects of the queried compound on gene expression
  output$plotChemicalGene <- renderPlot({
    
    if( chem$searchCompound == FALSE) return()    
    if( chem$plot){
      plot_chemical_gene_interaction( chem$ctd_chem, chem$ctd_chemical)
    }
      
  })
  
  
  ## Render the bar plot that links the queried compound to certain diseases
  output$plotDiseases <- renderPlot({
    
    if( chem$searchCompound == FALSE) return()    
    if( chem$plot){
      plot_diseases( chem$ctd_chem, chem$ctd_chemical, chem$ctd_cas)
    }
    
  })
  
  
  ## Render the bar plot that displays the results of a KEGG and REACTOME pathway enrichment
  output$plotPathways <- renderPlot({
    
    if( chem$searchCompound == FALSE) return()    
    if( chem$plot){
      plot_pathways( chem$ctd_chem, chem$ctd_chemical, chem$ctd_cas)
    }
  })
  
  
  ## Render the header of the transcriptome NCBI GEO results
  output$headerTrans <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()

    cat( "NCBI Geo data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })

  
  ## Render the information panel of NCBI GEO errors that might have occurred
  output$headerTransError <- renderPrint( {

    if( chem$searchCompound == FALSE ) return()
    if( length( chem$transcriptome[1]$msg) > 0 ){
      cat( "Please note that error(s) occurred during data set collection:<br>")
      cat( paste(chem$transcriptome[1]$msg, collapse = "<br>"), "<br>")
    } else{
      cat("<br>")
    }
    
  })
    
  
  ## Render the transcriptome results for NCBI GEO 
  output$tableTrans <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$transcriptome[2]$df },
     escape = FALSE
  )

  
  ## Render the header of the transcriptome ArrayExpress results
  output$headerTransAE <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    
    cat( "ArrayExpress data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })
  
  
  ## Render the information panel for ArrayExpress errors that might have occurred
  output$headerTransAEError <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    if( length( chem$transcriptomeAE[1]$msg) > 0 ){
      cat( "Please note that error(s) occurred during data set collection:<br>")
      cat( paste(chem$transcriptomeAE[1]$msg, collapse = "<br>"), "<br>")
    } else{
      cat("<br>")
    }
    
  })
  
  
  ## Render the transcriptome results for ArrayExpress
  output$tableTransAE <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$transcriptomeAE[2]$df },
    escape = FALSE
  )
  
  
  ## Render the header for the proteome results 
  output$headerProt <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    
    cat( "Proteomic data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })

  
  ## Render the information panel for proteome errors that might have occurred
  output$headerProtError <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    if( length( chem$proteome[1]$msg) > 0 ){
      cat( "Please note that error(s) occurred during data set collection:<br>")
      cat( paste(chem$proteome[1]$msg, collapse = "<br>"), "<br>")
    } else{
      cat("<br>")
    }
    
  })
  
  
  ## Render the results table for PRIDE proteome queries
  output$tableProt <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$proteome[2]$df },
    escape = FALSE
  )  
  
  
  ## Render the header for the metabolome results tab
  output$headerMeta <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    
    cat( "MetabolomeXchange data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })
  
  
  ## Render the information panel for metabolome errors that might have occrued
  output$headerMetaError <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    if( length( chem$metabolome[1]$msg) > 0 ){
      cat( "Please note that error(s) occurred during data set collection:<br>")
      cat( paste(chem$metabolome[1]$msg, collapse = "<br>"), "<br>")
    } else{
      cat("<br>")
    }
    
  })
  
  ## Render the results table for metabolomeXchange queries
  output$tableMeta <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$metabolome[2]$df },
    escape = FALSE
  )

}



#' This function removes specified tabs from the current session
#' @param session The current session
#' @param displayedTabs A list of tabs that should be removed
removeTabs <- function( session, displayedTabs){
  
  ## remove compound information tab
  removeTab( inputId = "tabset", target = "Compound", session = session)
  
  ## remove tabs displaying search results in case 'search' was hit again
  ## otherwise, the old tabs will also be displayed
  if( "CTDbase" %in% displayedTabs){
    removeTab( inputId = "tabset", target = "CTDbase", session = session)
  }
  if( "trans" %in% displayedTabs){
    removeTab( inputId = "tabset", target = "Transcriptome", session = session)
  }
  if( "prot" %in% displayedTabs){
    removeTab( inputId = "tabset", target = "Proteome", session = session)
  }
  if( "meta" %in% displayedTabs){
    removeTab( inputId = "tabset", target = "Metabolome", session = session)
  }
  
  return( NULL)
  
}


# Run the application 
shinyApp(ui = ui, server = server)

