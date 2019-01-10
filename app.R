#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library( shiny)
library( tidyverse)
library( CTDquerier)
library( httr)
library( GEOmetadb)
library( DT)
library( xml2)
library( plyr)
library( stringr)

## source some code
source( "shiny_metabolomics.R", local = TRUE)
source( "shiny_proteomics.R", local = TRUE)
source( "shiny_transcriptomics.R", local = TRUE)
source( "shiny_general.R", local = TRUE)
source( "shiny_compound.R", local = TRUE)



######
# CONSTANTS
pubchem_compound <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
pubmed_url <- "https://www.ncbi.nlm.nih.gov/pubmed/"
geo_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="

######
### try catch for server error warnings!!!!!
######



# to query proteomexchange database
#library( rpx)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   #titlePanel("Search for multi-omics data sets in various public sources."),
   titlePanel(
     fluidRow( column( 6, h1( "MOD-Finder"), HTML( "<h3>Search for <b>M</b>ulti-<b>O</b>mics <b>D</b>ata sets in various public sources.</h3>")),
               column( 6, img( height = 120, src="UFZ_logo.jpg", class="pull-right"))
               ),
     windowTitle = "MOD-Finder"
   ),
   
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     
      sidebarPanel(
        
        width = 3,
        
        h4( "Step I - Enter Compound"),
        
        textInput( inputId = "compoundName", label = "Compound Name", value = "Enter Compound", width = '400px'),

        checkboxGroupInput( inputId = "omicsLayer", choices = c( "Transcriptome" = "trans", "Proteome" = "prot", "Metabolome" = "meta"), label = "Omics Layer" ),
        
        actionButton( inputId = "refineCompound", label = "Refine"),
        
        hr(),
        
        h4( "Step II - Filter the exact Compound"),
        
        selectInput( inputId = "compoundList", choices = NULL,
                     label = "Specify the compound to search for:"),
        
        radioButtons( inputId = "splitCompound", 
                    label = "Name splitting", choices = c( "Search with refined compound")),
        
        actionButton( inputId = "startSearch", label = "Search")
                
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        
        
        tabsetPanel( id = "tabset"
        
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

    withProgress( message = "Search similar compound names", value = 0.5, {
      compound_list <- get_cids_by_name( chem$compound)
      if( nrow( compound_list ) >= 1 ){
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
      
      
      ## initialize a progress bar showing information about the current step
      withProgress( message="Collecting compound information", value=0, {
      
        if( !is.null( chem$cid)){ 
          
        ## collect information from pubchem using PUG REST
        ## first: CID and synonyms
        incProgress( amount = 0.2, detail = "Collecting Pubchem Information")
        pubchem <- get_synonyms_by_cid( chem$cid)
        synonyms <- pubchem$synonyms


        ## second: general information, such as SMILE and descriptions
        chem$results <- get_compound_information( chem$cid)

        ## get pubchem ID mappings
        pubchem <- get_description_by_cid( chem$cid)
              
        incProgress( amount = 0.2, detail = "Collecting CTD Information")

        ## get information about the compound from the CTDbase
        ## in order to guarantee that the correct compound is found,
        ## use the MeSH ID retrieved from Pubchem to search the CTDbase
        ## in case this is not successful,
        ## use the query_CTD_chem approach which searches for compounds with a specified edit-distance
        chem$ctd <- NULL
        if( "MeSH" %in% pubchem$sources){
          
          mesh <- strsplit( pubchem$functional_link[ grep( "MeSH", pubchem$sources)], split = ", ")[[1]]
          if( length( mesh) > 1) chem$ctd <- query_ctd_by_id( mesh[2])

        } 
        
        if( is.null( chem$ctd)){

          res <- try( { 

            ctd <- query_ctd_chem( terms = chem$exactCompound, max.distance = 1)
            chem$ctd$Chemical.Name <- ctd@terms@listData$ChemicalName
            chem$ctd$Chemical.ID <- ctd@terms@listData$ChemicalID
            if( ctd@terms@listData$CasRN != "") chem$ctd$CasRN <- ctd@terms@listData$CasRN
            chem$ctd$Synonyms <- ctd@terms@listData$Synonyms
            if( ctd@terms@listData$DrugBankIDs != "") chem$ctd$DrugBankIDs <- ctd@terms@listData$DrugBankIDs

          }, silent = TRUE)
          
        }


        ## create a data frame which stores all IDs collected during the runtime
        chem$ids <- data.frame( Database = c( "Pubchem"), IDs = c( sprintf( '<a href=\"https://pubchem.ncbi.nlm.nih.gov/compound/%s\" target=\"_blank\">%s</a>', chem$cid, chem$cid)), stringsAsFactors = FALSE)
        
        if( "DrugBank" %in% pubchem$sources & !is.null( chem$ctd$CasRN)) {
          chem$ids <- rbind( chem$ids, data.frame( Database = c( "CAS ID"), IDs = c( chem$ctd$CasRN), stringsAsFactors = FALSE))
        } else if( !is.null( chem$ctd$CasRN) & !is.null( chem$ctd$DrugBankIDs)) {
          chem$ids <- rbind( chem$ids, data.frame( Database = c( "CAS ID", "DrugBank"), IDs = c( chem$ctd$CasRN, chem$ctd$DrugBankIDs), stringsAsFactors = FALSE))
        } else if( !is.null( chem$ctd$CasRN)){
          chem$ids <- rbind( chem$ids, data.frame( Database = c( "CAS ID"), IDs = c( chem$ctd$CasRN), stringsAsFactors = FALSE))
        } else if( !is.null( chem$ctd$DrugBankIDs)){
          chem$ids <- rbind( chem$ids, data.frame( Database = c( "DrugBank"), IDs = c( chem$ctd$DrugBankIDs), stringsAsFactors = FALSE))
        }
        chemical_name <- chem$ctd$Chemical.Name
      
      
        ## if the synonyms list has more than 10 entries, show only the first ten
        ## keep the complete list in chem$synonyms
        ## and save the output list as a single string in chem$synonyms_list
        if( !is.null( chem$ctd$Synonyms)){
          chem$synonyms <- strsplit( as.character( chem$ctd$Synonyms), split = "\\|")[[1]]
        }
        chem$synonyms <- unique( tolower( c( chem$synonyms, synonyms)))

        if( length( chem$synonyms > 10)){
          chem$synonyms_list <- paste0( paste( chem$synonyms[1:10], collapse = " | "), " ... (10 of ", length( chem$synonyms), ")" )
        } else{
          chem$synonyms_list <- paste( chem$synonyms, collapse = " | ")
        }
      
        ## combine the information into two tables
        ## one focusing on the chemical properties and descriptions
        ## and the other focusing on the IDs conversions and mappings
        chem$results <- rbind( "Chemical Name" = chemical_name, chem$results)
        chem$results <- rbind( chem$results, "Synonyms" = chem$synonyms_list)
        chem$ids <- rbind( chem$ids, data.frame( Database = pubchem$sources, IDs = pubchem$functional_link))
          
        } else{
          
          chem$error <- "No compound information found for initial input name!"
          
        }
      
        ## add the compound information panel
        appendTab( inputId = "tabset",
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
        
        
        
        ## check which omics layer ist checked and add the respective tab in the panel
        ## search for data sets only in case the respective omics layer is selected
        progress_time = 0.6 / length( input$omicsLayer)
        if( "trans" %in% input$omicsLayer ){

          ## collect potential transcriptome data sets
          ## use GeoMetadb to search for uploaded data sets were the title contains the compound name
          incProgress( amount = progress_time, detail = "Collecting Transcriptome Data Sets")
          chem$transcriptome <- get_transcriptome_by_list( chem$exactCompoundList)
          chem$transcriptomeAE <- get_arrayExpress_by_list( chem$exactCompoundList)
          
          ## create a dynamic tab to display the previously gathered information about transcriptomics

          lines <- nrow( chem$transcriptome) + nrow( chem$transcriptomeAE)
          appendTab( inputId = "tabset",
                     session = session,
                     tabPanel( value = "Transcriptome",
                               title = paste0( "Transcriptome (", lines, ")"),
                               h3( textOutput( outputId = "headerTrans")),
                               DT::dataTableOutput( outputId = "tableTrans"),
                               br(),
                               h3( textOutput( outputId = "headerTransAE")),
                               DT::dataTableOutput( outputId = "tableTransAE"))
          )
          
        }
      
        if( "prot" %in% input$omicsLayer){
        
          ## collect potential proteome data sets
          ## use the rpx package to search amongst proteinXchange-uploaded data sets
          incProgress( amount = progress_time, detail = "Collecting Proteome Data Sets")
          chem$proteome <- get_proteome_by_list( chem$exactCompoundList)
        
          ## create a dynamic tab to display the previously gathered information about proteomics
          appendTab( inputId = "tabset",
                     session = session,
                     tabPanel( value = "Proteome",
                               title = paste0( "Proteome (", nrow( chem$proteome), ")"),
                               h3(textOutput( outputId = "headerProt")),
                               DT::dataTableOutput( outputId = "tableProt"))
          )
          
        }
        
        if( "meta" %in% input$omicsLayer){
        
          ## collect potential metabolome data sets
          ## use ?? to search Metabolights and Metabolomics Workbench
          incProgress( amount = progress_time, detail = "Collecting Metabolome Data Sets")
          chem$metabolome <- get_metabolomeXchange_by_list( chem$exactCompoundList)
          
          hits <- nrow( chem$metabolome)
          ## create a dynamic tab to display the previously gathered information about metabolomics
          appendTab( inputId = "tabset",
                     session = session,
                     tabPanel( value = "Metabolome",
                               title = paste0( "Metabolome (", hits, ")"),
                               h3(textOutput( outputId = "headerMeta")),
                               DT::dataTableOutput( outputId = "tableMeta"))
          )
        }

      })

    }
    
    chem$displayedTabs <- input$omicsLayer
    
  })
  
  output$headerCompound <- renderPrint( {
    
    if( chem$searchCompound == FALSE) return()

    if( "error" %in% names(chem) & !is.null( chem$error)){
      cat( "An ERROR occurred: \n", chem$error)
    } else if( "exactCompound" %in% names(chem) & !is.null( chem$exactCompound)){
      cat( "Compound Information: ", chem$exactCompound)
    }else{
      cat( "Something really strange happened here!")
    }

  })
  
  output$tableCompound <- renderTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    if( !is.null( chem$results)){
      head( chem$results, n = nrow( chem$results))}},
  
      rownames = TRUE,
      colnames = FALSE
  )
  
  output$headerCompoundIDs <- renderPrint( {
    
    if( chem$searchCompound == FALSE) return()
    if( !is.null( chem$ids)){
      cat( "ID conversion and mapping")
    }
    
  })
  
  output$tableCompoundIDs <- DT::renderDataTable({
    
   if( chem$searchCompound == FALSE ) return()
    
    if( !is.null( chem$ids)){
      head( chem$ids, n = nrow( chem$ids))}
    },

    rownames = FALSE,
    escape = FALSE
  
  )
  
  
  output$headerTrans <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()

    cat( "NCBI Geo data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })
  
  output$tableTrans <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$transcriptome },
     escape = FALSE
  )

  
  output$headerTransAE <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    
    cat( "ArrayExpress data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })
  
  
  output$tableTransAE <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$transcriptomeAE },
    escape = FALSE
  )
  
  
  output$headerProt <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    
    cat( "Proteomic data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })
  
  output$tableProt <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$proteome },
    escape = FALSE
  )  
  
  output$headerMeta <- renderPrint( {
    
    if( chem$searchCompound == FALSE ) return()
    
    cat( "MetabolomeXchange data sets with ", paste( chem$exactCompoundList, collapse = ", "))
    
  })
  
  output$tableMeta <- DT::renderDataTable({
    
    if( chem$searchCompound == FALSE ) return()
    
    chem$metabolome },
    escape = FALSE
  )

}


removeTabs <- function( session, displayedTabs){
  
  ## remove compound information tab
  removeTab( inputId = "tabset", target = "Compound", session = session)
  
  ## remove tabs displaying search results in case 'search' was hit again
  ## otherwise, the old tabs will also be displayed
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

