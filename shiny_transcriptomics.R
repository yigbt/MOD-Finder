#' The function takes a list of chemicals of interest to initiate a NCBI GEO search with each individual entry.
#' @param CoI_list List of chemicals of interest.
#' @return A list containing an errorstring and the result dataframe.
get_transcriptome_by_list <- function( CoI_list){
  
  
  nested_list <- lapply( CoI_list, function(x) get_transcriptome_by_name(x))
  
  ## create a list of all errors or warnings that occured during the query
  msgs <- unlist( lapply( nested_list, '[[', 1))
  msg <- msgs[ which( msgs != "")]

  ## create a combined dataframe of all CoIs in the list
  res <- do.call( rbind, lapply( nested_list, '[[', 2))

  # return the list
  return( list( msg = msg, df = unique( res)))
  
}

#' This function creates an SQL statement with a given chemical of interest (CoI) to retrieve all current data sets that are associated with this CoI from NCBI GEO.
#' @param CoI The chemical of interest.
#' @return A list containing an errorstring and the result dataframe.
get_transcriptome_by_name <- function( CoI){
  
  
  ## Establish the connection to the SQL database.
  con <- tryCatch(
    
    dbConnect( SQLite(), "data/GEOmetadb.sqlite"),
    
    warning = function( cond){
      message( "WARNING - SQL conncection could not be established to query NCBI GEO.")
      return( list( msg = paste0( "WARNING - SQL conncection could not be established to query NCBI GEO. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    },
    error = function(cond){
      message( "ERROR - SQL conncection could not be established to query NCBI GEO.")
      return( list( msg = paste0( "ERROR - SQL conncection could not be established to query NCBI GEO. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    }
    
  )
  

  ## query NCBI GEO SQL database
  ## get all datasets with a title that contains the compound
  gse <- tryCatch(
    
    dbGetQuery( con, paste0("SELECT gse.title, gse.gse, gse.pubmed_id, gse.type, COUNT ( gse_gsm.gsm) as '#samples', COUNT( DISTINCT gse_gpl.gpl) as '#platforms', gpl.organism FROM gse JOIN gse_gsm on gse.gse=gse_gsm.gse JOIN gse_gpl on gse.gse=gse_gpl.gse JOIN gpl on gse_gpl.gpl=gpl.gpl WHERE gse.title LIKE '%", CoI, "%' GROUP BY gse.gse")),
    
    warning = function( cond){
      message( "WARNING - NCBI SQL query resulted in a warning.")
      return( list( msg = paste0( "WARNING - NCBI SQL query resulted in a warning. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    },
    error = function(cond){
      message( "ERROR - NCBI SQL query resulted in an error.")
      return( list( msg = paste0( "ERROR - NCBI SQL query resulted in an error. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    }
    
  )
  
  ## include some refinements in the output table to shorten things (strings) up
  gse$type <- gsub( "high throughput sequencing", "HTS", gse$type)
  gse$gse <- create_geo_link( gse$gse)
  notNA <- which( !is.na( gse$pubmed_id))
  gse$pubmed_id[ notNA] <- create_pubmed_link( gse$pubmed_id[ notNA])
  
  colnames( gse) <- c( "Title", "GSE", "Pubmed ID", "Experiment", "#Samples", "#Platforms", "Organism")
  
  gse$Experiment <- gsub( "RNA-seq of coding RNA", "RNA-seq", gse$Experiment)
  gse$Experiment <- gsub( "Expression profiling by array", "Microarray", gse$Experiment)
  
  return( list( msg = "", df = gse))
  
}


#' The function takes a list of chemicals of interest to initiate an ArrayExpress search with each individual entry.
#' @param CoI_list List of chemicals of interest.
#' @return A list containing an errorstring and the result dataframe.
get_arrayExpress_by_list <- function( CoI_list){
  
  nested_list <- lapply( CoI_list, function( x) get_arrayExpress_by_name( x))
  
  ## create a list of all errors or warnings that occured during the query
  msgs <- unlist( lapply( nested_list, '[[', 1))
  msg <- msgs[ which( msgs != "")]
  
  ## create a combined dataframe of all CoIs in the list
  res <- do.call( rbind, lapply( nested_list, '[[', 2))
  
  # return the list
  return( list( msg = msg, df = unique( res)))
  
}

#' The function queries the RESTful API of ArrayExpress with a given chemical of interest (CoI) to retrieve all data sets that are associated with this CoI.
#' @param CoI The chemical of interest.
#' @return A list containing an errorstring and the result dataframe.
get_arrayExpress_by_name <- function( CoI){
  
  ## create RESTful request and try to run it
  request <- sprintf( "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords=\"%s\"", CoI)
  getobj <- tryCatch( 

    GET( request),

    warning = function( cond){
      message( "WARNING - ArrayExpress RESTful request resulted in a warning.")
      return( list( msg = paste0( "WARNING - ArrayExpress RESTful request resulted in a warning. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    },
    error = function(cond){
      message( "ERROR - ArrayExpress RESTful request resulted in an error.")
      return( list( msg = paste0( "ERROR - ArrayExpress RESTful request resulted in an error. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    }
    
  )  
      
  ## check if the returned object is error free
  if( http_error( getobj)){
    
    return( list( msg = paste0( "ERROR - ArrayExpress query resulted in an Internal Server Error! (Chemical: ", CoI, ")"),
                  df = data.frame()))
    
  }
  
  ## extract the data object, get the content
  content <- httr::content( getobj)

  ## in case no data sets could be found, return an empty dataframe  
  if( length( content) == 0){
    return( list( msg = "", df = data.frame()))
  }
  
  # do some reformating 
  ids <- sapply( content$experiments$experiment, function(x) x$id)
  acc <- sapply( content$experiments$experiment, function(x) sprintf( '<a href=\"https://www.ebi.ac.uk/arrayexpress/experiments/%s\" target=\"_blank\">%s</a>', x$accession, x$accession))
  title <- sapply( content$experiments$experiment, function(x) x$name)
  organisms <- sapply( content$experiments$experiment, function(x) paste( unlist( x$organism), collapse = ", "))
  exp <- sapply( content$experiments$experiment, function(x) unlist( x$experimenttype))
  exp <- gsub( "RNA-seq of coding RNA", "RNA-seq", exp)
  exp <- gsub( "transcription profiling by array", "Microarray", exp)
  
  ## return the list that contains an empty error message and the dataframe with information about data sets 
  return( list( msg = "", df = data.frame( "Title" = title, "IDs" = ids, "Accession" = acc, "Experiment" = exp, "Organism" = organisms)))
  
}

