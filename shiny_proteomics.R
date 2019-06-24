#' This function takes a list of chemicals of interest to query the EBI PRIDE database.
#' @param compound_list The list with names of chemicals of interest.
#' @return A list containing an errorstring and the result dataframe.
get_proteome_by_list <- function( CoI_list){
  
  nested_list <- lapply( CoI_list, function(x) get_proteome_by_name(x))
  
  ## create a list of all errors or warnings that occured during the query
  msgs <- unlist( lapply( nested_list, '[[', 1))
  msg <- msgs[ which( msgs != "")]
  
  ## create a combined dataframe of all CoIs in the list
  res <- do.call( rbind, lapply( nested_list, '[[', 2))
  
  ## return the list
  return( list( msg = msg, df = unique( res)))
  
}

#' This function queries the EBI PRIDE proteomics database with a specific chemical of interest.
#' Therefore, the RESTful API interface is utilized.
#' The retrieved JSON object is parsed to return a well-formated dataframe containing all necessary information about the Metabolomics data set.
#' @param compound The chemical of interest to query EBI PRIDE.
#' @return A list containing an errorstring and the result dataframe.
get_proteome_by_name <- function( CoI){

  request <- sprintf( "https://www.ebi.ac.uk/pride/ws/archive/project/list?query=\"%s\"", CoI)

  ## retrieve the RESTful object from EBI PRIDE
  getobj <- tryCatch(
    
    GET( request),
    
    warning = function(cond){
      message( "Querying PRIDE resulted in a warning.")
      return( list( msg = paste0( "WARNING - EBI PRIDE RESTful request resulted in a warning. (Chemical: ", CoI, ")"),
                    df = data.frame()))
    },
    error = function(cond){
      message( "Querying PRIDE resulted in an error.")
      return( list( msg = paste0( "ERROR - EBI PRIDE RESTful request resulted in an ERROR. (Chemical: ", CoI,")"),
                    df = data.frame()))
    }
  )
  
  ## check if there are any Server problems
  if( http_error( getobj)){
    
    return( list( msg = paste0("ERROR - EBI PRIDE query resulted in an Internal Server Error. (Chemical: ", CoI, ")"), 
                  df = data.frame())
    )
    
  }
  
  ## extract the data object, get the content
  content <- httr::content( getobj)
  
  ## return an empty dataframe and error message in case no data set can be founnd with this CoI
  if( length( content[1]$list) == 0){
    return( list( msg = "", df = data.frame()))
  }
  
  ## do some reformating 
  ids <- sapply( content$list, function(x) x$accession)
  titles <- sapply( content$list, function(x) x$title)
  subType <- sapply( content$list, function(x) x$submissionType)
  numAssay <- sapply( content$list, function(x) x$numAssays)
  species <- unlist( sapply( content$list, function(x) paste( x$species, collapse = ", ")))
  species <- gsub( " \\(.*\\)", "", species)
  
  ## generate the functional URLs
  ids <- sprintf( '<a href=\"https://www.ebi.ac.uk/pride/archive/projects/%s\" target=\"_blank\">%s</a>', ids, ids)

  ## return the list that contains an empty error message and the dataframe with information about data sets 
  return( list( msg = "", df = data.frame( "Title" = titles, "IDs" = ids, "Submission Type" = subType, "Nr of Assays" = numAssay, "Organism" = species)))
  
}
