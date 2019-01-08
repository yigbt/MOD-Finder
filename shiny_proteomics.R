
## query EBI PRIDE by a list of compounds
get_proteome_by_list <- function( compound_list){
  
  res <- do.call( rbind, lapply( compound_list, function(x) get_proteome_by_name(x)))
  return( unique( res))
  
}


## use the API REST property of EBI PRIDE database to search for proteomics studies
## where the title and/or the description contain a certain compound name
get_proteome_by_name <- function( compound){
  
  request <- sprintf( "https://www.ebi.ac.uk/pride/ws/archive/project/list?query=\"%s\"", compound)
  content <- httr::content( GET( request))
  
  if( length( content) == 0) return( data.frame())
  
  ids <- sapply( content$list, function(x) x$accession)
  titles <- sapply( content$list, function(x) x$title)
  subType <- sapply( content$list, function(x) x$submissionType)
  numAssay <- sapply( content$list, function(x) x$numAssays)
  species <- unlist( sapply( content$list, function(x) paste( x$species, collapse = ", ")))
  species <- gsub( " \\(.*\\)", "", species)
  
  ids <- sprintf( '<a href=\"https://www.ebi.ac.uk/pride/archive/projects/%s\" target=\"_blank\">%s</a>', ids, ids)
  
  return( data.frame( "Title" = titles, "IDs" = ids, "Submission Type" = subType, "Nr of Assays" = numAssay, "Organism" = species))
  
}
