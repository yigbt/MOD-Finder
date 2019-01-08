
create_geo_link <- function( geo_id){
  
  return( sprintf( '<a href=\"%s%s\" target=\"_blank\">%s</a>', geo_url, geo_id, geo_id) )
  
}


create_pubmed_link <- function( pubmed_id){
  
  
  
  return( ifelse( pubmed_id != "", sprintf( '<a href=\"%s%s\" target=\"_blank\">%s</a>', pubmed_url, pubmed_id, pubmed_id), ""))
  
}


resolve_doi <- function( doi){
  
  request <- paste0( "https://doi.org/api/handles/", doi)
  
  content <- httr::content( GET( request))
  
  if( content$responseCode == 1){
    return( content$values[[1]]$data$value)
  } else{
    return( "")
  }
  
}