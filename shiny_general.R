
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


#' Check if the given CAS number is correct according to the CAS definition
#' Format fo a CAS number is: Z_n ... Z_4 Z_3 - Z_2 Z_1 - R
#' sum over product of position_i x Z_i  modulo 10 should be R
#' @param cas The CAS number.
#' @return Boolean
#'  
is_cas_correct <- function( cas){
  
  split <- unlist( strsplit( cas, split = "-"))
  remain <- as.integer( split[3])
  cas_numericals <- as.integer( c( unlist( strsplit( split[1], split = "")), unlist( strsplit( split[2], split = "")))) 

  sum_cas <- sum( cas_numericals * seq( length( cas_numericals), 1))
  
  if( sum_cas %% 10 == remain){
    return( TRUE)
  } else{
    return( FALSE)
  }
  
}
