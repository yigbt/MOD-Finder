
## query NCBI GEO by a list of compounds
get_transcriptome_by_list <- function( compound_list){
  
  res <- do.call( rbind, lapply( compound_list, function(x) get_transcriptome_by_name(x)))
  return( unique( res))
  
}


## query NCBI GEO by a single compound
get_transcriptome_by_name <- function( compound){
  
  
  
  con <- dbConnect( SQLite(), "data/GEOmetadb.sqlite")
  
  ## get all datasets with a title that contains the compound
  gse <- dbGetQuery( con, paste0("SELECT gse.title, gse.gse, gse.pubmed_id, gse.type, COUNT ( gse_gsm.gsm) as '#samples', COUNT( DISTINCT gse_gpl.gpl) as '#platforms', gpl.organism FROM gse JOIN gse_gsm on gse.gse=gse_gsm.gse JOIN gse_gpl on gse.gse=gse_gpl.gse JOIN gpl on gse_gpl.gpl=gpl.gpl WHERE gse.title LIKE '%", compound, "%' GROUP BY gse.gse"))
  
  ## include some refinements in the output table to shorten things (strings) up
  gse$type <- gsub( "high throughput sequencing", "HTS", gse$type)
  gse$gse <- create_geo_link( gse$gse)
  notNA <- which( !is.na( gse$pubmed_id))
  gse$pubmed_id[ notNA] <- create_pubmed_link( gse$pubmed_id[ notNA])
  
  colnames( gse) <- c( "Title", "GSE", "Pubmed ID", "Experiment", "#Samples", "#Platforms", "Organism")
  
  gse$Experiment <- gsub( "RNA-seq of coding RNA", "RNA-seq", gse$Experiment)
  gse$Experiment <- gsub( "Expression profiling by array", "Microarray", gse$Experiment)
  
  return( gse)
  
}


## query ArrayExpress by a list of compound
get_arrayExpress_by_list <- function( compound_list){
  
  res <- do.call( rbind, lapply( compound_list, function( x) get_arrayExpress_by_name( x)))
  return( unique( res))
  
}


## query ArrayExpress by a single compound
get_arrayExpress_by_name <- function( compound){
  
  request <- sprintf( "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords=\"%s\"", compound)
  content <- httr::content( GET( request))
  
  if( length( content) == 0) return( data.frame())
  
  ids <- sapply( content$experiments$experiment, function(x) x$id)
  acc <- sapply( content$experiments$experiment, function(x) sprintf( '<a href=\"https://www.ebi.ac.uk/arrayexpress/experiments/%s\" target=\"_blank\">%s</a>', x$accession, x$accession))
  title <- sapply( content$experiments$experiment, function(x) x$name)
  organisms <- sapply( content$experiments$experiment, function(x) paste( unlist( x$organism), collapse = ", "))
  exp <- sapply( content$experiments$experiment, function(x) unlist( x$experimenttype))
  exp <- gsub( "RNA-seq of coding RNA", "RNA-seq", exp)
  exp <- gsub( "transcription profiling by array", "Microarray", exp)
  
  return( data.frame( "Title" = title, "IDs" = ids, "Accession" = acc, "Experiment" = exp, "Organism" = organisms))
  
}

