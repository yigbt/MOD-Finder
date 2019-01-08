
## query metabolomXchange by a list of compounds
get_metabolomeXchange_by_list <- function( compound_list){

  res <- do.call( rbind, lapply( compound_list, function(x) get_metabolomeXchange_by_name( x)))
  return( unique( res))
  
}


## query metabolomXchange by a single compound
get_metabolomeXchange_by_name <- function( compound){
  
  
  ## how about compounds including a space ??
  compound <- gsub( " ", "%20", compound)
  request <- sprintf( "http://api.metabolomexchange.org/datasets/%s", compound)
  content <- httr::content( GET( request))
  
  if( length( content) == 0) return( data.frame())
  
  url <- sapply( content, function(x) x$url)
  title <- sapply( content, function(x) x$title)
  acc <- sapply( content, function(x) x$accession)
  provider <- sapply( content, function(x) x$provider)
  
  org <- sapply( content, function(x) paste( unlist(x$meta$organism), collapse = ", "))
  analysis <- sapply( content, function(x) paste( unlist(x$meta$analysis), collapse = ", "))
  platform <- sapply( content, function(x) paste( unlist(x$meta$platform), collapse = ", "))
  org <- sapply( content, function(x) paste( unlist(x$meta$organism), collapse = ", "))
  desc <- sapply( content, function(x) paste( unlist(x$description), collapse = "; "))
  
  ## there are at least 235 different entities for the 'analysis' term
  ## therefore i will skip the mapping of terms such as 'mass spectrometry' to 'MS'
  #  values <- c( "NMR spectroscopy", "mass spectrometry", "MS analysis", "NMR development for metabolomics", "Untargeted HRMS for cross-laboratory characterization of common reference material.", "Comparison of degradation kinetics", "Endpoint measurement", "GC-MS non-targeted metabolomic profiling", "HERITAGE (HEalth, RIsk factors, exercise Training And GEnetics) family study", "Regular", "Tomato Seed Metabolites Profiling (dry seed and 6 hour imbibed seeds comparision)", "Lipid analysis novel C18 fatty acid anologues in complex lipids", "Acyl-carnitine analysis (plasma)", "Cell labeling",  "High and low insulin with and without essential amino acids", "in vitro study/drug dosage", "Genotype", "Genotype treatment")
  #  replacement <- c( "NMR", "MS", "MS", "NMR", "HRMS", "Comparison of degradation kinetics", "Endpoint measurement", "GC-MS", "HERITAGE", "Regular", "", "Lipid analysis", "", "Cell labeling",  "", "", "", "")
  #  map_analysis <- plyr::mapvalues( analysis, from = values, to = replacement)  
  
  ## some rearrangements for convenience
  map_provider <- plyr::mapvalues( provider, from = c("mtbls", "mnote", "mwbs", "meryb"), to = c("MetaboLights", "Metabolonote", "Metabolomics Workbench", "MeRy-B"))
  link <- sprintf( '<a href=\"%s">%s</a>', url, acc)
  
  ## create doi list and pubmed list
  publication <- sapply( content, function(x) paste( unlist( x$publication), collapse = ", "))
  pmid <- gsub( "PMID:", "", stringr::str_extract_all( pattern = "PMID:\\d+", string = publication, simplify = TRUE))
  pmid <- apply( pmid, 1, function( row) {
    gsub( ",+$", "", 
          paste( create_pubmed_link( row), 
                 collapse = ","))
  })
  
  ## which descriptions and title do actually contain the given compound?
  rows <- unique( c( grep( pattern = compound, desc), grep( pattern = compound, title)))
  
  return( data.frame( "Title" = title, "Accession" = link, "PubmedID" = pmid, "Organism" = org, "Analysis" = analysis, "Source" = map_provider )[ rows, ])
  
}




#### DOWN BELOW ARE OLD AND CURRENTLY UNUSED FUNCTION
#### THEY WERE WRITTEN TO QUERY METABOLOMICS WORKBENCH AND METABOLIGHTS SEPARATELY
#### AND LATER ON THEY WERE REPLACED BY A FUNCTION THAT QUERIES METABOLOMEXCHANGE
#### THERE IS ALSO A FUNCTION TO DOWNLOAD AND PARSE THE METABLIGHTS XML FILE TO QUERY THIS SOURCE


## query metabolomics workbench by a list of compounds
get_metabolome_by_list <- function( compound_list){
  
  res <- do.call( rbind, lapply( compound_list, function( x) get_metabolome_by_name( x)))
  return ( unique( res))
  
}


## query metabolomics workbench by a single compound
get_metabolome_by_name <- function( compound){
  
  ## this function uses the REST api from metabolomics workbench to search for potential metabolomics data sets given a specific compound.
  request <- sprintf( "http://www.metabolomicsworkbench.org/rest/study/study_title/%s/summary/", gsub( " ", "%20", compound))
  content <- httr::content( GET( request))
  
  if( length( content) == 0) return( data.frame())
  
  
  ids <- sapply( content, function(x) x$study_id)
  titles <- sapply( content, function(x) x$study_title)
  types <- sapply( content, function(x) ifelse( is.null(x$study_type), "", x$study_type))
  organisms <- sapply( content, function(x) x$subject_species)
  
  ids <- sprintf( '<a href=\"http://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=%s\" target=\"_blank\">%s</a>', ids, ids)
  
  return( data.frame( "Title" = titles, "IDs" = ids, "Experiment Type" = types, "Organism" = organisms))
  
} 


# query metabolights by a list of compounds
get_metaboLights_by_list <- function( compound_list){
  
  res <- do.call( rbind, lapply( compound_list, function(x) get_metaboLights_by_name( x)))
  return( unique( res))
  
}


## query metabolights by a single compound
get_metaboLights_by_name <- function( compound){
  
  
  if( !file.exists( "eb-eye_metabolights_studies.RData")){
    parse_metaboLights_xml_file()
  }  
  load( "eb-eye_metabolights_studies.RData")
  
  
  ## search for the compound in Title, Decription, and Publication Info
  hits <- grep( compound, metaboLights$Title)
  hits <- c( hits, grep( compound, metaboLights$Description))
  hits <- c( hits, grep( compound, metaboLights$Publication.Info ))
  
  hits <- sort(unique( hits))
  
  metaboLights$Link <- sprintf( '<a href=\"%s\" target=\"_blank\">%s</a>', metaboLights$Link, metaboLights$ID)
  metaboLights$Pubmed.ID <- create_pubmed_link( metaboLights$Pubmed.ID)
  
  
  return( metaboLights[ hits, c( 2, 1, 4, 5, 12, 6, 10)])
  
}


## download and parse the whole METABOLIGHTS xml file containing metadata to all provided studies
parse_metaboLights_xml_file <- function(){
  
  
  inputxml <- "eb-eye_metabolights_studies.xml"
  outfile <- "eb-eye_metabolights_studies.RData"
  
  ## download xml if not present
  if( !file.exists( inputxml)){
    download.file( "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/eb-eye/eb-eye_metabolights_studies.xml", inputxml)
  }
  
  ## read the given xml file
  xml <- xml2::read_xml( inputxml)
  
  ## create a nodeset pointing of the 'entry' nodes in the xml file
  entries <- xml2::xml_find_all( xml, ".//entry")
  
  ## save the entry ids into a list of ids
  ids <- xml_attr( entries, "id")
  studies_list <- grep( "MTBLS", ids)
  
  ids <- ids[ studies_list]
  
  ## get lists of children nodes of each study related entry
  names <- xml_find_all( entries, ".//name")[ studies_list]
  desc <- xml_find_all( entries, ".//description")[ studies_list]
  cross <- xml_find_all( entries, ".//cross_references")[ studies_list]
  dates <- xml_find_all( entries, ".//dates")[ studies_list]
  addons <- xml_find_all( entries, ".//additional_fields")[ studies_list]
  
  ## these lists should have the same length!!!
  if( length( unique( c(length(names), length( desc), length( cross), length( dates), length( addons)))) != 1 ){
    cat( "ERROR: something's wrong with the xml file")
  }
  
  ## get all repos
  repo <- as.character( xml_contents( xml_find_all( addons, ".//field[@name='repository']")))
  
  ## get all omics types
  omics <- as.character( xml_contents( xml_find_all( addons, ".//field[@name='omics_type']")))
  
  ## get all publication status'
  status <- as.character( xml_contents( xml_find_all( addons, ".//field[@name='study_status']")))
  
  ## get comma-separated list of all organisms
  org <- sapply( addons, function(x) paste( unique( as.character( xml_contents( xml_find_all( x, ".//field[@name='Organism']")))), collapse = ","))
  org <- gsub( "NCBITAXON:|BAO:|reference compound,?", "", org)
  
  ## get comma-separated list of all organism body parts
  org_part <- sapply( addons, function(x) paste( unique( as.character( xml_contents( xml_find_all( x, ".//field[@name='Organism Part']")))), collapse = ","))
  
  ## get available publication-relevant metadata
  publication <- sapply( addons, function(x) paste( unique( as.character( xml_contents( xml_find_all( x, ".//field[@name='publication']")))), collapse = ","))
  
  ## extract the Pubmed IDs if present
  pmid <- ifelse( grepl("PMID:", publication), gsub( ".*PMID:(\\d+).*", "\\1", publication), "")
  
  ## get the DOI if present
  doi <- ifelse( grepl( "\\d+/[a-z0-9\\-\\.]+[\\s\\.]", publication), gsub( ".*\\s([0-9\\.]+/[a-z0-9\\.\\-]+)[\\.\\s].*", "\\1", publication), "")
  doi <- as.character( sapply( doi, 
                               function(x) {
                                 http <- resolve_doi( x)
                                 ifelse( http != "", sprintf( '<a href=\"%s\" target=\"_blank\">%s</a>', http, x), "")},
                               simplify = TRUE))
  
  ## link to the dataset
  link <- as.character( xml_contents( xml_find_all( addons, ".//field[@name='full_dataset_link']")))
  
  ## gather information about the experiment
  technology <- sapply( addons, function(x) paste( unique( as.character( xml_contents( xml_find_all( x, ".//field[@name='technology_type']")))), collapse = ","))
  technology <- gsub( "mass spectrometry", "MS", technology)
  technology <- gsub( "NMR spectroscopy", "NMR", technology)
  
  
  metaboLights <- data.frame( "ID" = ids, "Title" = as.character( xml_contents(names)), "Description" = as.character( xml_contents(desc)), "Pubmed ID" = pmid, "DOI" = doi, "Organism" = org, "Organism Part" = org_part, "Repository" = repo, "Omics Type" = omics, "Status" = status, "Publication Info" = publication, "Technology" = technology, "Link" = link)
  
  save( metaboLights, file = outfile)
  
}

