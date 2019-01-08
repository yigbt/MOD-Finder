
query_ctd_by_id <- function( input_id){
  
  ## check if the RData file is present
  ## if not download and preprocess the file
  if( !file.exists( "data/CTD_chemicals.RData")){
    prepare_CTD_database_file()
  }
  
  ## load the CTD file
  load( file="data/CTD_chemicals.RData")
  
  ## search for the compound ID in the 2nd column
  chem <- ctd[ which( ctd$Chemical.ID == paste0( "MESH:", input_id)), ]
  
  
  if( nrow( chem) > 0){
    
    ## eliminate empty ID fields  
    if( chem$DrugBankIDs == "") chem$DrugBankIDs <-NULL
    if( chem$CasRN == "") chem$CasRN <- NULL
    
    return( chem)
    
  }
  
  return( NULL)
}


get_synonyms_by_cid <- function( input_cid){
  
  ## create a http request which is used to retrieve all synonyms of a given compound
  ## and the CID
  request <- paste0( pubchem_compound, input_cid, "/synonyms/JSON")
  
  content <- httr::content( GET( request))
  cid <- content$InformationList$Information[[1]]$CID
  
  synonyms <- unlist( content$InformationList$Information[[1]]$Synonym)
  
  return( list( cid=cid, synonyms=synonyms))
  
}


get_description_by_cid <- function( input_cid){
  
  request <- paste0( pubchem_compound, input_cid, "/description/JSON")
  
  content <- httr::content( GET( request))
  
  url <- unlist(sapply( content$InformationList$Information, function(x) x$DescriptionURL))
  sources <- unlist(sapply( content$InformationList$Information, function(x) x$DescriptionSourceName))
  
  ## remove FDA IDs and URLs since these do not have parsable unique identifier
  if( length( grep( "FDA", sources)) > 0){
    url <- url[-grep( "FDA", sources)]
    sources <- sources[ -grep( "FDA", sources)]
  }
  
  ## some kind of cleaning 
  sources <- gsub( "Human Metabolome Database \\(HMDB\\)", "HMDB", sources)
  sources <- gsub( "CAMEO Chemicals", "CAMEO", sources)
  
  ## get the corresponding IDs from the urls
  ids <- gsub( ".*\\/", "", url)
  ids <- gsub( ".*=", "", ids)
  
  ## turn the url into a functional link which can be clicked in the output table
  functional_link <- sprintf( '<a href=\"%s\" target=\"_blank\">%s</a>', url, ids)
  
  ## try to retrieve the other MeSH identifier
  ## current MeSH ID looks like this: 68001241
  ## but MeSH IDs used by CTDbase are more like this: D001241
  ## therefore, use the RESTful API from NCBI to get a mapping
  if( "MeSH" %in% sources ){
    
    meshid <- ids[ grep( "MeSH", sources)]
    
    ## create request object for the current MeSH ID
    ## and get the json object
    request <- paste0( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=mesh&id=",
                       meshid,"&retmode=json")
    content <- httr::content( GET( request))
    
    ## check if another identifier is in the retrieved json object
    if( "ds_meshui" %in% names(content$result[[ meshid]])){
      
      second_meshid <- content$result[[ meshid]]$ds_meshui
      ## add the second identifier to the list of IDs, but this time without a functional link
      functional_link[ grep( "MeSH", sources)] <- paste0( functional_link[ grep( "MeSH", sources)], ", ", second_meshid)
      
    }
    
  }
  
  return ( list( sources = sources, urls = url, ids = ids, functional_link = functional_link))
  
}


get_cids_by_name <- function( name){
  
  ## create a http request which is used to retrieve all synonyms of a given compound
  ## and the CID
  request <- paste0( "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                     name, "/synonyms/JSON?name_type=word")
  content <- httr::content( GET( request))
  
  if( "Fault" %in% names(content)){
    return( data.frame())
  }
  
  cid <- sapply( content$InformationList$Information, function(x) x$CID)
  pattern <- paste0( "^\\D*\\d\\D*", name, "\\D*$", "|" ,"^\\D*", name, "\\D*\\d\\D*$", "|", "^\\D*", name, "\\D*$")
  
  synonyms_list <- sapply( content$InformationList$Information, function(x) x$Synonym)
  synonyms <- sapply( synonyms_list, function(x){ x[grep( pattern, x, ignore.case = TRUE)[1]] })
  
  ## remove the NULL elements from synonyms and, beforehand, the respective IDs from the cid list
  cid <- cid[!sapply( synonyms, is.null)]
  synonyms[sapply( synonyms, is.null)] <- NULL
  
  results <- data.frame( cid = cid, name = unlist( synonyms) )
  
  return( arrange( results, cid))
  
}


get_compound_information <- function( cid){
  
  
  ## get canonical information from the properties table
  request <- paste0( pubchem_compound, cid, "/property/MolecularFormula,CanonicalSMILES,InChI/JSON")
  
  content <- httr::content( GET( request))
  
  properties <- content$PropertyTable$Properties[[1]]
  properties$CID <- NULL
  properties$InChI <- gsub( "InChI=", "", properties$InChI)
  
  ## get also the description
  request <- paste0( pubchem_compound, cid, "/description/JSON")
  
  content <- httr::content( GET( request))
  desc_list <- content$InformationList$Information
  
  ## get a list of potential description sources
  sources <- lapply( desc_list, function(x){ x$DescriptionSourceName})
  
  ## check if one of the following sources provides a compound description: MeSH, ChEBI, DrugBank, HMDB, etc...
  reference_sources <- c( "MeSH", "DrugBank", "Human Metabolome Database (HMDB)", "CAMEO Chemicals", "NCIt", "ChEBI")
  
  ## get a mapping between the reference sources and the current sources
  mapping <- BiocGenerics::match( reference_sources, sources)
  
  ## extract the actual information, i.e., the description and the corresponding source
  position <- mapping[ !is.na(mapping)][1]
  
  if( !is.na( position) && position >= 1 && position <= length( desc_list)){
    
    properties$Description <- desc_list[[ position]]$Description
    properties$DescriptionSource <- desc_list[[ position]]$DescriptionSourceName
    
    return( data.frame( values = unlist( properties), row.names = c( "Molecular Formular", "Chanonical SMILE", "InChI", "Definition", "Definition Source"), stringsAsFactors = FALSE))
    
  } else{
    
    properties$Description <- NULL
    properties$DescriptionSource <- NULL
    return( data.frame( values = unlist( properties), row.names = c( "Molecular Formular", "Chanonical SMILE", "InChI"), stringsAsFactors = FALSE))
    
  }
  
}


## Download the CTD_chemicals file from CTDbase.org
## parse the file and create an RData file containing a data.frame with an appropriate header
prepare_CTD_database_file <- function(){
  
  download_ctd_chem()
  if( !file.exists( "data/CTD_chemicals.tsv.gz")){
    
    cat( "ERROR: Could not download CTD_Chemicals file!\n")
    return( NULL)  
  }
  
  ctd <- read.csv( "data/CTD_chemicals.tsv.gz", header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE )
  colnames( ctd) <- c( "Chemical.Name", "Chemical.ID", "CasRN", "Definition", "ParentIDs", "TreeNumbers", "ParentTreeNumbers", "Synonyms", "DrugBankIDs")
  
  save( ctd, file="data/CTD_chemicals.RData")
  
}
