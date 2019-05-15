
#' The function collects alle Pubchem-annotated synonyms for a given CID by means of a RESTful HTTP request.
#' @param cid The Pubchem CID of the compound.
#' @return A list containing the CIDs and the synonyms.
get_synonyms_by_cid <- function( cid){
  
  ## create a http request which is used to retrieve all synonyms of a given compound
  ## and the CID
  request <- paste0( pubchem_compound, cid, "/synonyms/JSON")
  
  content <- httr::content( GET( request))
  cid <- content$InformationList$Information[[1]]$CID
  
  synonyms <- unlist( content$InformationList$Information[[1]]$Synonym)
  
  return( list( cid=cid, synonyms=synonyms))
  
}


#' The function collects all IDs that are published at the Pubchem website.
#' By means of a RESTful HTTP request through the CID, alles descriptions, that are linked at the Pubchem are evaluated.
#' Based on these links, the IDs that are used at these sources can be retrieved.
#' @param cid Pubchem ID of the compound of interest.
#' @return A list containing the IDs, ID source, and a functional link to the compound page at this source.
get_description_by_cid <- function( cid){
  
  request <- paste0( pubchem_compound, cid, "/description/JSON")
  
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


#' The function returns all pubchem entries that are found with a given compound name.
#' By means of a RESTful HTTP request, all entries that are found with a string search are returned.
#' @param name Compound to be used for searching omics data sets.
#' @return A data.frame collecting CIDs and compound names.
get_cids_by_name <- function( name){
  
  ## create a http request which is used to retrieve all synonyms of a given compound
  ## and the CID
  
  name <- gsub( " ", "%20", name)
  
  request <- paste0( "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
                     name, "/synonyms/JSON?name_type=word")
  content <- httr::content( GET( request))
  
  if( "Fault" %in% names(content)){
    return( data.frame( cid = c(0000), name = c("EMPTY") ))
  }
  
  cid <- sapply( content$InformationList$Information, function(x) x$CID)
  pattern <- paste0( "^\\D*\\d\\D*", name, "\\D*$", "|" ,"^\\D*", name, "\\D*\\d\\D*$", "|", "^\\D*", name, "\\D*$")
  
  synonyms_list <- sapply( content$InformationList$Information, function(x) x$Synonym)
  synonyms <- sapply( synonyms_list, function(x){ x[grep( pattern, x, ignore.case = TRUE)[1]] })
  
  ## remove the NULL elements from synonyms and, beforehand, the respective IDs from the cid list
  cid <- cid[!sapply( synonyms, is.null)]
  synonyms[sapply( synonyms, is.null)] <- NULL
  
  if( length( cid) > 1){
    results <- data.frame( cid = cid, name = unlist( synonyms) )
  } else{
    results <- data.frame( cid = c(0000), name = c("EMPTY"))    
  }
  
  
  return( arrange( results, cid))
  
}


#' This function collects important compound-related information from pubchem based on the compound CID.
#' By means of a RESTful HTTP request retrieved information contain:
#' Molecular Formular, SMILE, InChI, and the Definition (including the source)
#' @param cid The Pubchem compound ID (CID)
#' @return Dataframe with two columns (property and value)
#' @examples
#' get_compound_information( cid = 702)
#' get_compound_information( 2244)
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



#' Search with a given compound name or compound ID in the Comptox database for overlapping hits to refine the search. 
#' @param compound The search string
#' @return Dataframe containing CIDs and Names of potentially matching compounds
get_comptox_by_input <- function( compound){
  
  ## check if the input name is a string or an ID
  ## in case no ID pattern matches, assume the compound is a string
  if( length( grep( "^DTX[SC]ID\\d+$", compound)) > 0){

        m <- unique( c( match( compound, comptox$DTXSID), match( compound, comptox$DTXCID)))
  
  } else if( length( grep( "^\\d+$", compound)) > 0){
    
      m <- match( compound, comptox$CID)
    
  } else if( length( grep( "^\\d+-\\d\\d-\\d$", compound)) > 0){
    
      ## check if cas number is correct
      if( !is_cas_correct( compound)){
        showModal( modalDialog(
          title = "CAS number ERROR",
          "The entered CAS number is not correct!",
          easyClose = TRUE
        ))
        
        return( data.frame( cid=c(0000), name=c("CAS")))
      }
      
      m <- match( compound, comptox$CAS)

  } else{
    
  
    ## assume the entered compound is a name
    ## -> pattern matching with columns NAME ans SYNONYM
    sids_comptox <- comptox$CTXSID[grep( compound, comptox$NAME, ignore.case = TRUE)]
    sids_synonyms <- unique( synonyms$DTXSID[ grep( compound, synonyms$SYNONYM, ignore.case = TRUE)])
    sids <- unique( c( sids_comptox, sids_synonyms)) 

    m <- match( sids, comptox$DTXSID)

  }

  m <- m[ !is.na( m)]
  if( length( m) > 0){
    df <- data.frame( cid = as.integer(comptox$CID[ m]), name = comptox$NAME[ m])
  } else{
    df <- data.frame( cid=c(0000), name=c("EMPTY"))
  }
  
  return( arrange( df, cid))
  
}



#' Retrieve all IDs from Comptox dashboard that are mapped to a given pubchem ID.
#' @param cid The input pubchem ID
#' @return A dataframe with 2 columns for ID source and ID.
#' @example get_comptox_ids_by_cid( 702)  
get_comptox_ids_by_cid <- function( cid){
  
  m <- match( cid, comptox$CID)
  m <- m[ !is.na( m)]
  
  if( length( m) == 0){
    return()
  }
  
  df <- data.frame( Database = c( "CAS-RN", "Comptox Dashboard"), IDs = c( comptox$CAS[ m], c( sprintf( '<a href=\"https://comptox.epa.gov/dashboard/%s\" target=\"_blank\">%s</a>', comptox$DTXSID[m], comptox$DTXSID[m]))))
  return( df)
  
}


  