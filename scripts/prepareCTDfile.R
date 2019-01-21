
library( CTDquerier)

args = commandArgs( trailingOnly = TRUE)
if( length( args) != 1){
  
  stop( "Error in prepareCTDfile.R: Exactly one argument has to be given!", call. = FALSE)
  
}

## Download the CTD_chemicals file from CTDbase.org
## parse the file and create an RData file containing a data.frame with an appropriate header
prepare_CTD_database_file <- function( outdir){
  
  filename <- paste0( outdir, "CTD_chemicals.tsv.gz")
  
  if( !file.exists( filename)){
    download_ctd_chem( filename = filename)
  }else if( Sys.time() - file.info( filename)$mtime > 6 ){
    file.remove( filename)
    download_ctd_chem( filename = filename)
  }
  
  if( !file.exists(filename)){
    
    cat( "ERROR: Could not download CTD_Chemicals file!\n")
    return( NULL)  
  }
  
  ctd <- read.csv( filename, header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE )
  colnames( ctd) <- c( "Chemical.Name", "Chemical.ID", "CasRN", "Definition", "ParentIDs", "TreeNumbers", "ParentTreeNumbers", "Synonyms", "DrugBankIDs")
  
  save( ctd, file = paste0( outdir, "CTD_chemicals.RData"))
  
}

## download the CTD file and save the file in an RData file
prepare_CTD_database_file( args[1])

