
args = commandArgs( trailingOnly = TRUE)
if( length( args) != 1){
  
  stop( "Error in prepareCTDfile.R: Exactly one argument has to be given!", call. = FALSE)
  
}

## read all preprocessed comptox files and adapt the headers
cas <- read.csv( file = paste0( args[1], "Cas_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames( cas) <- c("CAS", "DTXSID", "NAME")
sdf <- read.csv( file = paste0( args[1], "SDF_file.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames( sdf) <- c("DTXCID", "DTXSID", "URL")
pubchem <- read.csv( file = paste0( args[1], "Pubchem_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
inchi <- read.csv( file = paste0( args[1], "InChI_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames( inchi) <- c( "DTXSID", "InChI_String", "InChI_Key")
synonyms <- read.csv( file = paste0( args[1], "Synonyms.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames( synonyms) <- c( "DTXCID", "DTXSID", "CAS", "NAME", "SYNONYM")


## remove all entries without a possible mapping 
caspub <- match( pubchem$DTXSID, cas$DTXSID)
rm_pubchem <- which( is.na(caspub) == TRUE)
pubchem <- pubchem[ -rm_pubchem, ]
caspub <- caspub[ !is.na( caspub)]

## combine cas and pubchem to create a ID mapping table including
## CAS CID SID CTXSID NAME
comptox <- data.frame( CAS = cas$CAS[ caspub], CID = pubchem$CID, SID = pubchem$SID, DTXSID = pubchem$DTXSID, NAME = cas$NAME[ caspub], stringsAsFactors = FALSE)

## remove those entries where a synonym is present but no mapping to the df data table is possible
temp <- match( synonyms$DTXSID, comptox$DTXSID)
rm_syn_list <- which( is.na( temp))
synonyms <- synonyms[ -rm_syn_list, ]
synonyms <- synonyms[, -c(1,3,4)]

## remove those sdf mappings with no connection to the df table
temp <- match( sdf$DTXSID, comptox$DTXSID)
rm_sdf_list <- which( is.na( temp))
sdf <- sdf[ -rm_sdf_list, ]

## remove InChI entries without mapping to DTXSID
temp <- match( inchi$DTXSID, comptox$DTXSID)
rm_inchi_list <- which( is.na( temp))
inchi <- inchi[ -rm_inchi_list, ]


save( comptox, sdf, synonyms, inchi, file = paste0( args[1], "comptox.RData"))
