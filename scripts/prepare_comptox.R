
args = commandArgs( trailingOnly = TRUE)
if( length( args) != 1){
  
  stop( "Error in prepareCTDfile.R: Exactly one argument has to be given!", call. = FALSE)
  
}

#path <- "../data/"

## read all preprocessed comptox files and adapt the headers
dtxcid <- read.csv( file = paste0( args[1], "dtxcid.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
# dtxcid <- read.csv( file = paste0( path, "dtxcid.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")

dtxsid <- read.csv( file = paste0( args[1], "dtxsid.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
# dtxsid <- read.csv( file = paste0( path, "dtxsid.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")

url <- read.csv( file = paste0( args[1], "url.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
# url <- read.csv( file = paste0( path, "url.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")

cas <- read.csv( file = paste0( args[1], "cas.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
# cas <- read.csv( file = paste0( path, "cas.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")

names <- read.csv( file = paste0( args[1], "names.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")
# names <- read.csv( file = paste0( path, "names.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote = "")

sdf <- cbind( dtxcid, dtxsid, url, cas, names)
colnames( sdf) <- c("DTXCID", "DTXSID", "URL", "CAS", "NAME")

pubchem <- read.csv( file = paste0( args[1], "Pubchem_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
# pubchem <- read.csv( file = paste0( path, "Pubchem_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

inchi <- read.csv( file = paste0( args[1], "InChI_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
# inchi <- read.csv( file = paste0( path, "InChI_mapping.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
colnames( inchi) <- c( "DTXSID", "InChI_String", "InChI_Key")

synonyms <- read.csv( file = paste0( args[1], "Synonyms.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
# synonyms <- read.csv( file = paste0( path, "Synonyms.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")
colnames( synonyms) <- c( "DTXSID", "SYNONYM")


## remove all entries without a possible mapping 
comptoxpub <- match( pubchem$DTXSID, sdf$DTXSID)
rm_pubchem <- which( is.na(comptoxpub) == TRUE)
pubchem <- pubchem[ -rm_pubchem, ]
comptoxpub <- comptoxpub[ !is.na( comptoxpub)]

## combine cas and pubchem to create a ID mapping table including
## CAS CID SID CTXSID NAME
comptox <- data.frame( CAS = sdf$CAS[ comptoxpub], CID = pubchem$CID, SID = pubchem$SID, DTXSID = pubchem$DTXSID, NAME = sdf$NAME[ comptoxpub], stringsAsFactors = FALSE)

## remove those entries where a synonym is present but no mapping to the df data table is possible
temp <- match( synonyms$DTXSID, sdf$DTXSID)
rm_syn_list <- which( is.na( temp))
synonyms <- synonyms[ -rm_syn_list, ]

## remove those sdf mappings with no connection to the df table
temp <- match( sdf$DTXSID, comptox$DTXSID)
rm_sdf_list <- which( is.na( temp))
sdf <- sdf[ -rm_sdf_list, ]

## remove InChI entries without mapping to DTXSID
temp <- match( inchi$DTXSID, sdf$DTXSID)
rm_inchi_list <- which( is.na( temp))
inchi <- inchi[ -rm_inchi_list, ]


save( comptox, sdf, synonyms, inchi, file = paste0( args[1], "comptox.RData"))
