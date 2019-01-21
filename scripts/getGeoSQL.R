
args = commandArgs( trailingOnly = TRUE)
if( length( args) != 1){
  stop( "Exactly one argument (destination directory) had to be passed to the script.", call. = FALSE)
}

## get GeoMetaDB
library( GEOmetadb, quietly = TRUE)

## get the actual SQL file and store it at the given directory
getSQLiteFile( destdir = args[1])

