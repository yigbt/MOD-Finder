#!/usr/bin/env bash

data=$1
scripts=`echo $0 | sed 's/\/create_local_databases.sh/\//'`

echo $scripts

if [ -z "$data" ]
then
    echo "USAGE: ./create_local_databases.sh output-path"
    echo
    echo "ERROR: Please specify a data directory"
    echo 
    exit
fi

if [ ! -d $data ]
then
    echo "USAGE: ./create_local_databases.sh output-path"
    echo
    echo "ERROR: No such directory!"
    echo 
    exit
fi

## check if the provided data directory contains a slash at the end
## in case there is no slash - add one...
case "$data" in
    */)
	;;
    *)
	data=$data"/"
	;;
esac


echo `date`" Start updating Database information for MOD-Finder"

echo `date`" Start GEOmetadb download"
Rscript $scripts"getGeoSQL.R" $data
echo `date`" End GEOmetadb download"

echo `date`" Start CTD_chemicals download"
Rscript $scripts"prepareCTDfile.R" $data
echo `date`" End CTD_chemicals download"

echo `date`" Start Comptox download"
bash $scripts"prepare_comptox_local.sh" $scripts $data
echo `date`" End Comptox download"
