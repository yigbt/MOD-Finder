#!/usr/bin/env bash

dir=$1
cas="Cas_mapping"
pubchem="Pubchem_mapping"
inchi="InChI_mapping"


#######
## STEP 0
## download the comptox download page to extract the current downloadable files
#######

wget -O "comptox_download.html" https://comptox.epa.gov/dashboard/downloads



#######
## STEP 1
## Download necessary files from comptox
#######

## get the current ftp link of the cas number file
url=`grep "DSSTox Identifier to PubChem Identifier Mapping File" comptox_download.html | sed 's/^.*\(ftp:\/\/.*txt\).*$/\1/'`
wget -O $dir$pubchem".txt" $url
mv $dir$pubchem".txt" $dir$pubchem".tsv"


## get the current ftp link of the cas number file
url=`grep ">DSSTox identifiers mapped to CAS Numbers and Names File" comptox_download.html | sed 's/^.*\(ftp:\/\/.*xlsx\).*$/\1/'`
wget -O $dir$cas".xlsx" $url
ssconvert -O 'separator="	" format=raw' $dir$cas".xlsx" $dir$cas".txt"
rm $dir$cas".xlsx"
mv $dir$cas".txt" $dir$cas".tsv"
sed -i 's/"//g' $dir$cas".tsv"


## get the current ftp link of the synonyms file
url=`grep "DSSTox Synonyms File" $dir"comptox_download.html" comptox_download.html | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
wget $url
unzip $file
rm $file
mv ${file::-4}".sdf" $dir"Synonyms.sdf"


## get the current link of the SDF file
## get the filename from the link
url=`grep "DSSTox SDF File" $dir"comptox_download.html" comptox_download.html | sed 's/^.*\(ftp:\/\/.*sdf.gz\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
wget $url
gzip -d $file
mv ${file::-3} $dir"SDF_file.sdf"


## extract the current mapping file from the download page (ftp link)
## get the filename from the link
## produce the unzipped name of the file, which contains the date of its generation
url=`grep "DSSTox Mapping File" $dir"comptox_download.html" comptox_download.html | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
unzipped=`echo ${file::-4} | sed 's/\(.*_\)\(.*\)$/dsstox_\2.tsv/'`
wget $url
unzip $file
echo -e "dsstox_substance_id\tinchi_string\tinchi_key" > $dir$inchi".tsv"
cat $unzipped >> $dir$inchi".tsv"
rm $file $unzipped



#######
## STEP 2
## Extract the needed information from SDF files
#######
grep -A 1 PUBCHEM_EXT_DATASOURCE $dir"SDF_file.sdf" | grep DTXSID > dtxsid.txt
grep -A 1 PUBCHEM_EXT_SUBSTANCE_URL $dir"SDF_file.sdf" | grep https > url.txt
grep -A 1 DSSTox_Structure_Id $dir"SDF_file.sdf" | grep DTXCID > dtxcid.txt
# grep -A 1 PUBCHEM_SUBSTANCE_SYNONYM $dir"SDF_file.sdf" | grep -v PUBCHEM | grep -v "\-\-" > synonyms.txt
# grep -A 2 PUBCHEM_SUBSTANCE_SYNONYM ../data/SDF_file.sdf  | grep -P "^\d+-\d+" > cas.txt

echo -e "dsstox_structure_id\tdsstox_substance_id\turl" > $dir"SDF_file.tsv"
paste dtxcid.txt dtxsid.txt url.txt | column -s $'\t' >> $dir"SDF_file.tsv"
rm dtxsid.txt dtxcid.txt url.txt $dir"SDF_file.sdf"


#######
## STEP 3
## Extract the synonyms from synonyms.sdf file
#######
./sdf2matrix.py -i $dir"Synonyms.sdf" -o $dir"Synonyms.tsv"
rm $dir"Synonyms.sdf"
rm comptox_download.html


#######
## STEP 4
## Combine the mappings and save the data tables in an R readable-format
#######
Rscript prepare_comptox.R $dir
