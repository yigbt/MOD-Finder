#!/usr/bin/env bash

scripts=$1
data=$2

cas="Cas_mapping"
pubchem="Pubchem_mapping"
inchi="InChI_mapping"



#######
## STEP 0
## download the comptox download page to extract the current downloadable files
#######
echo `data`" Download comptox information file"
wget -O $data"comptox_download.html" https://comptox.epa.gov/dashboard/downloads



#######
## STEP 1
## Download necessary files from comptox
#######

## get the current ftp link of the cas number file
echo `date`" Process comptox-pubchem identifier"
url=`grep "DSSTox Identifier to PubChem Identifier Mapping File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*txt\).*$/\1/'`
wget -O $data$pubchem".txt" $url
mv $data$pubchem".txt" $data$pubchem".tsv"

## head Pubchem_mapping.tsv
## SID     CID     DTXSID
## 316388891       20404   DTXSID30873143
## 316388890       10142816        DTXSID70873142
## 316388889       50742127        DTXSID40873139
## 316388888       19073841        DTXSID20873137
## 316388887       11505215        DTXSID00873135


## get the current ftp link of the synonyms file
echo `date`" Process comptox synonyms file"
url=`grep "DSSTox Synonyms File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
wget $url
unzip $file
rm $file
mv ${file::-4}".sdf" $data"Synonyms.sdf"


## get the current link of the SDF file
## get the filename from the link
echo `date`" Process comptox SDF file"
url=`grep "DSSTox SDF File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
wget $url
unzip $file
cat ${file::-4}_* > $data"SDF_file.sdf"
rm $file ${file::-4}_*


## extract the current mapping file from the download page (ftp link)
## get the filename from the link
## produce the unzipped name of the file, which contains the date of its generation
echo `date`" Process comptox mapping file"
url=`grep "DSSTox Mapping File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
unzipped=`echo ${file::-4} | sed 's/\(.*_\)\(.*\)$/dsstox_\2.tsv/'`
wget $url
unzip $file
echo -e "dsstox_substance_id\tinchi_string\tinchi_key" > $data$inchi".tsv"
cat $unzipped >> $data$inchi".tsv"
rm $file $unzipped

## head $data"/InChI_mapping.tsv"
## dsstox_substance_id	inchi_string	inchi_key
## DTXSID7020001	InChI=1S/C11H9N3/c12-10-6-5-8-7-3-1-2-4-9(7)13-11(8)14-10/h1-6H,(H3,12,13,14)	FJTNLJLPLJDTRM-UHFFFAOYSA-N
## DTXSID5039224	InChI=1S/C2H4O/c1-2-3/h2H,1H3	IKHGUXGNUITLKF-UHFFFAOYSA-N
## DTXSID50872971	InChI=1S/C4H8N2O/c1-3-5-6(2)4-7/h3-4H,1-2H3/b5-3+	IMAGWKUTFZRWSB-HWKANZROSA-N
## DTXSID2020004	InChI=1S/C2H5NO/c1-2-3-4/h2,4H,1H3/b3-2+	FZENGILVLUJGJX-NSCUHMNNSA-N


#######
## STEP 2
## Extract the needed information from SDF files
#######
echo `date`" Process SDF file to get SID, CID, URL, and CAS mapping"
grep -A 1 '<DSSTox_Compound_id>' $data"SDF_file.sdf" | grep DTXCID > $data"dtxcid.txt"
grep -A 1 '<DSSTox_Substance_id>' $data"SDF_file.sdf" | grep DTXSID > $data"dtxsid.txt"
grep -A 1 '<Dashboard_URL>' $data"SDF_file.sdf" | grep http > $data"url.txt"
grep -A 1 '<CASRN>' $data"SDF_file.sdf" | grep -v CASRN | grep -v "\-\-" > $data"cas.txt"
grep -A 1 "<Preferred_name>" $data"SDF_file.sdf" | grep -v "<Preferred_name>" | grep -v "^\-\-" | sed 's/^\s$/NA/' > $data"names.txt" 
#echo -e "dsstox_structure_id\tdsstox_substance_id\turl\tcas-rn" > $data"SDF_file.tsv"
#paste dtxcid.txt dtxsid.txt url.txt cas.txt | column -s $'\t' >> $data"SDF_file.tsv"
#rm dtxsid.txt dtxcid.txt url.txt cas.txt $data"SDF_file.sdf"


#######
## STEP 3
## Extract the synonyms from synonyms.sdf file
#######
echo `date`" Process downloaded synonyms file to get tsv format"
python $scripts"sdf2matrix.py" -i $data"Synonyms.sdf" -o $data"Synonyms.tsv"

## head $data"Synonyms.tsv"
## dsstox_structure_id     dsstox_substance_id     casrn   preferred_name  synonym
## DTXCID101       DTXSID7020001   26148-68-5      A-alpha-C       Amino-alpha-carboline
## DTXCID101       DTXSID7020001   26148-68-5      A-alpha-C       2-Amino-alpha-carboline
## DTXCID101       DTXSID7020001   26148-68-5      A-alpha-C       BRN 0744369
## DTXCID101       DTXSID7020001   26148-68-5      A-alpha-C       UNII-P0GZ1ICS6X
## DTXCID101       DTXSID7020001   26148-68-5      A-alpha-C       2-Amino-9H-pyrido[2,3-b]indole
## DTXCID404       DTXSID2020004   107-29-9        Acetaldehyde oxime      Acetaldoxime
## DTXCID404       DTXSID2020004   107-29-9        Acetaldehyde oxime      4-01-00-03121
## DTXCID404       DTXSID2020004   107-29-9        Acetaldehyde oxime      Aldoxime
## DTXCID404       DTXSID2020004   107-29-9        Acetaldehyde oxime      BRN 1209252
## DTXCID404       DTXSID2020004   107-29-9        Acetaldehyde oxime      EINECS 203-479-6


#######
## STEP 4
## Combine the mappings and save the data tables in an R readable-format
#######
echo `date`" Process the gathered comptox information into RData file"
Rscript $scripts"prepare_comptox.R" $data


#######
## STEP 5
## Cleaning
#######
rm $data"dtxsid.txt" $data"dtxcid.txt" $data"url.txt" $data"cas.txt" $data"SDF_file.sdf" $data"names.txt"
rm $data"Synonyms.sdf"
rm $data"comptox_download.html"
