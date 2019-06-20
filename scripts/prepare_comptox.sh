#!/usr/bin/env bash

dir=$1
data=$dir"data/"
scripts=$dir"scripts/"
cas="Cas_mapping"
pubchem="Pubchem_mapping"
inchi="InChI_mapping"
log="/var/log/shiny-server.log"



#######
## STEP 0
## download the comptox download page to extract the current downloadable files
#######
echo `date`"Download comptox information file" >> $log
wget -O $data"comptox_download.html" https://comptox.epa.gov/dashboard/downloads



#######
## STEP 1
## Download necessary files from comptox
#######

## get the current ftp link of the cas number file
echo `date`"Process comptox-pubchem identifier" >> $log
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


### get the current ftp link of the cas number file
#echo `date`"Process comptox-CASnr identifier" >> $log
#url=`grep ">DSSTox identifiers mapped to CAS Numbers and Names File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*xlsx\).*$/\1/'`
#wget -O $data$cas".xlsx" $url
#ssconvert -O 'separator="	" format=raw' $data$cas".xlsx" $data$cas".txt"
#rm $data$cas".xlsx"
#mv $data$cas".txt" $data$cas".tsv"
#sed -i 's/"//g' $data$cas".tsv"

## head Cas_mapping.tsv --> no tabs inserted
## 60-35-5DTXSID7020005AcetamideDTXCID505InChI=1S/C2H5NO/c1-2(3)4/h1H3,(H2,3,4)DLFVBJFMPXGRIB-UHFFFAOYSA-N
## 616-91-1DTXSID5020021N-Acetyl-L-cysteineDTXCID4021InChI=1S/C5H9NO3S/c1-3(7)6-4(2-10)5(8)9/h4,10H,2H2,1H3
## 8052-16-2DTXSID4020030Actinomycin C
## 50-76-0DTXSID9020031Actinomycin DDTXCID5031InChI=1S/C62H86N12O16/c1-27(2)42-59(84)73-23-17-19-36(73)57(8


## get the current ftp link of the synonyms file
echo `date`"Process comptox synonyms file" >> $log
url=`grep "DSSTox Synonyms File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
wget $url
unzip $file
rm $file
mv ${file::-4}".sdf" $data"Synonyms.sdf"


## get the current link of the SDF file
## get the filename from the link
echo `date`"Process comptox SDF file" >> $log
url=`grep "DSSTox SDF File" $data"comptox_download.html" | sed 's/^.*\(ftp:\/\/.*zip\).*$/\1/'`
file=`echo $url | sed 's/.*\///'`
wget $url
unzip $file
cat ${file::-4}_* > $data"SDF_file.sdf"
rm $file ${file::-4}_*


## extract the current mapping file from the download page (ftp link)
## get the filename from the link
## produce the unzipped name of the file, which contains the date of its generation
echo `date`"Process comptox mapping file" >> $log
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
echo `date`"Process SDF file to get SID, CID, URL, and CAS mapping" >> $log
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
echo `date`"Process downloaded synonyms file to get tsv format" >> $log
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
echo `date`"Process the gathered comptox information into RData file" >> $log
Rscript $scripts"prepare_comptox.R" $data


#######
## STEP 5
## Cleaning
#######
rm $data"dtxsid.txt" $data"dtxcid.txt" $data"url.txt" $data"cas.txt" $data"SDF_file.sdf"
rm $data"Synonyms.sdf"
rm $data"comptox_download.html"
