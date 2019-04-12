#!/usr/bin/env python

import re
import sys
import argparse as P

parse = P.ArgumentParser( description = "Extract relevant information from SDF file and write an R readable comma-separated file.")
parse.add_argument( '--input', '-i', metavar=str, help='Specify the input SDF file.')
parse.add_argument( '--output', '-o', metavar=str, help='Specify the output tsv file.')


args = parse.parse_args()
count = 0
name = False
name_value = ""
cid = False
cid_value = ""
sid = False
sid_value = ""
cas = False
cas_value = ""
syn = False

outfile = open( args.output, "w")
outfile.write( "dsstox_structure_id\tdsstox_substance_id\tcasrn\tpreferred_name\tsynonym\n")

with open( args.input) as f:

    for line in f:

        if re.search( "DSSTox_Substance_Id", line):
            sid = True
            continue
        if sid:
            sid_value = line.rstrip()
            sid = False
            continue
        if re.search( "DSSTox_Structure_Id", line):
            cid = True
            continue
        if cid:
            cid_value = line.rstrip()
            cid = False
            continue
        if re.search( "<CAS-RN>", line):
            cas = True
            continue
        if cas:
            cas_value = line.rstrip()
            cas = False
            continue
        if re.search( "<Preferred_Name>", line):
            name = True
            continue
        if name:
            name_value = line.rstrip()
            name = False
            continue
        if re.search( "<Synonyms>", line):
            syn = True
            continue
        if syn:
            if line.rstrip() == "":
                syn = False
                continue
            outfile.write( cid_value + "\t" + sid_value + "\t" + cas_value + "\t" + name_value + "\t" + line.rstrip() + "\n")
        

outfile.close()
