#!/usr/bin/env python

import re
import sys
import argparse as P

parse = P.ArgumentParser( description = "Extract relevant information from SDF file and write an R readable comma-separated file.")
parse.add_argument( '--input', '-i', metavar=str, help='Specify the input SDF file.')
parse.add_argument( '--output', '-o', metavar=str, help='Specify the output tsv file.')


args = parse.parse_args()
count = 0
sid_value = ""
syn = False

outfile = open( args.output, "w")
outfile.write( "dsstox_substance_id\tsynonym\n")

with open( args.input) as f:

    for line in f:

        if line.startswith("DTXS"):
            sid_value = line.rstrip()
            continue
        if line.startswith(">  <Syn"):
            syn = True
            continue
        if syn:
            if line.rstrip() == "":
                syn = False
                continue
            outfile.write( sid_value + "\t" + line.rstrip() + "\n")
        

outfile.close()
