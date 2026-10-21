#!/usr/bin/env python

# extract the species name and taxid for viral sequences in genbank format
# MDS 11/3/2025

import sys

from Bio import GenBank
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq

# genbank format file - first command-line arg
input_handle  = open(sys.argv[1], "r")
# second command line, a static segment character string that will be output in a "segment" column
segment       = sys.argv[2]

# print header line for tsv output
print ("accession", "organism", "taxid", "segment", "description", "isolate", "strain", sep="\t")

# for each sequence in file
for rec in SeqIO.parse(input_handle, "genbank") :

   if rec.features:
      for feature in rec.features :

         accession = rec.name
         description = rec.description
         organism = "NA"
         taxid = "NA"
         isolate = "NA"
         strain = "NA"

         if feature.type == "source": 
           organism = feature.qualifiers.get("organism")[0]
           taxid    = feature.qualifiers.get("db_xref")[0].replace("taxon:", "")
           if 'isolate' in feature.qualifiers.keys() and feature.qualifiers['isolate']:
              isolate  = feature.qualifiers.get("isolate")[0]
           if 'strain' in feature.qualifiers.keys() and feature.qualifiers['strain']:
              strain   = feature.qualifiers.get("strain")[0]
           # output accession -> organism map
           print (accession, organism, taxid, segment, description, isolate, strain, sep="\t")

#output_handle.close()
input_handle.close()

