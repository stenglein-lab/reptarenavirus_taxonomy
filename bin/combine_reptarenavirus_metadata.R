library(tidyverse)

# combine reptarenavirus metadata from difference sources:
#
# 1. from genbank files of nucleotide sequences: accession, organism name, sequnce taxid, segment
# 2. from E-Utils: accession, species name, species taxid
#
# Mark Stenglein 11/6/2025

genbank_derived_map <- read.delim("../metadata/genbank_metadata.txt")
eutils_derived_map  <- read.delim("../metadata/accession_species_taxid_map.txt")

head(genbank_derived_map)
head(eutils_derived_map)

# merge tables
reptarenavirus_accession_metadata <- left_join(genbank_derived_map, eutils_derived_map)

# output
write.table(reptarenavirus_accession_metadata, "../metadata/reptarenavirus_accession_metadata.txt",
            sep="\t", row.names = F, quote=F)
