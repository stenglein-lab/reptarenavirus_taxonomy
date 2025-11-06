library(rentrez)
library(tidyverse)

# this script uses the NCBI E-Utils interface, via the rentrez package interface, 
# to identify nucleotide sequences associated with the reptarenavirus taxids
# 
# Mark Stenglein, 11/6/2025


# current ICTV reptarenavirus species and their NCBI taxids 
# could read this in from a file, but this is simple
reptarenavirus_species <- 
  tribble(~species_name,                  ~species_taxid,
          "Reptarenavirus aurei",         3052723,
          "Reptarenavirus californiae",   3052724,
          "Reptarenavirus commune",       3052725,
          "Reptarenavirus giessenae",     3052726,
          "Reptarenavirus rotterdamense", 3052727,
          "unclassified Reptarenavirus",  1654822)
          
# this function fetches the NCBI nucleotide IDs associated with a particular taxid
fetch_taxid_nucleotide_accession <- function(taxid, max_ids = 1000) {
  
  nt_search_results <- 
    entrez_search(db="nucleotide", term=paste0("txid", taxid, "[Organism:exp]"), retmax=max_ids, idtype="acc")
  
  # return the IDS as a character vector
  nt_search_results$ids
}

# download accessions for each taxid
nt_accessions_list <- lapply(reptarenavirus_species$species_taxid, fetch_taxid_nucleotide_accession)
# name the resulting list items
names(nt_accessions_list) <- reptarenavirus_species$species_taxid

# use dplyr to convert list of vectors to a dataframe
taxid_accessions_map <- 
  enframe(nt_accessions_list, name = "species_taxid", value = "accession") %>%
  unnest(accession)

# convert taxids to integers 
taxid_accessions_map$species_taxid = as.numeric(taxid_accessions_map$species_taxid)

# ESearch returns accessions as accession.version
# keep just accession and get rid of version
taxid_accessions_map$accession <-
  str_match(taxid_accessions_map$accession, "(\\S+)\\.([0-9])+")[,2]

# join in species names
taxid_accessions_map <- left_join(taxid_accessions_map, reptarenavirus_species) 
              
# write out the map
write.table(taxid_accessions_map, "../metadata/accession_species_taxid_map.txt",
            quote=F, sep="\t", row.names=F)
