library(tidyverse)

# combine reptarenavirus metadata from difference sources:
#
# 1. from genbank files of nucleotide sequences: accession, organism name, sequnce taxid, segment
# 2. from E-Utils: accession, species name, species taxid
#
# Mark Stenglein 11/6/2025

genbank_derived_map <- read.delim("../metadata/genbank_metadata.txt")
eutils_derived_map  <- read.delim("../metadata/accession_species_taxid_map.txt")

# read in spreadsheet containing table 1 from Hepojoki et al (2015), PMID: 26041290, re: multiply infected snakes
Hepojoki_table_1 <- read_excel("../metadata/Hepojoki_2015_table_1.xlsx")
Hepojoki_snakes <- select(Hepojoki_table_1, accession, Hepojoki_snake_number) 

head(genbank_derived_map)
head(eutils_derived_map)

# merge tables
reptarenavirus_accession_metadata <- left_join(genbank_derived_map, eutils_derived_map)


# make some labeling variables
reptarenavirus_accession_metadata <- reptarenavirus_accession_metadata %>% 
  mutate(acc_org = paste0(accession, "_", organism),
         acc_spp = paste0(accession, "_", species_name))

# make a boolean to indicate whether accession is classified
reptarenavirus_accession_metadata <- reptarenavirus_accession_metadata %>% 
  mutate(classified = case_when(
    species_name %in% c("unclassified Reptarenavirus") ~ F,
    .default = T) )

# parse out snake # (a la Stenglein et al 2015, doi:10.1371/journal.ppat.1004900) from sequence description
reptarenavirus_accession_metadata$Stenglein_snake_number <- 
  str_match(reptarenavirus_accession_metadata$description, "/snake(\\d+)/")[,2]

# join in multiple infection data from Hepojoki et al 2015, PMID 26041290
reptarenavirus_accession_metadata <- left_join(reptarenavirus_accession_metadata, Hepojoki_snakes)
reptarenavirus_accession_metadata <- mutate(reptarenavirus_accession_metadata, 
                                            snake_number = case_when(
                                              !is.na(Hepojoki_snake_number) ~ Hepojoki_snake_number,
                                              !is.na(Stenglein_snake_number) ~ paste0("Stenglein_snake_", Stenglein_snake_number),
                                              .default = NA
                                            ))


# how many sequences per snake?
n_seqs_per_snake <- reptarenavirus_accession_metadata %>% group_by(snake_number) %>% summarize(n=n()) %>% arrange(-n)
n_seqs_per_snake

# output metadata table
write.table(reptarenavirus_accession_metadata, "../metadata/reptarenavirus_accession_metadata.txt",
            sep="\t", row.names = F, quote=F)

saveRDS(reptarenavirus_accession_metadata, "../metadata/reptarenavirus_accession_metadata.RDS")


###################################
# New proposed species 
###################################

# new proposed species: 



###################################