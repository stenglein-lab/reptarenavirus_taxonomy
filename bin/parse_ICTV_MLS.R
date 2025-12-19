library(tidyverse)
library(readxl)


# read in ICTV Master Species List from 2025/10/13
# downloaded from: https://ictv.global/msl
msl <- read_excel("../input/VMR_MSL40.v2.20251013.xlsx", sheet = "VMR MSL40")
reptarenavirus_isolates <- filter(msl, Genus == "Reptarenavirus")
reptarenavirus_species_names <- reptarenavirus_isolates %>% group_by(Species) %>% summarize() %>% pull(Species)

# Pat Schloss's oxford comma function
# see e.g.: https://github.com/SchlossLab/Schloss_Rarefaction_mSphere_2024/blob/main/submission/manuscript.Rmd
oxford_comma <- function(x) {
  
  if(length(x) < 2){
    x
  } else if(length(x) == 2){
    paste(x, collapse = " and ")
  } else {
    paste(paste(x[-length(x)], collapse=", "), x[length(x)], sep=", and ")
  }
}

species_list <- oxford_comma(reptarenavirus_species_names)
species_abbr_list <- str_replace_all(species_list, "Reptarenavirus", "R.")
reptarenavirus_species <- list(
 n_spp = length(reptarenavirus_species_names),
 spp_names = species_list,
 spp_abbr_names = species_abbr_list)

saveRDS(reptarenavirus_species, "../RDS/reptarenavirus_species.RDS")

