library(tidyverse)
library(seqinr)
library(MSA2dist)

# Read in assignment of accession (tip label) -> NCBI-defined taxon
metadata <- readRDS("../metadata/reptarenavirus_accession_metadata.RDS")

# read in alignments
L_MSA   <- read.alignment("../make_alignments_and_trees/results/alignments/L_CDS.mafft.fasta", format = "fasta")
NP_MSA  <- read.alignment("../make_alignments_and_trees/results/alignments/NP_CDS.mafft.fasta", format = "fasta")

# convert to pairwise distances using dist.alignment
L_dist   <- as.matrix(dist.alignment(L_MSA, matrix="identity"))
# square each element of the matrix because dist.alignment returns square root of distance
L_dist   <- L_dist ^ 2 
NP_dist  <- as.matrix(dist.alignment(NP_MSA, matrix="identity"))
NP_dist  <- NP_dist ^ 2 

# calculate evolutionary distances (instead of p-distances) using MSA2dist
L_MSA_ds    <- aln2dnastring(L_MSA)
L_dist_TN93 <- dnastring2dist(L_MSA_ds, model="TN93", threads=8)

# how close is a sequence to the closest classified sequence?
closest_classified_sequence <- function(an_accession, dist_matrix) {
  
  debug <- 0
  if (debug) {
    an_accession <- "MW091469"
    dist_matrix  <- L_dist
    # dist_matrix  <- L_dist_TN93$distSTRING
  }
  
  # from dist.alignment documentation: 
  # The resulting matrix contains the squared root of the pairwise distances. 
  # For example, if identity between 2 sequences is 80 the squared root of (1.0 - 0.8) i.e. 0.4472.
  this_acc_distances  <- dist_matrix [an_accession, ]
  
  # convert to long format and square root distance to identity
  this_acc_identities <- 
    tibble(pairwise_identity = as.numeric(1 - this_acc_distances),
           accession         = names(this_acc_distances))
  
  # join in metadata
  this_acc_identities <- left_join(this_acc_identities, metadata)
  
  # find closest classified
  nearest_classified <- filter(this_acc_identities,  accession != an_accession & classified == T) %>% 
    arrange(-pairwise_identity)  
  
  # return the first (nearest)
  nearest_classified %>% filter(row_number() == 1)
}

closest_MN567045_relative <- closest_classified_sequence("MN567045", L_dist)
closest_MN567045_relative$this_accession <- "MN567045"

closest_MW091469_relative <- closest_classified_sequence("MW091469", L_dist)
closest_MW091469_relative$this_accession <- "MW091469"

saveRDS(closest_MN567045_relative, "../RDS/closest_MN567045_relative.RDS")
saveRDS(closest_MW091469_relative, "../RDS/closest_MW091469_relative.RDS")
