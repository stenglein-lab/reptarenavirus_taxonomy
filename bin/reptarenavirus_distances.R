library(tidyverse)
library(seqinr)
library(MSA2dist)
library(patchwork)

# Read in assignment of accession (tip label) -> NCBI-defined taxon
metadata <- readRDS("../metadata/reptarenavirus_accession_metadata.RDS")

# read in alignments
# nt alignments
L_MSA   <- read.alignment("../make_alignments_and_trees/results/alignments/L_CDS.mafft.fasta", format = "fasta")
NP_MSA  <- read.alignment("../make_alignments_and_trees/results/alignments/NP_CDS.mafft.fasta", format = "fasta")

# protein alignments
L_prot_MSA   <- read.alignment("../make_alignments_and_trees/results/alignments/L_prot.mafft.fasta", format = "fasta")
NP_prot_MSA  <- read.alignment("../make_alignments_and_trees/results/alignments/NP_prot.mafft.fasta", format = "fasta")

# convert to pairwise distances using dist.alignment
# square each element of the matrix because dist.alignment returns square root of distance

# note that these are just basic pairwise distances, which under-report true evolutionary distance
# especially for distantly related nucleotide sequences, because of saturating mutations

# nt alignments
L_dist   <- as.matrix(dist.alignment(L_MSA, matrix="identity"))
L_dist   <- L_dist ^ 2 
NP_dist  <- as.matrix(dist.alignment(NP_MSA, matrix="identity"))
NP_dist  <- NP_dist ^ 2 

# prot alignments
L_prot_dist   <- as.matrix(dist.alignment(L_prot_MSA, matrix="identity"))
L_prot_dist   <- L_prot_dist ^ 2 
NP_prot_dist  <- as.matrix(dist.alignment(NP_prot_MSA, matrix="identity"))
NP_prot_dist  <- NP_prot_dist ^ 2 


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


# distances between example reassortant sequences 
example_reassortant_S <- c(
NP_dist["KP071610", "KP071539"],
NP_dist["KP071610", "KP071606"],
NP_dist["KP071539", "KP071606"])

# what is the minimum shared identity between example S segs
min_example_S_dist <- min(1-example_reassortant_S)
max_example_S_dist <- max(1-example_reassortant_S)
min_example_S_dist
max_example_S_dist

example_reassortant_L <- c(
L_dist["KP071611", "KP071540"],
L_dist["KP071611", "KP071607"],
L_dist["KP071540", "KP071607"])

min_example_L_dist <- min(1-example_reassortant_L)
max_example_L_dist <- max(1-example_reassortant_L)
min_example_L_dist
max_example_L_dist

saveRDS(min_example_S_dist, "../RDS/min_example_S_dist.RDS")
saveRDS(min_example_S_dist, "../RDS/max_example_S_dist.RDS")
saveRDS(min_example_L_dist, "../RDS/min_example_L_dist.RDS")
saveRDS(max_example_L_dist, "../RDS/max_example_L_dist.RDS")


# Distances in new species

# accessions of R. califronae sequences
R_californae_L_accessions <- c( "KP071531", "KP071533", "NC_018484", "KP071529" )
R_californae_S_accessions <- c( "NC_018481", "KP071528", "KP071530", "KP071532" )
R_californae_accessions   <- c( R_californae_L_accessions, R_californae_S_accessions )

# turn protein distance matrices into tidy tibbles

# NA out duplicate and diagonal pairwise distance values to avoid repeating and self comparison
NP_dist_tri <- NP_prot_dist
NP_dist_tri[upper.tri(NP_dist_tri, diag=T)] <- NA
NP_dist_long <- 
  NP_dist_tri %>% 
  as_tibble(rownames = "acc_1")  %>%
  pivot_longer(cols = -acc_1, names_to = "acc_2", values_to = "distance") %>%
  filter(!is.na(distance))

L_dist_tri <- L_prot_dist
L_dist_tri[upper.tri(L_dist_tri, diag=T)] <- NA
L_dist_long <- 
  L_dist_tri %>% 
  as_tibble(rownames = "acc_1")  %>%
  pivot_longer(cols = -acc_1, names_to = "acc_2", values_to = "distance") %>%
  filter(!is.na(distance))

# assign accessions to new species
NP_dist_long <- mutate(
  NP_dist_long,
  acc_1_species = if_else(acc_1 %in% R_californae_accessions, "Reptarenavirus_californae", "Reptarenavirus_giessenae" ),
  acc_2_species = if_else(acc_2 %in% R_californae_accessions, "Reptarenavirus_californae", "Reptarenavirus_giessenae" ),
  same_species = acc_1_species == acc_2_species,
  same_species_label = if_else(same_species, "Same species", "Different species")
)

L_dist_long <- mutate(
  L_dist_long,
  acc_1_species = if_else(acc_1 %in% R_californae_accessions, "Reptarenavirus_californae", "Reptarenavirus_giessenae" ),
  acc_2_species = if_else(acc_2 %in% R_californae_accessions, "Reptarenavirus_californae", "Reptarenavirus_giessenae" ),
  same_species = acc_1_species == acc_2_species,
  same_species_label = if_else(same_species, "Same species", "Different species")
)

NP_dist_same_spp <- filter(NP_dist_long, same_species)
NP_dist_diff_spp <- filter(NP_dist_long, !same_species)

L_dist_same_spp <- filter(L_dist_long, same_species)
L_dist_diff_spp <- filter(L_dist_long, !same_species)

min(NP_dist_diff_spp$distance)
max(NP_dist_diff_spp$distance)
min(NP_dist_same_spp$distance)
max(NP_dist_same_spp$distance)

min(L_dist_diff_spp$distance)
max(L_dist_diff_spp$distance)
min(L_dist_same_spp$distance)
max(L_dist_same_spp$distance)

L_dist_p <- 
  ggplot(L_dist_long) +
  geom_histogram(aes(x=100*(1-distance), fill=same_species), binwidth=1) +
  facet_wrap(~same_species_label, ncol=1, scales="free_y") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("coral3", "slateblue")) +
  xlab("Pairwise L protein amino acid identity (%)") +
  ylab("Number pairwise comparisons")
L_dist_p

NP_dist_p <- 
  ggplot(NP_dist_long) +
  geom_histogram(aes(x=100*(1-distance), fill=same_species), binwidth=1) +
  facet_wrap(~same_species_label, ncol=1, scales="free_y") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("coral3", "slateblue")) +
  xlab("Pairwise NP amino acid identity (%)") +
  ylab("Number pairwise comparisons")
NP_dist_p

dist_p <- 
  L_dist_p + NP_dist_p + plot_layout(ncol = 1) + plot_annotation(tag_levels = 'A')
dist_p

ggsave(dist_p, file=paste0(output_dir, "pairwise_distance_figure.pdf"), 
       units="in", width=7.5, height=10)
ggsave(dist_p, file=paste0(output_dir, "pairwise_distance_figure.png"), 
       units="in", width=7.5, height=10, device=grDevices::png)

