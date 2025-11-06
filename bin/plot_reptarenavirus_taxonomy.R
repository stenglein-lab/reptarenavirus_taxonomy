library(tidyverse)
library(ggtree)
library(phytools)
# library(ape)
library(treeio)

# read in trees
L_tree  <- read.newick("../make_alignments_and_trees//results/trees/L_CDS.mafft.fasta.treefile")
NP_tree <- read.newick("../make_alignments_and_trees//results/trees/NP_CDS.mafft.fasta.treefile")

# midpoint root trees
L_midpoint  <- midpoint_root(L_tree)
NP_midpoint <- midpoint_root(NP_tree)

# Read in assignment of accession (tip label) -> NCBI-defined taxon
metadata <- read.delim("../metadata/reptarenavirus_accession_metadata.txt")

# make some labeling variables
metadata <- metadata %>% 
  mutate(acc_org = paste0(accession, "_", organism),
         acc_spp = paste0(accession, "_", species_name))

# make a boolean to indicate whether accession is classified
metadata <- metadata %>% 
  mutate(classified = case_when(
           species_name %in% c("unclassified Reptarenavirus") ~ F,
           .default = T) )

# how many sequences for of each species / segment?
metadata <- metadata %>% group_by(species_taxid, segment) %>% mutate(n_representatives=n())

# how many classified / unclassified?
classified_counts <- metadata %>% group_by(segment, classified) %>% summarize(n=n())
saveRDS(classified_counts, "../RDS/classified_counts.RDS")

L_metadata <- filter(metadata, segment == "L")
S_metadata <- filter(metadata, segment == "S")

# a shared graphical theme for plots
shared_tree_theme <- function() {
  theme_tree(base_size = 11) +
    theme()
}

# make a plot of all classified/unclassified sequences
L_classified_tree <- 
  ggtree(L_midpoint, ladderize = F) %<+% L_metadata +
  geom_tippoint(aes(fill = classified), 
                size=1.5, color="black", stroke=0.2, shape=21) +
  scale_fill_manual(values=c("white", "red")) +
  shared_tree_theme() +
  theme(legend.position = "none") +
  geom_treescale(linesize=0.5, x = 0, y = -15, fontsize=3, offset=2) +
  ggtitle("L")


L_classified_tree 

# make an NP tree
NP_classified_tree <- 
  ggtree(NP_midpoint, ladderize = F) %<+% S_metadata +
  geom_tippoint(aes(fill = classified), 
                size=1.5, color="black", stroke=0.2, shape=21) +
  scale_fill_manual(values=c("white", "red")) +
  shared_tree_theme()  +
  geom_treescale(linesize=0.5, x = 0, y = -5, fontsize=3, offset=2) +
  ggtitle("NP")

NP_classified_tree 

# make a combined figure
classified_unclassified_fig <- L_classified_tree + NP_classified_tree
classified_unclassified_fig

# save as an RDS and a pdf
saveRDS(classified_unclassified_fig, file="../figures/classified_unclassified_figure.RDS")
ggsave(classified_unclassified_fig, file="../figures/extras/fig_x_classified_unclassified_figure.pdf", 
       units="in", width=7.5, height=5)

# make supplemental figure versions with tip labels
supplemental_L_classified_unclassified_fig <- 
  L_classified_tree +
  geom_tiplab(aes(label=acc_org), size=1) +   
  coord_cartesian(clip="off") 
  
supplemental_S_classified_unclassified_fig <- 
  NP_classified_tree +
  geom_tiplab(aes(label=acc_org), size=1)  +
  coord_cartesian(clip="off")

# save as RDS and pdf
saveRDS(supplemental_L_classified_unclassified_fig, "../figures/extras/Supplemental_figure_X_L_classified_unclassified_tree.RDS")
saveRDS(supplemental_S_classified_unclassified_fig, "../figures/extras/Supplemental_figure_X_S_classified_unclassified_tree.RDS")
ggsave("../figures/extras/Supplemental_figure_X_L_classified_unclassified_tree.pdf", 
       supplemental_L_classified_unclassified_fig,
       units="in", width=7.5, height=10)
ggsave( "../figures/extras/Supplemental_figure_X_S_classified_unclassified_tree.pdf", 
        supplemental_S_classified_unclassified_fig,
        units="in", width=7.5, height=10)
  

# calculate and save some values for the paper
classified_unclassified_summary <- metadata %>% 

# Aurora borealis virus 3, taxid 2447924
  
# accession_taxon_map <- accession_taxon_map %>% mutate(ugv_label = if_else(taxid == 1672385, acc_org, ""))



# get node 
MRCA()
# University of Helsinki virus, taxid 1382279
ggtree(L_midpoint, ladderize = F) %<+% accession_taxon_map +
  geom_tiplab(aes(label=if_else(taxid==1382279, acc_org, "")), size=1) 
  # geom_tippoint(aes(fill = unclassified), 
                # size=1.5, color="black", stroke=0.2, shape=21) 
  
# Keijut pohjoismaissa virus 1, taxid 2447927
ggtree(L_midpoint, ladderize = F) %<+% accession_taxon_map +
  geom_tiplab(aes(label=if_else(taxid==2447927, acc_org, "")), size=1)
  

L_t <- L_tree_with_taxon <- L_tree %<+% accession_taxon_map 

full_join(L_tree, accession_taxon_map)
