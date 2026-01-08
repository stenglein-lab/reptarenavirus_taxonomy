library(tidyverse)
library(readxl)
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(seqinr)
library(RColorBrewer)
library(MSA2dist)
library(patchwork)


#
# This script contains code used for analyses related to the paper:
# "A proposal for a simplified reptarenavirus taxonomy based on reassorment compatibility"
#
# Mark Stenglein 11/6/2025
#

# define output directories
output_dir             = "../paper/"
rds_output_dir         = paste0(output_dir, "RDS/")
figure_output_dir      = paste0(output_dir, "figures/")
supp_figure_output_dir = paste0(figure_output_dir, "supplemental/")
table_output_dir       = paste0(output_dir, "tables/")


# Handle sequence metadata  -------------------------------------------------------------------

# combine reptarenavirus metadata from difference sources:
#
# 1. from genbank files of nucleotide sequences: accession, organism name, sequence taxid, segment
# 2. from E-Utils: accession, species name, species taxid
# 3. A table corresponding to table 1 from Hepojoki et al (2015), PMID: 26041290, re: multiply infected snakes
# 4. The master species list (MLS) from the ICTV website
#

genbank_derived_map <- read.delim("../metadata/genbank_metadata.txt")
eutils_derived_map  <- read.delim("../metadata/accession_species_taxid_map.txt")

# read in spreadsheet containing table 1 from Hepojoki et al (2015), PMID: 26041290, re: multiply infected snakes
Hepojoki_table_1 <- read_excel("../metadata/Hepojoki_2015_table_1.xlsx")
Hepojoki_snakes <- select(Hepojoki_table_1, accession, Hepojoki_snake_number) 

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
write.table(reptarenavirus_accession_metadata, 
            file = paste0(table_output_dir, "reptarenavirus_accession_metadata.txt"),
            sep="\t", row.names = F, quote=F)

metadata <- reptarenavirus_accession_metadata


# Reptarenavirus species information from the ICTV MSL ----------------------------------------

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

saveRDS(reptarenavirus_species, paste0(rds_output_dir, "reptarenavirus_species.RDS"))


# Classified / unclassified stats and figures -------------------------------------------------

# read in trees
L_tree  <- read.newick("../make_alignments_and_trees/results/trees/L_CDS.mafft.fasta.treefile")
NP_tree <- read.newick("../make_alignments_and_trees/results/trees/NP_CDS.mafft.fasta.treefile")

# midpoint root trees
L_midpoint  <- midpoint_root(L_tree)
NP_midpoint <- midpoint_root(NP_tree)

# how many classified / unclassified?
classified_counts <- metadata %>% group_by(segment, classified) %>% summarize(n=n())

# how many S & L
total_sequence_counts <- classified_counts %>% group_by(segment) %>% summarize(total_seq = sum(n))
L_sequences = total_sequence_counts %>% filter(segment == "L") %>% pull(total_seq)
S_sequences = total_sequence_counts %>% filter(segment == "S") %>% pull(total_seq)

# how many classified 
classified_counts_total = classified_counts %>% group_by(classified) %>% summarize(n = sum(n))
classified_total        = classified_counts_total %>% filter(classified) %>% pull(n)
unclassified_total      = classified_counts_total %>% filter(!classified) %>% pull(n)
fraction_classified     = classified_total   / (classified_total + unclassified_total)
fraction_unclassified   = unclassified_total / (classified_total + unclassified_total)
 

classified_L     = classified_counts %>% filter(segment == "L" & classified) %>% pull(n) 
classified_S     = classified_counts %>% filter(segment == "S" & classified) %>% pull(n) 
fraction_classified_L = classified_L / L_sequences
fraction_classified_S = classified_S / S_sequences

# make a list of data values to be included in paper
classified_values <- list(
  L_sequences = L_sequences,
  S_sequences = S_sequences,
  classified_L = classified_L,
  classified_S = classified_S,
  fraction_classified_L = fraction_classified_L,
  fraction_classified_S = fraction_classified_S,
  fraction_classified   =  fraction_classified,
  fraction_unclassified =  fraction_unclassified
)

# these values will be used in paper
saveRDS(classified_values, paste0(rds_output_dir, "classified_values.RDS"))

L_metadata <- filter(metadata, segment == "L")
S_metadata <- filter(metadata, segment == "S")

# Functions for making regular trees (non-tanglegrams) 

# a shared graphical theme 
shared_tree_theme <- function() {
  theme_tree(base_size = 11) +
    theme(text=element_text(family="Helvetica"))
}

# plot one tree
plot_tree <- function(tree_to_plot, 
                      node_depth_label_cutoff = 0.97,
                      x_axis_scalar = 1.25,
                      tip_font_size = 3,
                      tree_title = NA) {
  
  debug <- 0
  if (debug) {
    tree_to_plot <- L_midpoint
    node_depth_label_cutoff = 0.97
  }
  
  # find max root to tip distance to create reasonable X limits for plotting
  root_to_tip_distances <- node.depth.edgelength(tree_to_plot)
  tip_distances <- root_to_tip_distances[1:Ntip(tree_to_plot)]
  max_root_to_tip_distance <- max(tip_distances)
  
  # root edge length, based on total tree length
  root_edge_length = max_root_to_tip_distance * 0.025
  
  # Programmatically decide which support values to show:
  # only show support values for internal nodes not near the tips 
  # (based on fraction of depth)
  # 
  # alternative schemes to decide which nodes to label could include:
  # - labeling nodes with a certain number of tip descendants
  # - labeling nodes that have support values above a certain cutoff
  #
  
  # if this tree has internal node labels
  if ( length(tree_to_plot$node.label) > 0 ) {
    
    # depths of internal nodes
    internal_node_depths <-  root_to_tip_distances[-(1:Ntip(tree_to_plot))]
    # store in data frame and keep track of existing labels (support values for internal nodes)
    node_labels <- tibble(node_depth = internal_node_depths,
                          node_labels = tree_to_plot$node.label)
    # calculate a fractional depth: how far along the total tree length is this internal node?
    node_labels <- mutate(node_labels, fractional_depth = node_depth / max_root_to_tip_distance)
    # create new labels for nodes sufficiently far away from tips
    # also don't label root node support value
    node_labels <- mutate(node_labels, new_label = 
                            if_else(fractional_depth > 0 & 
                                      fractional_depth <= node_depth_label_cutoff, 
                                    node_labels, ""))
    
    # relabel internal nodes 
    tree_to_plot$node.label <- pull(node_labels, new_label)
  }
  
  # create a metadata object that encodes whether tips are from one of our new sequences
  # parse out accession from tip labels: characters after final underscore
  our_tips <- tibble(tip_label = tree_to_plot$tip.label)
  
  # make the ggtree plot
  tree_plot <- 
    ggtree(tree_to_plot, ladderize=FALSE) %<+% 
    # merge in metadata for coloring tips
    metadata +
    # tip labels, colored by whether our seq or not
    geom_tippoint(aes(fill = classified), 
                  size=1.5, color="black", stroke=0.2, shape=21) +
    # node support values
    geom_nodelab(size=2.5, nudge_x = 0, color="grey50") +
    # scale bar
    geom_treescale(linesize=0.5, x = 0, y = -15, fontsize=3, offset=2) +
    # root edge
    geom_rootedge(root_edge_length) +
    # color scale for our seqs vs genbank 
    scale_fill_manual(values=c("white", "red")) +
    theme_tree() +
    # make space for root edge and labels
    xlim(-root_edge_length, max_root_to_tip_distance * x_axis_scalar) +
    # don't show legend
    theme(legend.position = "none") +
    shared_tree_theme() 
  
  # add optional titlte
  if (!is.na(tree_title))  {
    tree_plot <- tree_plot + ggtitle(tree_title)
  }
  
  tree_plot
}


# create classified / unclassified tree figures

# L tree
L_classified_tree <- plot_tree(L_midpoint, 0, 1, tree_title="L")
L_classified_tree 

# NP tree
NP_classified_tree <- plot_tree(NP_midpoint, 0, 1, tree_title="NP")  + theme(legend.position = "right")
NP_classified_tree 

# make a combined figure
classified_unclassified_fig <- L_classified_tree + NP_classified_tree
classified_unclassified_fig

# save as pdf and png 
ggsave(classified_unclassified_fig, 
       file=paste0(figure_output_dir, "classified_unclassified_figure.pdf"), 
       units="in", width=7.5, height=5)
ggsave(classified_unclassified_fig, 
       file=paste0(figure_output_dir, "classified_unclassified_figure.png"), 
       units="in", width=7.5, height=5, device=grDevices::png)

# make supplemental figure versions with tip labels
supplemental_L_classified_unclassified_fig <- 
  plot_tree(L_midpoint, 0.62, 1.2, 3, tree_title="L") + 
  geom_tiplab(aes(label=acc_org), size=1.5) 

supplemental_L_classified_unclassified_fig 
  
supplemental_S_classified_unclassified_fig <- 
  plot_tree(NP_midpoint, 0.75, 1.3, 3, tree_title="NP") + 
  geom_tiplab(aes(label=acc_org), size=2)  

supplemental_S_classified_unclassified_fig 

# save as png and pdf
ggsave(
  paste0(supp_figure_output_dir, "Supplemental_figure_X_L_classified_unclassified_tree.pdf"), 
  supplemental_L_classified_unclassified_fig,
  units="in", width=7.5, height=10)
ggsave( 
  paste0(supp_figure_output_dir, "Supplemental_figure_X_S_classified_unclassified_tree.pdf"), 
  supplemental_S_classified_unclassified_fig,
  units="in", width=7.5, height=9)
ggsave( 
  paste0(supp_figure_output_dir, "Supplemental_figure_X_L_classified_unclassified_tree.png"), 
  supplemental_L_classified_unclassified_fig,
  units="in", width=7.5, height=10, device=grDevices::png)
ggsave( 
  paste0(supp_figure_output_dir, "Supplemental_figure_X_S_classified_unclassified_tree.png"), 
  supplemental_S_classified_unclassified_fig,
  units="in", width=7.5, height=9, device=grDevices::png)
  

# Tanglegrams ---------------------------------------------------------------------------------------

# Make tanglegrams showing multiply-infected snakes

# current tree tips
L_tips  <- tibble(accession = L_midpoint$tip.label)
NP_tips <- tibble(accession = NP_midpoint$tip.label)

# map of acccessions -> snake #s
# snake #s are defined in Stenglein et al (2015) and Hepojoki et al (2015) and are part
# of metadata processed above
acc_snake <- metadata %>% select(accession, snake_number)
L_tips  <- left_join(L_tips, acc_snake) %>% filter(!is.na(snake_number))
NP_tips <- left_join(NP_tips, acc_snake) %>% filter(!is.na(snake_number))

# how many links should their be in the tree?
L_per_spp <- L_tips %>% group_by(snake_number) %>% summarize(nL=n())
NP_per_spp <- NP_tips %>% group_by(snake_number) %>% summarize(nNP=n())
n_links <- left_join(L_per_spp, NP_per_spp) %>% mutate(per_spp_links = nL * nNP)

# make the association matrix
# The association matrix for phytools::cophylo()
# from: https://rdrr.io/cran/phytools/man/cophylo.html
# matrix containing the tip labels in tr1 to match to the tip labels in tr2. 
# Note that not all labels in either tree need to be included; and, furthermore, 
# one label in tr1 can be matched with more than one label in tr2, or vice versa.
L_NP_assoc <- inner_join(L_tips, NP_tips, by=join_by(snake_number), relationship = "many-to-many") 

# full cophylogeny
L_NP_assoc_mat <- L_NP_assoc %>% select(-snake_number)
L_NP_assoc_mat <- as.matrix(L_NP_assoc_mat)
colnames(L_NP_assoc_mat) <- c("L_acc", "NP_acc")

# check the join did what we wanted it too
stopifnot(sum(n_links$per_spp_links, na.rm = T) == nrow(L_NP_assoc_mat))

# how many segments per snake?
seqs_per_snake_per_segment <- metadata %>% group_by(snake_number, segment) %>% summarize(n_seqs=n())
seqs_per_snake             <- metadata %>% group_by(snake_number) %>% summarize(n_seqs=n())

# make a jitter scatter plot of segments per snake
seqs_per_snake_per_segment_wide <- 
  seqs_per_snake_per_segment %>% pivot_wider(names_from = segment, values_from = n_seqs) %>%
  mutate(paper = if_else(str_detect(snake_number, "Hepojoki"), "Hepojoki_2015", "Stenglein_2015"))

ggplot(filter(seqs_per_snake_per_segment_wide, !is.na(snake_number))) +
  geom_jitter(aes(x=S, y = L, fill=paper), 
              shape=21, size=3, color="black", 
              stroke=0.25, width=0.05, height=0.05,
              alpha=0.75) +
  scale_fill_manual(values=c("coral2", "slateblue")) +
  theme_bw() + 
  xlab("Number of S segments per snake") +
  ylab("Number of L segments per snake") 

ggsave(paste0(figure_output_dir, "segments_per_snake.pdf"), units="in", width=7, height=6)
  
# these snakes only have 2 segments
singly_infected_snakes <- filter(seqs_per_snake, n_seqs == 2) %>% pull(snake_number)

# these snakes have > 2 segments, i.e. are multiply infected
multiply_infected_snakes <-  filter(seqs_per_snake, n_seqs > 2 & !is.na(snake_number)) %>% pull(snake_number)

# make a cophylogeny for multiple infection plots 
# this makes the object but doesn't plot it yet
full_cophy <- cophylo(L_midpoint, NP_midpoint, assoc = L_NP_assoc_mat)

# this function makes a tanglegram showing links between co-infecting sequences in one or more snake
make_multiply_infected_tanglegram <- 
  function(snakes_to_highlight   = multiply_infected_snakes, 
           cophy                 = full_cophy,
           highlight_colors      = rep ("grey90",length(snakes_to_highlight)),
           filename_prefix       = "multiply_infected_tanglegram",
           label_tips            = T )  {
    
    debug <- 0
    if (debug){
      snakes_to_highlight   = color_by_snake$snake_number
      cophy                 = full_cophy
      highlight_colors      = color_by_snake$color
      filename_prefix       = "multiply_infected_tanglegram_debug"
    }
    
    # check that enough colors specified
    stopifnot(length(highlight_colors) == length(snakes_to_highlight))
    
    # this will map tip labels (accession)
    snake_num_color_map <- tibble(snake_number = snakes_to_highlight,
                                  color = highlight_colors)
    
    link_color_map <- left_join(L_NP_assoc, snake_num_color_map)
    
    # the link.col parameter in plot.cophylo is a vector whose *length and order* matches the assoc matrix
    link_colors_vec <- pull(link_color_map, color)
    
    plot_multiply_infected_tanglegram <- function(label_tips = T){
      plot.cophylo(cophy,
                   link.type = "curved",
                   link.lwd  = 1,
                   link.lty  = "solid",
                   link.col  = link_colors_vec,
                   fsize     = c(0.25, 0.3),
                   pts       = F,
                   scale.bar = c(0.5, 0.5))
    }
    
    # plot and save as PDF
    pdf(file = paste0(figure_output_dir, filename_prefix, ".pdf"), width=8.5, height=11)
    plot_multiply_infected_tanglegram()
    dev.off()
}

# tanglegrams showing linked sequences from all multiply-infected snakes

# make a color map with different colors for snakes from different papers
# this is just singly-infected snakes
color_by_paper_singly <- tibble(snake_number = singly_infected_snakes,
                         color = case_when(
                           is.na(snake_number)                   ~ NA,
                           str_detect(snake_number, "Hepojoki")  ~ "coral2",
                           .default = "slateblue" )) 

# make a color map with a unique color for each snake
big_palette <- brewer.pal(9, "Set1") 
big_palette <- colorRampPalette(big_palette)(length(multiply_infected_snakes))
# pie(rep(1, length(big_palette)), col = big_palette , main="") 

# some colorblind friendly color palettes 
# from: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# pie(rep(1, 8), col = colorBlindBlack8)
colorBlind2  <- c("#0072B2", "#D55E00")
colorBlind3  <- c("#000000", "#0072B2", "#D55E00")

color_by_snake <- tibble(snake_number = multiply_infected_snakes,
                         color = big_palette)

# all singly-infected snakes
make_multiply_infected_tanglegram(snakes_to_highlight = color_by_paper_singly$snake_number,
                                  filename_prefix     = "singly_infected_snakes", 
                                  highlight_colors    = color_by_paper_singly$color)

# exporting this figure as PDF and will manually add annotation (snake #s, tree labels) in Affinity Designer
make_multiply_infected_tanglegram(snakes_to_highlight = c("Stenglein_snake_26", "Stenglein_snake_33", "Hepojoki_2015_snake_8"), 
                                  filename_prefix     = "snake_26_33_Hepojoki_snake_8_multiple_infection", 
                                  highlight_colors    = colorBlind3)

# exporting this figure as PDF and will manually add annotation (snake #s, tree labels) in Affinity Designer
# these are all the snakes from Stenglein 2015 and Hepojoki 2015
all_connected_snakes <-  filter(seqs_per_snake, !is.na(snake_number)) %>% pull(snake_number)
make_multiply_infected_tanglegram(snakes_to_highlight = all_connected_snakes,
                                  filename_prefix     = "all_connected_snakes", 
                                  highlight_colors    = rep("grey90", length(all_connected_snakes)))



# read in alignments
# nt alignments
L_MSA   <- read.alignment("../make_alignments_and_trees/results/alignments/L_CDS.mafft.fasta", format = "fasta")
NP_MSA  <- read.alignment("../make_alignments_and_trees/results/alignments/NP_CDS.mafft.fasta", format = "fasta")

# read in protein alignments
L_prot_MSA   <- read.alignment("../make_alignments_and_trees/results/alignments/L_prot.mafft.fasta", format = "fasta")
NP_prot_MSA  <- read.alignment("../make_alignments_and_trees/results/alignments/NP_prot.mafft.fasta", format = "fasta")

# convert to pairwise distances using dist.alignment

# note that these are just basic pairwise distances, which under-report true evolutionary distance
# especially for distantly related nucleotide sequences, because of saturating mutations

# nt sequence alignments
# square each element of the matrix because dist.alignment returns square root of distance
L_dist   <- as.matrix(dist.alignment(L_MSA, matrix="identity"))
L_dist   <- L_dist ^ 2 
NP_dist  <- as.matrix(dist.alignment(NP_MSA, matrix="identity"))
NP_dist  <- NP_dist ^ 2 

# protein sequence alignments
# square each element of the matrix because dist.alignment returns square root of distance
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

# values to be used in paper
saveRDS(closest_MN567045_relative, paste0(rds_output_dir, "closest_MN567045_relative.RDS"))
saveRDS(closest_MW091469_relative, paste0(rds_output_dir, "closest_MW091469_relative.RDS"))

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

# these RDS objects will be imported and use to populate values in paper
saveRDS(min_example_S_dist, paste0(rds_output_dir, "min_example_S_dist.RDS"))
saveRDS(min_example_S_dist, paste0(rds_output_dir, "max_example_S_dist.RDS"))
saveRDS(min_example_L_dist, paste0(rds_output_dir, "min_example_L_dist.RDS"))
saveRDS(max_example_L_dist, paste0(rds_output_dir, "max_example_L_dist.RDS"))

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

# do same for L distances
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

# create histograms of pairwise distances between sequences from within and between species 
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

# make a 2-paneled figure of pairwise distance histograms
dist_p <- 
  L_dist_p + NP_dist_p + plot_layout(ncol = 2) + plot_annotation(tag_levels = 'A')
dist_p

ggsave(dist_p, file=paste0(figure_output_dir, "pairwise_distance_figure.pdf"), 
       units="in", width=7.5, height=5)
ggsave(dist_p, file=paste0(figure_output_dir, "pairwise_distance_figure.png"), 
       units="in", width=7.5, height=5, device=grDevices::png)


# create a supplemental table listing the proposed new species 
# assignment for all available reptarenavirus sequences

new_species_assignments <- tibble(accession = metadata$accession)
new_species_assignments <- new_species_assignments %>% mutate(
  proposed_species = case_when(
    accession %in% R_californae_accessions ~ "Reptarenavirus californae",
    .default = "Reptarenavirus giessenae"
  )
)

# write out the table
write.table(arrange(new_species_assignments, proposed_species) , 
            file=paste0(table_output_dir, "Supplemental_table_accession_proposed_species_map.txt"),
            quote=F,
            sep="\t",
            row.names = F)
