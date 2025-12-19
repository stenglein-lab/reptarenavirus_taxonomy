library(tidyverse)
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(seqinr)

# read in trees
L_tree  <- read.newick("../make_alignments_and_trees//results/trees/L_CDS.mafft.fasta.treefile")
NP_tree <- read.newick("../make_alignments_and_trees//results/trees/NP_CDS.mafft.fasta.treefile")

# midpoint root trees
L_midpoint  <- midpoint_root(L_tree)
NP_midpoint <- midpoint_root(NP_tree)

# Read in assignment of accession (tip label) -> NCBI-defined taxon
metadata <- readRDS("../metadata/reptarenavirus_accession_metadata.RDS")

# how many sequences for of each species / segment?
metadata <- metadata %>% group_by(species_taxid, segment) %>% mutate(n_representatives=n())

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

saveRDS(classified_values, "../RDS/classified_values.RDS")

L_metadata <- filter(metadata, segment == "L")
S_metadata <- filter(metadata, segment == "S")

# a shared graphical theme for plots
shared_tree_theme <- function() {
  theme_tree(base_size = 11) +
    theme()
}

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

# L tree
L_classified_tree <- plot_tree(L_midpoint, 0, 1, tree_title="L")
L_classified_tree 

# NP tree
NP_classified_tree <- plot_tree(NP_midpoint, 0, 1, tree_title="NP")  + theme(legend.position = "right")
NP_classified_tree 

# make a combined figure
classified_unclassified_fig <- L_classified_tree + NP_classified_tree
classified_unclassified_fig

# save as an RDS and a pdf
saveRDS(classified_unclassified_fig, file="../paper/figures/classified_unclassified_figure.RDS")
ggsave(classified_unclassified_fig, file="../paper/figures/classified_unclassified_figure.pdf", 
       units="in", width=7.5, height=5)
ggsave(classified_unclassified_fig, file="../paper/figures/classified_unclassified_figure.png", 
       units="in", width=7.5, height=5)

# make supplemental figure versions with tip labels
supplemental_L_classified_unclassified_fig <- 
  plot_tree(L_midpoint, 0.62, 1.2, 3, tree_title="L") + 
  geom_tiplab(aes(label=acc_org), size=1.5) 

supplemental_L_classified_unclassified_fig 
  
supplemental_S_classified_unclassified_fig <- 
  plot_tree(NP_midpoint, 0.75, 1.3, 3, tree_title="NP") + 
  geom_tiplab(aes(label=acc_org), size=2)  

supplemental_S_classified_unclassified_fig 

# save as RDS and pdf
# saveRDS(supplemental_L_classified_unclassified_fig, "../figures/extras/Supplemental_figure_X_L_classified_unclassified_tree.RDS")
# saveRDS(supplemental_S_classified_unclassified_fig, "../figures/extras/Supplemental_figure_X_S_classified_unclassified_tree.RDS")
ggsave("../paper/figures/supplemental/Supplemental_figure_X_L_classified_unclassified_tree.pdf", 
       supplemental_L_classified_unclassified_fig,
       units="in", width=7.5, height=10)
ggsave( "../paper/figures/supplemental/Supplemental_figure_X_S_classified_unclassified_tree.pdf", 
        supplemental_S_classified_unclassified_fig,
        units="in", width=7.5, height=9)
ggsave("../paper/figures/supplemental/Supplemental_figure_X_L_classified_unclassified_tree.png", 
       supplemental_L_classified_unclassified_fig,
       units="in", width=7.5, height=10)
ggsave( "../paper/figures/supplemental/Supplemental_figure_X_S_classified_unclassified_tree.png", 
        supplemental_S_classified_unclassified_fig,
        units="in", width=7.5, height=9)
  

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
  


