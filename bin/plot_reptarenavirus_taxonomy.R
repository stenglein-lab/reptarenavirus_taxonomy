library(tidyverse)
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(seqinr)
# library(ggsci)
library(RColorBrewer)
# library(ggExtra)
# library(ragg)

# output directory
output_dir = "../paper/figures/"

# read in trees
L_tree  <- read.newick("../make_alignments_and_trees/results/trees/L_CDS.mafft.fasta.treefile")
NP_tree <- read.newick("../make_alignments_and_trees/results/trees/NP_CDS.mafft.fasta.treefile")

# midpoint root trees
L_midpoint  <- midpoint_root(L_tree)
NP_midpoint <- midpoint_root(NP_tree)

# Read in assignment of accession (tip label) -> NCBI-defined taxon
metadata <- readRDS("../metadata/reptarenavirus_accession_metadata.RDS")

# how many sequences for of each species / segment?
metadata <- metadata %>% group_by(species_taxid, segment) %>% mutate(n_representatives=n()) %>% ungroup()

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

###############################
# Functions to make trees
###############################
###############################

# a shared graphical theme for plots
shared_tree_theme <- function() {
  theme_tree(base_size = 11) +
    theme(text=element_text(family="Helvetica"))
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


###################################
# classified / unclassified trees
###################################
# L tree
L_classified_tree <- plot_tree(L_midpoint, 0, 1, tree_title="L")
L_classified_tree 

# NP tree
NP_classified_tree <- plot_tree(NP_midpoint, 0, 1, tree_title="NP")  + theme(legend.position = "right")
NP_classified_tree 

# make a combined figure
classified_unclassified_fig <- L_classified_tree + NP_classified_tree
classified_unclassified_fig

# A function to save PNGs using base functions
# savePNG <- function(plot, file, width=7.5, height=10, units="in", res=300) {
#   png(file, width = width, height = height, units = units, res = res)
#   print(plot)
#   dev.off()
# }

# save as an RDS and a pdf
saveRDS(classified_unclassified_fig, file=paste0(output_dir,"classified_unclassified_figure.RDS"))
ggsave(classified_unclassified_fig, file=paste0(output_dir, "classified_unclassified_figure.pdf"), 
       units="in", width=7.5, height=5)
ggsave(classified_unclassified_fig, file=paste0(output_dir, "classified_unclassified_figure.png"), 
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
  
###################################

###############################
# Tanglegrams
###############################

# Make association for labels

# current tree tips
L_tips  <- tibble(accession = L_midpoint$tip.label)
NP_tips <- tibble(accession = NP_midpoint$tip.label)

# get rid of unclassified 
acc_taxid <- metadata %>% filter(classified)  %>% select(accession, species_taxid)

L_tips  <- left_join(L_tips, acc_taxid) %>% filter(!is.na(species_taxid))
NP_tips <- left_join(NP_tips, acc_taxid) %>% filter(!is.na(species_taxid))

# how many links should their be in the tree?
L_per_spp <- L_tips %>% group_by(species_taxid) %>% summarize(nL=n())
NP_per_spp <- NP_tips %>% group_by(species_taxid) %>% summarize(nNP=n())
n_links <- left_join(L_per_spp, NP_per_spp) %>% mutate(per_spp_links = nL * nNP)

# make the association matrix
# The association matrix for phytools::cophylo()
# from: https://rdrr.io/cran/phytools/man/cophylo.html
# matrix containing the tip labels in tr1 to match to the tip labels in tr2. 
# Note that not all labels in either tree need to be included; and, furthermore, 
# one label in tr1 can be matched with more than one label in tr2, or vice versa.
L_NP_assoc <- inner_join(L_tips, NP_tips, by=join_by(species_taxid), relationship = "many-to-many") 

L_NP_assoc_mat <- L_NP_assoc %>% select(-species_taxid)
L_NP_assoc_mat <- as.matrix(L_NP_assoc_mat)
colnames(L_NP_assoc_mat) <- c("L_acc", "NP_acc")

# check the join did what we wanted it too
stopifnot(sum(n_links$per_spp_links) == nrow(L_NP_assoc_mat))
  
# make cophylogeny (this makes the object but doesn't plot it yet)
cophy <- cophylo(L_midpoint, NP_midpoint, assoc = L_NP_assoc_mat)

# make color schemes for different species   
colors = brewer.pal(5, "Dark2")

spp_taxids <- filter(metadata, classified) %>% group_by(species_taxid) %>% summarize() %>% pull(species_taxid)
n_taxid = length(spp_taxids)

taxid_color_map <- tibble(species_taxid = spp_taxids, color = rev(colors))

# color tanglegram lines by species classification for classified seqs
# or grey for unclassified seqs
# the link.col parameter in plot.cophylo is a vector whose *length and order* matches the assoc matrix
link_colors <- left_join(L_NP_assoc, taxid_color_map)
link_colors_vec <- pull(link_colors, color)

test_link_colors = rep("lightgrey", nrow(L_NP_assoc_mat))
test_link_colors[length(test_link_colors)] <- "coral3"


make_classified_tanglegram <- function(){
  # plot tanglegram
  plot.cophylo(cophy,
               link.type = "curved",
               link.lwd=1,
               link.lty="solid",
               link.col=link_colors_vec,
               fsize=0.2,
               pts=F,
               scale.bar = c(0.5, 0.5))
  
  # this doesn't plot node support values, but it could.
  # see: http://blog.phytools.org/2015/10/node-edge-tip-labels-for-plotted.html
}

# save as PDF
pdf(file = "../paper/figures/classified_tanglegram.pdf", width=8.5, height=11)
make_classified_tanglegram()
dev.off()

# png(file = "../paper/figures/L_NP_taxid_tanglegram.png", width=8.5, height=11)
# make a color-map of Stenglein vs Hepojoki snakes
color_by_paper <- tibble(snake_number = multiply_infected_snakes,
                         color = case_when(
                           is.na(snake_number)                   ~ NA,
                           str_detect(snake_number, "Hepojoki")  ~ "coral2",
                           .default = "slateblue" )) 

# png(file = "../paper/figures/classified_tanglegram.png", res=300, width = 2250, height=2850)
# make_classified_tanglegram()
# dev.off()


###############################

########################################################
# Make tanglegrams showing multiply-infected snakes
########################################################

# current tree tips
L_tips  <- tibble(accession = L_midpoint$tip.label)
NP_tips <- tibble(accession = NP_midpoint$tip.label)

# map of acccessions -> snake #s
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
ggsave(paste0(output_dir, "segments_per_snake.pdf"), units="in", width=7, height=6)
  
# these snakes only have 2 segments
singly_infected_snakes <- filter(seqs_per_snake, n_seqs == 2) %>% pull(snake_number)

# these snakes only > 2 segments
multiply_infected_snakes <-  filter(seqs_per_snake, n_seqs > 2 & !is.na(snake_number)) %>% pull(snake_number)

# make a cophylogeny for multiple infection plots: this makes the object but doesn't plot it yet
full_cophy <- cophylo(L_midpoint, NP_midpoint, assoc = L_NP_assoc_mat)

# make a tanglegram showing links between co-infecting sequences in one or more snake
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
    pdf(file = paste0(output_dir, filename_prefix, ".pdf"), width=8.5, height=11)
    plot_multiply_infected_tanglegram()
    dev.off()
}

# tanglegrams showing linked sequences from all multiply-infected snakes

# make a color map with different colors for snakes from different papers
color_by_paper <- tibble(snake_number = multiply_infected_snakes,
                         color = case_when(
                           is.na(snake_number)                   ~ NA,
                           str_detect(snake_number, "Hepojoki")  ~ "coral2",
                           .default = "slateblue" )) 

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

make_multiply_infected_tanglegram(snakes_to_highlight = color_by_paper$snake_number,
                                  filename_prefix     = "multiply_infected_snakes_color_by_paper",
                                  highlight_colors    = color_by_paper$color)

make_multiply_infected_tanglegram(snakes_to_highlight = color_by_paper$snake_number,
                                  filename_prefix     = "multiply_infected_snakes_color_by_snake",
                                  highlight_colors    = color_by_paper$color)

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



########################################################


  

