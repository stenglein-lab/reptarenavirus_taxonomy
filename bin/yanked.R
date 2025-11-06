# yanked code from reptarenavirus taxonomy analysis

get_accessions <- function(taxon_map, this_taxid) {
  filter(taxon_map, taxid == this_taxid) %>% pull(accession)
}
ABV3_acc <- get_accessions(L_acc_taxon_map, 1382279)
is_monophyletic(L_tree, ABV3_acc)

ABV3_acc <- get_accessions(L_acc_taxon_map, 1223561)
is_monophyletic(L_tree, ABV3_acc)

# quantify monophyleticness
# is 
is_monophyletic <- function(tree, accessions) {
  # get the MRCA
  this_mrca <- MRCA(tree, accessions) 
  
  # get the descendants of the MRCA: does it include other accessions?
  mrca_descendant_tips <- offspring(tree, this_mrca, type = "tips")
  
  # check if they are the same and return TRUE/FALSE
  identical(mrca_descendant_tips, accessions)
}



is_monophyletic(L_tree)
