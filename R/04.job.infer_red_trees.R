################################################################################
#Estimate clade counts from marker gene trees: Bacteria
#
#Script does the following things:
# - load tree and tree data
# - iterate through alternative outgroup rootings and calculate RED values for each permutation
# - cut full tree at pre-defined RED values (from GTDB) and count predicted number of clades as well as their size
#
#Usage:
#Rscript cut.clades.Bacteria.R [marker_id]
################################################################################
#Load libraries
library(tidyverse, quietly = T)
library(castor, quietly = T)
library(ape, quietly = T)
library(phangorn, quietly = T)
library(phytools, quietly = T)
################################################################################


################################################################################
#Get command line argument
args = commandArgs(trailingOnly=TRUE)

marker_id <- args[1]

curr.file.tree_data <- paste0(marker_id, ".tree_data.Rdata")

if (!file.exists(curr.file.tree_data)) {warning(paste("WARNING: file", curr.file.tree_data, "does not exist. Quitting.")); q()}

load(curr.file.tree_data)

#Set folder where stuff will be stored
#folder.store <- [insert_output_folder]

#Load sample metadata
load("data.sample.rarefy.Rdata")

#Load gene clustering data
writeLines(paste(date(), marker_id, "=> loading gene clustering data."))
load(paste0(marker_id, ".spire_gtdb.clustered.965.Rdata"))

#Set reference clade counts
#=> these are taken from https://gtdb.ecogenomic.org/stats/r220
phylumGTDB = 175 
classGTDB = 538
orderGTDB = 1840
familyGTDB = 4870
genusGTDB = 23112
################################################################################


################################################################################
#Get (temporary) midpoint rooting
writeLines(paste(date(), marker_id, "=> rooting the tree."))
rootedByMidPoint <- midpoint(tree.raw)

#Get node labels
nodeLabsRaw <- as.data.frame(rootedByMidPoint$node.label)
rootedByMidPointLabeled <- makeNodeLabel(rootedByMidPoint)
nodeLabsNew <- as.data.frame(rootedByMidPointLabeled$node.label)
################################################################################


################################################################################
#Get current relevant phyla
use.phyla <- dat.tree.tips %>%
  filter(!is.na(phylum)) %>%
  group_by(phylum) %>%
  tally %>%
  filter(n >= 10) %>%
  pull(phylum)

#use.phyla
#Get marker-gene-specific under-detection of GTDB r220 clades
convert.clade_counts <- dat.tree.tips %>%
  filter(source != "spire.unbinned") %>%
  select(raw.tip.label, phylum, class, order, family, genus) %>%
  pivot_longer(phylum:genus) %>%
  filter(!is.na(value)) %>%
  group_by(name, value) %>%
  tally %>%
  group_by(name) %>%
  tally %>%
  pivot_wider(values_from = n) %>%
  #These are numbers for archaea from the GTDB website...
  mutate(
    phylum = phylum / phylumGTDB,
    class = class / classGTDB,
    order = order / orderGTDB,
    family = family / familyGTDB,
    genus = genus / genusGTDB
  )

#Preallocate results collector
collect.red <- tibble()

#Iterate through phyla, re-root tree, calculate and store RED
collect.red <- tibble()
writeLines(paste(date(), marker_id, "=> calculating RED."))
for (phy in use.phyla) {
  #Get current tip set for outgroup
  curr.tips <- dat.tree.tips %>%
    filter(phylum == phy) %>%
    pull(raw.tip.label)
  
  #Get MRCA of outgroup (=> this will become the root)
  curr.mrca <- get_mrca_of_set(rootedByMidPoint, curr.tips)
  
  #Re-root tree at MRCA of current phylum's tips
  tree.rerooted.mrca <- root_at_node(rootedByMidPointLabeled, curr.mrca - Ntip(tree.raw))
  
  #Get RED values
  curr.red <- tibble(
    node = tree.rerooted.mrca$node.label,
    red = get_reds(tree.rerooted.mrca)
  )
  
  #Store RED values
  collect.red <- bind_rows(
    collect.red,
    curr.red
  )
}
writeLines(paste(date(), marker_id, "=> RED calculations done."))

#Get average RED values per node (except when they're the root)
avg.red <- collect.red %>%
  filter(red > 0) %>%
  group_by(node) %>%
  summarise(red = median(red, na.rm = T))
################################################################################


################################################################################
#Estimate clade counts from full tree (including unbinned contigs)

#Convert edge list to extract parent node for each internal node
curr.ancestors <- tibble(
  node = c(
    rootedByMidPointLabeled$tip.label,
    rootedByMidPointLabeled$node.label
  )[rootedByMidPointLabeled$edge[,2]],
  ancestor = c(
    rootedByMidPointLabeled$tip.label,
    rootedByMidPointLabeled$node.label
  )[rootedByMidPointLabeled$edge[,1]]
) %>%
  filter(grepl("Node", node)) %>%
  arrange(ancestor) %>%
  group_by(ancestor) %>%
  mutate(
    sibling = rev(node)
  ) %>%
  ungroup

#Convert edge list into descendants list
writeLines(paste(date(), marker_id, "=> get pre-populated node data (THIS CAN TAKE TIME!!!)"))
curr.nodes.full <- tibble(
  node = c(
    rootedByMidPointLabeled$tip.label,
    rootedByMidPointLabeled$node.label
  )[rootedByMidPointLabeled$edge[,1]],
  descendant = c(
    rootedByMidPointLabeled$tip.label,
    rootedByMidPointLabeled$node.label
  )[rootedByMidPointLabeled$edge[,2]]
) %>%
  #Add ancestor info
  left_join(curr.ancestors, by = "node") %>%
  #Add RED values
  left_join(avg.red, by = "node") %>%
  #Add RED values of descendants
  left_join(avg.red, by = c("descendant" = "node"), suffix = c("", ".desc")) %>%
  #Add RED values of ancestor
  left_join(avg.red, by = c("ancestor" = "node"), suffix = c("", ".anc")) %>%
  #Add REF values of sibling
  left_join(avg.red, by = c("sibling" = "node"), suffix = c("", ".sib")) %>%
  arrange(red) %>%
  group_by(node) %>%
  mutate(
    red.desc.mean = mean(red.desc),
    red.sib.mean = mean(c(red, red.sib))
  )
writeLines(paste(date(), marker_id, "=> done with pre-populated node data."))

#Iterate through tax levels and process full tree at specified RED values
#=> these values are taken from Parks et al 2018 
red.cutoffs <- tibble(
  level = c("phylum", "class", "order", "family", "genus"),
  red = c(0.32, 0.46, 0.62, 0.77, 0.93)
)

collect.clade_sizes <- tibble()
collect.h5_counts.bac <- collect.source_counts.bac <- tibble()
for (lev in c("phylum", "class", "order", "family", "genus")) {
  writeLines(paste(date(), marker_id, "=> full tree,", lev))
  
  curr.cutoff <- red.cutoffs %>%
    filter(level == lev) %>%
    pull(red)
  
  #Cut tree
  writeLines(paste(date(), marker_id, "=>", lev, "cut tree (THIS WILL TAKE TIME!!!)"))
  cut.nodes <- curr.nodes.full %>%
    #Filter by RED value
    #=> exclude edge cases (RED == 0 for root and RED == 1 for tips)
    filter(red > 0 & red < 1) %>%
    #Get absolute RED difference from cutoff for preferential choice of nodes
    mutate(
      diff_from_cut = abs(red - curr.cutoff),
      diff_from_cut.desc = abs(red.desc.mean - curr.cutoff),
      diff_from_cut.anc = abs(red.anc - curr.cutoff),
      diff_from_cut.sib = abs(red.sib.mean - curr.cutoff)
    ) %>%
    #Filter node set
    #=> keep nodes that are closer to RED cutoff than their ancestor AND the mean of their descendants
    mutate(
      #Is this node closer to the cutoff than (the mean of) its descendants?
      closer_than_desc = diff_from_cut < diff_from_cut.desc | is.na(diff_from_cut.desc),
      #Is this node at least as close to the cutoff as its ancestor?
      closer_than_anc = diff_from_cut <= diff_from_cut.anc | is.na(diff_from_cut.anc),
      #Is the mean distance of this node and its sibling at least as close to the cutoff as their common ancestor?
      closer_than_anc_with_sib = diff_from_cut.sib <= diff_from_cut.anc & !is.na(diff_from_cut.sib),
      #Is this node's sibling a tip?
      sibling_is_tip = !grepl("Node", sibling),
      #Is this nodes's sibling within the current RED cutoff range?
      sibling_within_range = abs(red.sib - curr.cutoff) <= 0.1
    ) %>%
    filter(closer_than_desc & closer_than_anc_with_sib)
  writeLines(paste(date(), marker_id, "=> done cutting tree for", lev))
  
  if (nrow(cut.nodes) == 0) {next()}
  
  #Filter node set to remove ancestral nodes that could be (accidentally) left over
  use.nodes.clean <- cut.nodes %>%
    select(node) %>%
    distinct() %>%
    #Filter node set to remove ancestral nodes (to avoid checking multiple times)
    rowwise() %>%
    mutate(
      ancestors = list(rootedByMidPointLabeled$node.label[Ancestors(rootedByMidPointLabeled, node, type = "all") - Ntip(rootedByMidPointLabeled)])
    ) %>%
    unnest(cols = c(ancestors)) %>%
    ungroup() %>%
    filter(!node %in% ancestors) %>%
    pull(node) %>%
    unique()
  
  cut.nodes.clean <- cut.nodes %>%
    filter(node %in% use.nodes.clean)
  
  if (nrow(cut.nodes.clean) == 0) {next()}
  
  #Get descendant tips for each selected node
  writeLines(paste(date(), marker_id, "=>", lev, "get descendant tips for selected nodes (THIS WILL TAKE TIME!!!)"))
  map.tips <- cut.nodes.clean %>%
    select(node) %>%
    distinct() %>%
    #Get descendants (tips only) for each retained node
    #=> this "Descendants" function seems to be slow-ish and acting up. Maybe there are alternatives?
    rowwise() %>%
    mutate(
      tmp.tips = Descendants(rootedByMidPointLabeled, node, "tips")
    ) %>%
    unnest(cols = c(tmp.tips))
  writeLines(paste(date(), marker_id, "=> done getting descendant tips for", lev))
  
  #Find how many tips are not covered by current node set and
  #treat "dropped" tips as singleton "nodes" => they form their own group
  dropped.tips <- tibble(
    node = rootedByMidPointLabeled$tip.label[!1:Ntip(rootedByMidPointLabeled) %in% map.tips$tmp.tips],
    tmp.tips = (1:Ntip(rootedByMidPointLabeled))[!1:Ntip(rootedByMidPointLabeled) %in% map.tips$tmp.tips]
  )
  
  #Check if there are residual duplicate tips
  duplicate.tips <- map.tips %>% group_by(tmp.tips) %>% filter(n() > 1) %>% ungroup() %>% select(tmp.tips) %>% distinct()
  
  #Merge structures, add tip source
  map.tips.clean <- bind_rows(
    map.tips,
    dropped.tips
  ) %>%
    rowwise() %>%
    mutate(raw.tip.label = rootedByMidPointLabeled$tip.label[tmp.tips]) %>%
    left_join(dat.tree.tips %>% select(raw.tip.label, source), by = "raw.tip.label")
  
  #Compute clade sizes
  curr.clade_sizes <- map.tips.clean %>%
    #Get clade sizes
    group_by(node, source) %>%
    tally() %>%
    mutate(level = lev) %>%
    select(level, node, source, n)
  
  #Store clade sizes
  collect.clade_sizes <- bind_rows(
    collect.clade_sizes,
    curr.clade_sizes
  )
  
  #Add tip sample name, then pull in sample data
  tmp.h5_counts <- map.tips.clean %>%
    rowwise() %>%
    mutate(raw.tip.label = rootedByMidPointLabeled$tip.label[tmp.tips]) %>%
    left_join(dat.tree.tips %>% select(raw.tip.label, sample_name, cluster), by = "raw.tip.label") %>%
    left_join(dat.clustering.by_gene %>% select(cluster, sample_name.full = sample_name, source.full = source), by = "cluster") %>%
    left_join(dat.sample.rarefy %>% select(sample_id, metalog_sample_id, ena_ers_sample_id, study_id, h5, h5.cat), by = c("sample_name.full" = "sample_id"))
  
  #Tabulate habitat counts
  curr.h5_counts <- tmp.h5_counts %>%
    group_by(node, source, h5) %>%
    tally %>%
    # group_by(source, h5) %>%
    # tally %>%
    mutate(
      marker_id = marker_id,
      level = lev
    ) %>%
    relocate(marker_id, level, source)
  
  #Tabulate gene source counts
  curr.source_counts <- tmp.h5_counts %>%
    group_by(node, source.full) %>%
    tally %>%
    mutate(
      marker_id = marker_id,
      level = lev
    ) %>%
    rename(source = source.full) %>%
    relocate(marker_id, level, source)
  
  #Store habitat counts
  collect.h5_counts.bac <- bind_rows(
    collect.h5_counts.bac,
    curr.h5_counts
  )
  
  #Store source counts
  collect.source_counts.bac <- bind_rows(
    collect.source_counts.bac,
    curr.source_counts
  )
  
  #Export
  write.table(curr.clade_sizes, file=paste0(folder.store, marker_id, "_", lev, "_curr.clade_sizes_bac.tab"), sep="\t", row.names=F, col.names=T, quote=F)
  save(collect.h5_counts.bac, collect.source_counts.bac, file = paste0(folder.store, marker_id, ".h5_counts.Rdata"))
}
################################################################################


################################################################################
#Store data

write.table(collect.clade_sizes, file=paste0(folder.store, marker_id,"_collect.clade_sizes_bac.tab"), sep="\t", row.names=F, col.names=T, quote=F)

write.table(convert.clade_counts, file=paste0(folder.store, marker_id,"_convert.clade_counts_bac.tab"), sep="\t", row.names=F, col.names=T, quote=F)

# save the workspace 
save.image(file=paste0(folder.store, marker_id,'_workspaceTaxonomyPredictions_bac.RData'))
################################################################################

####
