#!/usr/bin/Rscript
#
#R script to perform rarefactions
#
#Input:
#[marker_id] [cutoff]
#
#Output:
#rarefaction results for different types of rarefaction runs
#
#***
#
#**NOTE**: To run code below, downlaod data from Zenodo (DOI 10.5281/zenodo.17482698) and store locally in the according folder structure.
#
#***
################################################################################
################################################################################


################################################################################
#Get input command line arguments
args <- commandArgs(trailingOnly = T)

#Get marker gene ID
mid <- args[1]
cutoff <- as.character(args[2])
################################################################################


################################################################################
#Load libraries
library("tidyverse", quietly = T)
library("tidyr", quietly = T)
library("stringr", quietly = T)

#Make R behave marginally less moronic
options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
################################################################################


################################################################################
#Set folder and file names, load data
#folder.base <- [insert_base_folder_here]
folder.data <- paste0(folder.base, "data/")

#Sample data
file.sample_data <- paste0(folder.base, "metadata/data.sample.clean.Rdata")
load(file.sample_data)

#Global clustering and conversion coefficient data
file.clustering_summary <- paste0(folder.data, "/00_SUMMARY.spire_gtdb.clustered.rarefaction.Rdata")
load(file.clustering_summary)

#HMM clustering data
file.gene_clustering <- paste0(folder.data, "/clustering/", mid, ".spire_gtdb.clustered.", cutoff, ".Rdata")
if (!file.exists(file.gene_clustering)) {writeLines(paste(date(), "=> FILE NOT FOUND: ", file.gene_clustering)); q()}
load(file.gene_clustering)
################################################################################


################################################################################
#Pre-process sample data
order.by_h5 <- c(
  "human gut (adult, healthy)",
  "human gut (adult, diseased)",
  "human gut (non-industrialized)",
  "human gut (infant)",
  "human gut (child)",
  "human gut (elderly)",
  "pig gut",
  "bovine gut",
  "mammalian gut (other)",
  "rumen",
  "bird gut",
  "animal gut (other)",
  "human oral",
  "human skin",
  "human airways",
  "human urogenital",
  "human (other)",
  "built environment",
  "wastewater",
  "rhizosphere",
  "phyllosphere",
  "plant host (other)",
  "soil",
  "wetland",
  "marine",
  "fresh water",
  "groundwater",
  "hot spring",
  "hydrothermal vent",
  "aquatic (other)",
  "air",
  "other"
)

dat.sample.rarefy <- dat.sample.clean

#Precompute cumulative sums
dat.sample.rarefy.cumsum <- dat.sample.rarefy %>%
  #Account for incomplete mappings
  group_by(h5) %>%
  summarise(
    size = n()
  ) %>%
  mutate(
    cumsum = cumsum(size)
  ) %>%
  mutate(
    use.cumsum = lag(cumsum),
    use.cumsum = ifelse(is.na(use.cumsum), 0, use.cumsum)
  )

dat.sample.rarefy.incremental_steps <- dat.sample.rarefy %>%
  group_by(h5) %>%
  summarise(
    step = unique(c(0, round(seq(0, 0.1 * n(), by = max(10, n() / 100))),  round(seq(0.1 * n(), n(), by = min(500, n() / 20))), n()))
  ) %>%
  left_join(dat.sample.rarefy.cumsum %>% select(h5, use.cumsum)) %>%
  ungroup() %>%
  filter(step > 0) %>%
  mutate(
    n.samples = step + use.cumsum
  )
################################################################################


################################################################################
#Do the actual stuff
################################################################################
#Preallocate
collect.rarefaction.by_sample <- collect.rarefaction.by_habitat_increment <- collect.baseline <- tibble()

#Define helper function to compute effective number of clusters found
get.q1 <- function(x) {
  x <- x[x > 0]
  x.rel <- x / sum(x)
  exp(sum( - x.rel * log(x.rel)))
}

#Report
writeLines(paste(date(), "=>", mid, cutoff))

############################################################################
#Tally up clusters to hierarchically define baseline levels
collect.baseline <- dat.clustering.by_cluster %>%
  filter(source != "spire.unbinned") %>%
  pivot_wider(names_from = "source", values_from = "n", values_fill = 0) %>%
  summarise(
    marker_id = mid,
    n.pg3 = sum(pg3 > 0),
    n.gtdb_r220 = sum(gtdb_r220 > 0),
    n.gtdb_r220.unique = sum(gtdb_r220 > 0 & pg3 == 0),
    n.spire.mag = sum(spire.mag > 0),
    n.spire.mag.unique = sum(spire.mag > 0 & gtdb_r220 == 0 & pg3 == 0)
  ) %>%
  mutate(
    marker_id = mid,
    cutoff = cutoff
  ) %>%
  relocate(marker_id, cutoff)
############################################################################


############################################################################
#Rarefy the n of clusters discovered per sample added
writeLines(paste(date(), "=> performing incremental rarefactions"))
############################################################################
for (i in 1:5) {
  curr.rarefaction.incremental <- tibble()
  #Randomzie sample order, but only wihtin groups
  curr.samples.randomized <- dat.sample.rarefy %>%
    arrange(h5) %>%
    group_by(h5) %>%
    mutate(sample_id.shuffled = sample(sample_id)) %>%
    ungroup() %>%
    pull(sample_id.shuffled)
  
  #Pre-subset cluster data to (maximally) relevant subset
  rarefy.genes.incremental <- dat.clustering.by_gene %>%
    filter(sample_name %in% curr.samples.randomized) %>%
    mutate(source = ifelse(is.na(bin_id), "spire.unbinned", "spire.mag"))
  
  #Loop through pre-defined steps]
  for (curr.step in dat.sample.rarefy.incremental_steps$n.samples) {
    if (curr.step == 0) {next()}
    
    #Get info on current step
    curr.h5 <- dat.sample.rarefy.incremental_steps %>% filter(n.samples == curr.step) %>% slice_min(step) %>% pull(h5)
    curr.step.h5 <- dat.sample.rarefy.incremental_steps %>% filter(n.samples == curr.step) %>% slice_min(step) %>% pull(step)
    
    #Report
    writeLines(paste(date(), "=>", mid, cutoff, "incremental runs: iteration", i, curr.h5, curr.step.h5, "samples,", curr.step, "samples total."))
    
    #Get current samples
    curr.samples.i <- curr.samples.randomized[1:curr.step]
    
    #Pull current clusters
    curr.clusters <- rarefy.genes.incremental %>%
      filter(sample_name %in% curr.samples.i) %>%
      group_by(cluster, source) %>%
      tally
    
    #Get q1 diversity estimate on full clustering "discovered"
    #=> irrespective of whether genes were in MAGs or not
    curr.q1 <- curr.clusters %>%
      group_by(cluster) %>%
      summarise(n = sum(n)) %>%
      summarise(q1.total = get.q1(n)) %>%
      pull(q1.total)
    
    #Get more summarised data on clusters
    #=> n of "hierarchically" discovered clusters (unique in the unbinned relative to the MAG'ed fraction)
    curr.clusters.summarised <- curr.clusters %>%
      pivot_wider(names_from = "source", values_from = "n", values_fill = 0) %>%
      ungroup() %>%
      #add dummy column if spire.mag column doesn't exist
      mutate(
        spire.mag = if(! "spire.mag" %in% colnames(.)) 0 else spire.mag,
        spire.unbinned = if(! "spire.unbinned" %in% colnames(.)) 0 else spire.unbinned
      ) %>%
      summarise(
        n.total = n(),
        q1.total = curr.q1,
        n.spire.mag = sum(spire.mag > 0),
        q1.spire.mag = get.q1(spire.mag),
        n.spire.unbinned = sum(spire.unbinned > 0),
        q1.spire.unbinned = get.q1(spire.unbinned),
        n.spire.unbinned.unique = sum(spire.mag == 0 & spire.unbinned > 0),
        n.spire.unbinned.non_singleton = sum(spire.unbinned > 1),
        n.spire.unbinned.non_singleton.unique = sum(spire.mag == 0 & spire.unbinned > 1)
      )
    
    #Retrieve additional info on pulled clsuters
    curr.clusters.summarised.glo <- dat.clustering.by_cluster %>% 
      filter(cluster %in% curr.clusters$cluster) %>%
      pivot_wider(names_from = "source", values_from = "n", values_fill = 0) %>%
      #Add dummy columns if data happened to not exist
      mutate(
        pg3 = if(! "pg3" %in% colnames(.)) 0 else pg3,
        gtdb_r220 = if(! "gtdb_r220" %in% colnames(.)) 0 else gtdb_r220,
        spire.mag = if(! "spire.mag" %in% colnames(.)) 0 else spire.mag,
        spire.unbinned = if(! "spire.unbinned" %in% colnames(.)) 0 else spire.unbinned
      ) %>%
      summarise(
        glob.pg3 = sum(pg3 > 0),
        glob.gtdb_r220 = sum(gtdb_r220 > 0),
        glob.gtdb_r220.unique = sum(gtdb_r220 > 0 & pg3 == 0),
        glob.spire.mag = sum(spire.mag > 0),
        glob.spire.mag.unique = sum(spire.mag > 0 & gtdb_r220 == 0 & pg3 == 0),
        glob.spire.unbinned = sum(spire.unbinned > 0),
        glob.spire.unbinned.unique = sum(spire.unbinned > 0 & spire.mag == 0 & gtdb_r220 == 0 & pg3 == 0),
        glob.spire.unbinned.non_singleton = sum(spire.unbinned > 1),
        glob.spire.unbinned.non_singleton.unique = sum(spire.unbinned > 1 & spire.mag == 0 & gtdb_r220 == 0 & pg3 == 0)
      )
    
    tmp.pass <- bind_cols(
      curr.clusters.summarised,
      curr.clusters.summarised.glo
    ) %>%
      mutate(
        marker_id = mid,
        cutoff = cutoff,
        h5 = curr.h5,
        step = curr.step.h5,
        step.glob = curr.step,
        i = i
      ) %>%
      relocate(marker_id, cutoff, h5, step, step.glob, i)
    
    curr.rarefaction.incremental <- bind_rows(
      curr.rarefaction.incremental,
      tmp.pass
    )
  }
  
  ############################################################################
  #Pass results
  collect.rarefaction.by_habitat_increment <- bind_rows(
    collect.rarefaction.by_habitat_increment,
    curr.rarefaction.incremental
  )
  
  #Save
  #NOTE: rarefaction results for all marker genes (files stored in line below) are available via Zenodo as 'rarefaction.by_habitat.incremental.tar.gz'
  save(collect.rarefaction.by_habitat_increment, collect.baseline, file = paste0(folder.data, "/rarefaction.by_habitat.incremental.", mid, ".", cutoff, ".Rdata"))
  ############################################################################
}
############################################################################



############################################################################
#Rarefy the n of clusters discovered per sample added
writeLines(paste(date(), "=> performing per-habitat rarefactions"))
############################################################################
for (h5.rarefy in c(unique(as.character(dat.sample.rarefy$h5)), "all")) {
  #Subset sample data to relevant samples
  if (h5.rarefy == "all") {
    rarefy.samples <- dat.sample.rarefy %>%
      #Account for incomplete mappings
      #filter(sample_id %in% dat.clustering.by_gene$sample_name) %>%
      pull(sample_id)
  } else {
    rarefy.samples <- dat.sample.rarefy %>% filter(h5 == h5.rarefy) %>%
      #Account for incomplete mappings
      #filter(sample_id %in% dat.clustering.by_gene$sample_name) %>%
      pull(sample_id)
  }
  
  #Preallocate
  curr.rarefaction <- tibble()
  
  #Pre-subset cluster data to (maximally) relevant subset
  rarefy.genes <- dat.clustering.by_gene %>%
    filter(sample_name %in% rarefy.samples) %>%
    mutate(source = ifelse(is.na(bin_id), "spire.unbinned", "spire.mag"))
  
  #Predefine rarefaction steps
  curr.rarefaction.steps <- c(
    seq(10, 100, by = 10),
    seq(10^2, 10^3, by = 100),
    seq(10^3, 10^4, by = 1000),
    seq(10^4, 10^5, by = 10000),
    seq(round(0.1 * length(rarefy.samples)), length(rarefy.samples), by = round(0.1 * length(rarefy.samples))),
    length(rarefy.samples)
  )
  curr.rarefaction.steps <- unique(sort(curr.rarefaction.steps[curr.rarefaction.steps <= length(rarefy.samples)]))
  
  ############################################################################
  for (step in curr.rarefaction.steps) {
    writeLines(paste(date(), "=>", mid, cutoff, "per-habitat runs:", h5.rarefy, step))
    for (i in 1:5) {
      #Pull current clusters
      curr.clusters <- rarefy.genes %>%
        filter(sample_name %in% sample(rarefy.samples, step)) %>%
        group_by(cluster, source) %>%
        tally
      
      #Get q1 diversity estimate on full clustering "discovered"
      #=> irrespective of whether genes were in MAGs or not
      curr.q1 <- curr.clusters %>%
        group_by(cluster) %>%
        summarise(n = sum(n)) %>%
        summarise(q1.total = get.q1(n)) %>%
        pull(q1.total)
      
      #Get more summarised data on clusters
      #=> n of "hierarchically" discovered clusters (unique in the unbinned relative to the MAG'ed fraction)
      curr.clusters.summarised <- curr.clusters %>%
        pivot_wider(names_from = "source", values_from = "n", values_fill = 0) %>%
        ungroup() %>%
        #add dummy column if spire.mag column doesn't exist
        mutate( spire.mag = if(! "spire.mag" %in% colnames(.)) 0 else spire.mag) %>%
        mutate( spire.unbinned = if(! "spire.unbinned" %in% colnames(.)) 0 else spire.unbinned) %>%
        summarise(
          n.total = n(),
          q1.total = curr.q1,
          n.spire.mag = sum(spire.mag > 0),
          q1.spire.mag = get.q1(spire.mag),
          n.spire.unbinned = sum(spire.unbinned > 0),
          q1.spire.unbinned = get.q1(spire.unbinned),
          n.spire.unbinned.unique = sum(spire.mag == 0 & spire.unbinned > 0),
          n.spire.unbinned.non_singleton = sum(spire.unbinned > 1),
          n.spire.unbinned.non_singleton.unique = sum(spire.mag == 0 & spire.unbinned > 1)
        )
      
      #Retrieve additional info on pulled clsuters
      curr.clusters.summarised.glo <- dat.clustering.by_cluster %>% 
        filter(cluster %in% curr.clusters$cluster) %>%
        pivot_wider(names_from = "source", values_from = "n", values_fill = 0) %>%
        #Add dummy columns if data happened to not exist
        mutate(
          pg3 = if(! "pg3" %in% colnames(.)) 0 else pg3,
          gtdb_r220 = if(! "gtdb_r220" %in% colnames(.)) 0 else gtdb_r220,
          spire.mag = if(! "spire.mag" %in% colnames(.)) 0 else spire.mag,
          spire.unbinned = if(! "spire.unbinned" %in% colnames(.)) 0 else spire.unbinned
        ) %>%
        summarise(
          glob.pg3 = sum(pg3 > 0),
          glob.gtdb_r220 = sum(gtdb_r220 > 0),
          glob.gtdb_r220.unique = sum(gtdb_r220 > 0 & pg3 == 0),
          glob.spire.mag = sum(spire.mag > 0),
          glob.spire.mag.unique = sum(spire.mag > 0 & gtdb_r220 == 0 & pg3 == 0),
          glob.spire.unbinned = sum(spire.unbinned > 0),
          glob.spire.unbinned.unique = sum(spire.unbinned > 0 & spire.mag == 0 & gtdb_r220 == 0 & pg3 == 0),
          glob.spire.unbinned.non_singleton = sum(spire.unbinned > 1),
          glob.spire.unbinned.non_singleton.unique = sum(spire.unbinned > 1 & spire.mag == 0 & gtdb_r220 == 0 & pg3 == 0)
        )
      
      tmp.pass <- bind_cols(
        curr.clusters.summarised,
        curr.clusters.summarised.glo
      ) %>%
        mutate(
          marker_id = mid,
          cutoff = cutoff,
          h5 = h5.rarefy,
          step = step,
          i = i
        ) %>%
        relocate(marker_id, cutoff, h5, step, i)
      
      curr.rarefaction <- bind_rows(
        curr.rarefaction,
        tmp.pass
      )
    }
  }
  ############################################################################
  
  
  ############################################################################
  #Pass results
  collect.rarefaction.by_sample <- bind_rows(
    collect.rarefaction.by_sample,
    curr.rarefaction
  )
  
  #Save
  #NOTE: rarefaction results for all marker genes (files stored in line below) are available via Zenodo as 'rarefaction.by_habitat.tar.gz'
  save(collect.rarefaction.by_sample, collect.baseline, file = paste0(folder.data, "/rarefaction.by_habitat.", mid, ".", cutoff, ".Rdata"))
  ############################################################################
}
################################################################################


################################################################################

