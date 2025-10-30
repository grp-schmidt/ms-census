# Analysis Code: a census of hidden and discoverable microbial diversity beyond genome-centric approaches

This repository contains analysis code for the manuscript on "A census of hidden and discoverable microbial diversity beyond genome-centric approaches" by [Prasoodanan, Maistrenko, et al](https://www.biorxiv.org/content/10.1101/2025.06.26.661807v1).

## Abstract

Cataloguing Earth’s biodiversity remains one of biology’s most formidable challenges, and the greatest diversity is expected to reside among the smallest organisms: microbes. Yet the ongoing census of microbial life is hampered by disparate sampling of Earth’s habitats, challenges in isolating uncultivated organisms, limited resolution in taxonomic marker gene amplicons, and incomplete recovery of metagenome-assembled genomes (MAGs). Here, we quantified discoverable bacterial and archaeal diversity in a comprehensive, curated cross-habitat dataset of 92,187 metagenomes. Clustering 502M sequences of 130 marker genes, we detected 705k bacterial and 27k archaeal species-level clades, the vast majority of which was hidden among ‘unbinned’ contigs, beyond current genome-centric approaches. At deeper taxonomic levels, we estimate that 10 archaeal and 145 bacterial novel phyla and around 80k novel genera are discoverable in current data. We identified soils and aquatic environments as novel lineage recovery hotspots, yet predict that discovery will remain in full swing across habitats as more data accrues.

## Description

### Data sources

The analyses are based on assembled contigs and genomes (based on cultivated isolates or as Metagenome-Assembled Genomes) from three sources:

proGenomes v3 (**PG3**), Fullam et al (Nucl Acid Res, 2023). [paper](https://academic.oup.com/nar/article/51/D1/D760/6835361). [website](https://progenomes.embl.de).

Genome Taxonomy Database (**GTDB**) r220, Parks et al (Nucl Acid Res, 2022). [paper](https://academic.oup.com/nar/article/50/D1/D785/6370255). [website](https://gtdb.ecogenomic.org/stats/r220).

Searchaeble Planetary-scale mIcrobiome REsource (**SPIRE**) v1.1, Schmidt et al (2024). [paper](https://academic.oup.com/nar/article/52/D1/D777/7332059). [website](spire.embl.de).

From the preprint:

> Genomes from all three datasets were downloaded and taxonomically (re-)classified against the GTDB r220 using GTDB-tk v2.4.0, and consensus taxonomies for species-level clusters were inferred based on adjusted majority votes, as described previously. Compared to SPIRE v1.0, 6,959 metagenomic samples were excluded for the present analysis based on data type and provenance (e.g., excluding samples explicitly enriched for viruses) or insufficient assembly size. For the remaining 92,187 samples, habitat annotations and contextual data were updated by manually curating against an extended microntology v0.3.0 encompassing 103 terms and categories [...] and further organised into 32 higher-level categories.

---

### Extraction of marker genes.

> The majority of analyses presented in the main text are based on 168 near-universal taxonomic marker genes as established by the GTDB: 120 bacterial (‘bac120’, [Parks et al, 2017](https://www.nature.com/articles/s41564-017-0012-7)) and 53 archaeal (‘arc53’, [Dombrowski et al, 2020](https://pubmed.ncbi.nlm.nih.gov/32770105/)) markers, with an overlap of five genes. We downloaded profile Hidden Markov Models (HMMs) for these marker sets as part of the GTDB-tk v2.4.0 database and used the HMMER v3.4 hmmsearch routine to identify and extract marker gene sequences among predicted Open Reading Frames (ORFs).

The corresponding HMM models were downloaded via the GTDB-tk release v2. They remain, to our knowledge, unchanged in the GTDB since at least release GTDB r207 (current at time of writing is r226). Metadata on both sets is available online:

`arc53`: [ar53_msa_marker_info_r207.tsv](https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/auxillary_files/ar53_msa_marker_info_r207.tsv)

`bac120`: [bac120_msa_marker_info_r207.tsv](https://data.ace.uq.edu.au/public/gtdb/data/releases/release207/207.0/auxillary_files/bac120_msa_marker_info_r207.tsv)

All extracted marker gene sequences are available for download [here](spire.embl.de/downloads)

---

### Clustering of marker genes.

Extracted marker gene sequences were clustered at various similarity thresholds ranging from 95% to 99% using the [MMseqs2]([https://www.nature.com/articles/s41467-018-04964-5](https://mmseqs.com/latest/userguide.pdf)) `cascading clustering` routine. For further analyses, the 96.5% threshold was selected (see preprint).

Commands used:

```
mmseqs creatdb [marker_id].faa [marker_id].db
mmseqs cluster [marker_id].db [marker_id].clustered tmp.clustering --min-seq-id [cutoff] --threads 24 --split-memory-limit 200G
mmseqs createtsv [marker_id].db [marker_id].db [marker_id].clustered [marker_id].clustered.tsv
```

_Note_: `[marker_id].faa` contain the extracted sequences for *one* marker gene at a time, from all three data sources (PG3, GTDB r220 and SPIRE).

---

### Parse clustering data and infer cluster_count-to-species_count conversion factors

Rationale (from the preprint):

> To convert the number of marker gene sequence clusters into a corresponding number of species-level clusters, we estimated conversion factors as follows. First, we generated marker gene cluster discovery curves via iterative logarithmic rarefaction, i.e. we downsampled the number of considered gene sequences along a logarithmic scale (10, 20,…, 100 sequences; 200, 300,…, 1000 sequences; etc) with 10 iterations at each step. At each rarefaction point, we recorded the number of ‘discovered’ marker gene sequence clusters and the number of represented species or species-level genome clusters in proGenomes3, GTDB r220 and among SPIRE MAGs, considering each data source individually. We then used linear regression models of the type number_of_species ∼ number_of_gene_clusters with a forced intercept at 0 along these rarefactions to estimate gene cluster to species conversion factors (i.e., the number of newly discovered species per newly discovered marker gene cluster). Based on benchmarks of marker gene sequence similarity cutoffs ranging from 95% to 99.5%, we based further analyses on 96.5% clusters as these showed very robust linear fits (with standard errors in the range of 10^-3) with conversion factors closest to identity (i.e., roughly one species discovered per marker gene cluster discovered), and high consistency across the different species-level reference clusterings in the underlying datasets (based on 40 specI marker genes in proGenomes3, 95% whole-genome ANI in GTDB and a combination of both approaches among SPIRE MAGs).

The corresponding code can be found in [01.parse_clustering.Rmd](https://github.com/grp-schmidt/ms-census/blob/main/R/01.parse_clustering.Rmd).

The resulting estimated conversion factors approximate real species counts well:

![image](https://github.com/user-attachments/assets/e8c88322-f4b2-4942-b647-dc7583d9dcf8)

---

### Perform habitat-resolved and habitat-incremental rarefaction analyses.

Rationale (from the preprint):

> Inference of global and habitat-resolved species discovery curves
> 
> For each of 32 broadly defined habitat categories (ranging from 267 to 19,659 samples per group [...]), as well as for the combined set of all 92,187 metagenomes in SPIRE v1.1, we generated bacterial and archaeal species discovery curves. We iteratively subsampled the number of considered metagenomes along a logarithmic rarefaction scale, with 5 random permutations per step. At each rarefaction point and for each considered marker gene, we inferred the number of discovered species based on the number of discovered marker gene clusters, using gene-specific conversion factors as described above. Rarefaction permutations were averaged per marker gene and then summarised per marker gene set (bac120 for Bacteria and arc53 for Archaea) as median predicted species counts at each rarefaction step across marker genes.
> 
> In an independent approach, we performed an ‘incremental’ rarefaction, stratified by habitats: we sequentially added samples from each of the 32 broadly defined habitat categories, starting with ‘adult, healthy human gut’ samples, randomized discovery order within each habitat block and then tracked the globally discovered species (see Figure 2). In other words, we tracked the contribution of each habitat to overall discovered species diversity, beyond the diversity discovered in previously considered habitats.

The corresponding code can be found in [02.job.perform_rarefactions.R](https://github.com/grp-schmidt/ms-census/blob/main/R/02.job.perform_rarefactions.R)

---

### Analyse species discovery curves (including species discovery coefficients

Rationale (from the preprint):

> Counts of discovered marker gene clusters (and, by proxy, inferred discovered species) were hierarchically stratified by data source to account for overlap between datasets as follows: all clusters containing at least one sequence originating from a genome in proGenomes3 were labelled as ‘proGenomes3’, irrespective of the origin of other sequences within the same cluster (data series marked in dark blue in main figures 1 & 2); clusters containing sequences from the GTDB r220 (and any other data source except proGenomes3) were labelled as ‘GTDB’ (light blue in figures 1 & 2); clusters containing sequences from SPIRE MAGs (but neither from proGenomes3 nor GTDB) were labelled as ‘SPIRE MAG’ (orange data series); clusters containing only sequences from unbinned SPIRE contigs were labelled as ‘unbinned, non-singleton’ if they contained at least two sequences from different contigs (dark yellow data series), and ‘unbinned, all’ otherwise (light yellow data series). Thus, clusters labelled as ‘GTDB’ represent sequence diversity contained in GTDB r220 beyond that already represented in proGenomes3, while a significant subset of ‘proGenomes3’-labelled clusters also contained sequences originating from GTDB, SPIRE MAGs or unbinned contigs. Similarly, ‘SPIRE MAG’ clusters corresponded to diversity not already covered by proGenomes3 and GTDB r220; and ‘unbinned’ clusters encompassed diversity not represented in any genome or MAG in the dataset. Moreover, by design we only tracked marker gene clusters in metagenomic samples, meaning that a non-negligible subset of genome-based clusters from proGenomes3 and GTDB r220 were not represented in rarefaction runs and calculations, as corresponding sequences were not detected in any of the 92,187 considered metagenomic assemblies.

Further, on `species discovery coefficients`:

> Calculation of species discovery coefficients
> 
> We quantified the rate at which novel species-level clusters are discovered in each habitat (and globally, across all habitats) using the equation
>
> S = k * N ^ -γ
> 
>where S is the number of newly ‘discovered’ species per N samples added to the survey; k is a proportionality constant; and γ is a saturation coefficient. An analogous formula is commonly used to describe a (microbial) pangenome’s openness 68,69, i.e. the degree to which novel genes are discovered as more genomes from the same species are considered. Moreover, the approach is conceptually and mathematically related to the Heridan-Heaps law in linguistics which describes the number of distinct words in a document as a function of that document’s length 70.
> 
> We fit equation (1) to each habitat-specific and global rarefaction curve described above, stratified by data source, resulting in estimates for γ for each considered marker gene. For more intuitive interpretability, we calculated a species discovery coefficient α as
> 
> α = 1 - γ
> 
> and summarised values across marker genes within each domain-specific set as median values. Thus defined, α scales on [-∞, 1], although only mildly negative values are expected to be observed in practice. If α ≤ 0, species discovery in a habitat is saturated, meaning that adding more samples of the same type is not expected to add novel species to the survey – analogous to a ‘closed’ microbial pangenome where additional strains do not add novel genes. Values for α in [0, 1] correspond to unsaturated species discovery curves where additional samples continue to add novel species to the survey (analogous to ‘open’ pangenomes); lower α values indicate a more pronounced ‘flattening off’ in the rarefaction curve, indicating a more pronounced slowdown in novel species discovery; higher α values indicate a less pronounced decrease in the rate of species discovery. For α -> 1, species discovery is fully unsaturated, meaning that each newly added sample adds novel species to the survey, with no discernible ‘flattening off’ in the species discovery curve.

The corresponding code can be found in [03.analyze_rarefactions.Rmd](https://github.com/grp-schmidt/ms-census/blob/main/R/03.analyze_rarefactions.Rmd). The code was used to generate the plots underlying figures 1, 2 & S1-S3 in the manuscript.

---

### Inference of marker gene phylogenies.

Rationale (from the preprint):

> We first re-aligned cluster representative amino acid sequences against the respective marker gene HMMs using HMMER v3.4’s `hmmalign` routine, trimmed the resulting pseudo-multiple sequence alignments to informative columns using `clipkit`  in `kpi-gappy` mode and removed sequences with >70% gaps among the remaining columns. Based on the resulting alignments, we inferred phylogenetic trees using `FastTree2` under the Whelan-Goldman model.

Example code (run for each marker gene). Alignment, pruning & filtering:

```
hmmalign --amino --trim --outformat pfam [marker_id].HMM [marker_id].rep.faa |
  grep -v "^#=G" |
  perl -ane 'if (/^#/) {print} else {($acc, $seq) = split(/\\s+/); $seq =~ s/[a-z]/./g; print \"$acc $seq\\n\"}' > [marker_id].rep.faa.stk

clipkit -m kpi-gappy [marker_id].rep.faa.stk
cat [marker_id].rep.faa.stk.clipkit | grep -v "^#=G" |
  perl -ane 'if (/#/) {print} else {($acc, $seq) = split /\\s+/; $count = $seq =~ tr/[A-Z,-]//; $count_char = $seq =~ tr/[A-Z]//; if ($count_char / $count > 0.7) {print \"$acc $seq\\n\"}}' > [marker_id].rep.faa.stk.clipkit.filt
```

Phylogeny inference:

```
cat [marker_id].rep.faa.stk.clipkit.filt |
  grep -v "^#" |
  grep -v "#/" | 
  perl -ane 'print \">\"; s/ /\\n/; print' > [marker_id].rep.faa.stk.clipkit.filt.fa
fasttree -wag -nopr [marker_id].rep.faa.stk.clipkit.filt.fa > [marker_id].nwk
```

The resulting phylogenies (for a filtered subset of marker genes) and corresponding tip-level metadata were uploaded to the EBI BioStudies repository (see data availability statement).

---

### Relative Evolutionary Divergence (RED)-based inference and analysis of deeper clades

Rationale (from the preprint):

> We calculated relative evolutionary divergence (RED) values (Parks et al, 2018) for each node in the resulting trees using the castor package v1.8.3 in R. To infer marker-gene specific relative evolutionary divergence (RED) cutoffs at phylum, class, order, family, and genus levels, we subset phylogenies to sequences originating from fully taxonomically classified genomes in proGenomes3, the GTDB r220 or among SPIRE MAGs, cut the trees at incremental RED values (with a RED tolerance of ±0.1) and compared the resulting partitions (i.e., tip sets) versus assigned taxonomic labels, quantifying partition consistency as adjusted mutual information (AMI). These steps were performed iteratively on trees re-rooted using each recognized phylum in turn as an outgroup and then summarised across alternative root positions (ignoring the current outgroup) to obtain average RED values per internal node, in a workflow adapted from (Parks et al, 2018).
> 
> After confirming that the resulting RED cutoffs at peak AMI (i.e., best capturing the distribution of established taxonomic labels) well approximated those reported in (Parks et al, 2018) across individual archaeal marker genes, we used GTDB’s reference cutoffs (with a ±0.1 tolerance interval) for bacterial marker gene trees and cut the full trees at the respective levels. The number of resulting clusters (i.e., nodes where the tree was cut) and their composition (i.e., each node’s descendant tips) were then used as estimates of clade-level groups at each taxonomic level, for each marker gene. We manually selected genes that showed the most consistent and robust profiles across iterations (and rootings) and whose phylogenies were not ostensibly impacted by paralogs, bringing the final sets of considered genes to 29 archaeal and 53 bacterial markers.

The corresponding code can be found in

* [04.job.infer_red_trees.R](https://github.com/grp-schmidt/ms-census/blob/main/R/04.job.infer_red_trees.R) for the inference of deeper clades based on RED cutoff, per marker gene (script works specifically for Bacteria, hard-coded defaults need to be adapted for Archaea)
*  [05.analyze_phylogenies.Rmd](https://github.com/grp-schmidt/ms-census/blob/main/R/05.analyze_phylogenies.Rmd) to generate analyses and plots underlying Figure 4 in the manuscript.

---

## Availability of underlying data

Metagenomic assemblies, Metagenome-Assembled Genomes (MAGs), gene calls and corresponding annotations are available via the SPIRE db [downloads page](spire.embl.de/downloads).

Inferred marker gene phylogenies with annotations, as well as pre-generated tree visualizations for archaeal markers are available via the EBI BioStudies repository the following accessions:

| Accession | Description |
| ---- | ---- |
| [S-BSST2111](https://www.ebi.ac.uk/biostudies/studies/S-BSST2111) | 28 pre-generated plots of archaeal marker gene phylogenies, as PDF |
| [S-BSST2112](https://www.ebi.ac.uk/biostudies/studies/S-BSST2112) | 53 bacterial marker gene phylogenies in `nwk` format |
| [S-BSST2113](https://www.ebi.ac.uk/biostudies/studies/S-BSST2113) | Tip-level annotations and metadata for 53 bacterial marker gene phylogenies |
| [S-BSST2116](https://www.ebi.ac.uk/biostudies/studies/S-BSST2116) | Tip-level annotations and metadata for 28 bacterial marker gene phylogenies |
| [S-BSST2117](https://www.ebi.ac.uk/biostudies/studies/S-BSST2117) | 28 archaeal marker gene phylogenies in `nwk` format |

Additional data, e.g. sample-level annotations, are available as supplementary material to the [preprint](https://www.biorxiv.org/content/10.1101/2025.06.26.661807v1.supplementary-material) on bioRxiv.

## Contributors

Lead authors:

Vishnu Prasoodanan P K, APC Microbiome & School of Medicine, University College Cork, Ireland. [ORCID: 0009-0002-1872-5790](https://orcid.org/0009-0002-1872-5790)

Oleksandr M Maistrenko, Department of Marine Microbiology & Biogeochemistry, Royal Netherlands Institute for Sea Research(NIOZ), The Netherlands. [ORCID: 0000-0003-1961-7548](https://orcid.org/0000-0003-1961-7548).


Correspondence:

Thomas S B Schmidt, APC Microbiome & School of Medicine, University College Cork, Ireland. [ORCID: 0000-0001-8587-4177](https://orcid.org/0000-0001-8587-4177).

sebastian _dot_ schmidt _at_ ucc _dot_ ie

