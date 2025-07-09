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

The corresponding code can be found in []().

The resulting estimated conversion factors approximate real species counts well:

![image](https://github.com/user-attachments/assets/e8c88322-f4b2-4942-b647-dc7583d9dcf8)


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

