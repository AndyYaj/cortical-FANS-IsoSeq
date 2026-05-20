# IsoSeq-FANS-Human-Cortex

Authors: Andy Yang, Miguel Rodriguez de los Santos, Alexy Kozenkov, Stella Dracheva, Jack Humphrey, Michael S. Breen

Affiliation: Icahn School of Medicine at Mount Sinai, New York, NY 

Correspondance: michael.breen@mssm.edu

## Cell type-resolved long-read transcriptomics of the human cortex reveals pervasive alternative splicing and disease-relevant isoforms

Alternative splicing generates extensive transcriptomic diversity in the human brain, but the full cell type-resolved landscape of isoform variation remains unresolved due to the constraints of short-read sequencing. Here, we integrated fluorescence-activated nuclei sorting with long-read (PacBio Iso-Seq) and short-read (Illumina) RNA sequencing to generate isoform-resolved transcriptomes of five major cortical cell types, including MGE-derived GABAergic neurons, glutamatergic neurons, oligodendrocytes, astrocytes, and microglia from adult dorsolateral and orbitofrontal cortex. 

<img width="9921" height="8315" alt="Figure_01_revised" src="https://github.com/user-attachments/assets/d4c1c751-2651-4803-930f-0e88013c0281" />

We identified more than 220,000 unique full-length isoforms, 35–56% of which were previously unannotated depending on cell type and region. These novel isoforms were longer, contained more exons, and displayed greater coding potential than annotated transcripts. This resource provides a comprehensive, cell type–resolved atlas of isoform diversity in the human cortex and establishes a foundation for mechanistic studies of RNA regulation and disease vulnerability in the brain.

This repo contains all custom code and resources used to process, analyze and visualize data for this study. The raw PacBio IsoSeq and Illumina Short-Read RNA-seq data will be made available shortly on X. Additionally, processed and summarized data can be found on Zonedo and on the supplemental section of our pre-print. Lastly, cell type-resolved isoforms can be visualized on the [UCSC genome browser](https://genome.ucsc.edu/) using our [BED files](cell_type_resolved_isoforms.zip).

## Pipeline Overview

Data analysis largely consists of three main phases: 
1. [Data Processing](#1-data-processing)
2. [Constructing and Annotating Refined Isoforms](#2-isoform-annotations)
3. [Running Downstream Analyses](#3-running-downstream-analyses)
  
<img width="9921" height="7323" alt="FIGURE_02_revised" src="https://github.com/user-attachments/assets/3a427732-cb9d-4db5-bda4-749497cf320c" />

## 1. Code to produce the manuscript figures

| Figure | Description | Path |
|--------|-------------|------|
| 01 | Description of figure 1 | `path/to/figure1/` |
| 02 | Description of figure 2 | `path/to/figure2/` |
| 03 | Description of figure 3 | `path/to/figure3/` |
| 04 | Description of figure 4 | `path/to/figure4/` |
| 05 | Description of figure 5 | `path/to/figure5/` |
| 06 | Description of figure 6 | `path/to/figure6/` |
| 07 | Description of figure 7 | `path/to/figure7/` |
| 08 | Description of figure 8 | `path/to/figure8/` |

## 2. Data accessibility

| Data | Description | Reference |
|--------|-------------|------|
| IsoSeq Long Reads | Raw IsoSeq long reads | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
| DLPFC Short Reads | Complementay raw Illumina short read | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
| OFC Short Reads | Complementay raw Illumina short read | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
| FANS NueN+/- Short Reads | Splice junction support | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
