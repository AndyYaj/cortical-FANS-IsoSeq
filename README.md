# IsoSeq-FANS-Human-Cortex

Authors: Andy Yang, Miguel Rodriguez de los Santos, Alexy Kozenkov, Stella Dracheva, Jack Humphrey, Michael S. Breen

Affiliation: Icahn School of Medicine at Mount Sinai, New York, NY 

Correspondance: michael.breen@mssm.edu

## Cell type-resolved long-read transcriptomics of the human cortex reveals pervasive alternative splicing and disease-relevant isoforms

Alternative splicing generates extensive transcriptomic diversity in the human brain, but the full cell type-resolved landscape of isoform variation remains unresolved due to the constraints of short-read sequencing. Here, we integrated fluorescence-activated nuclei sorting with long-read (PacBio Iso-Seq) and short-read (Illumina) RNA sequencing to generate isoform-resolved transcriptomes of five major cortical cell types, including MGE-derived GABAergic neurons, glutamatergic neurons, oligodendrocytes, astrocytes, and microglia from adult dorsolateral and orbitofrontal cortex. 

<img width="9921" height="8315" alt="Figure_01_revised" src="https://github.com/user-attachments/assets/d4c1c751-2651-4803-930f-0e88013c0281" />

We identified more than 220,000 unique full-length isoforms, 35–56% of which were previously unannotated depending on cell type and region. These novel isoforms were longer, contained more exons, and displayed greater coding potential than annotated transcripts. Glial populations, particularly oligodendrocytes and microglia, exhibited the greatest isoform diversity, with more than half of all isoforms displaying strong cell type–specific expression. Differential transcript usage revealed pervasive cell-specific splicing, including in genes central to neuronal and glial function such as SLC5A6 and TWF1. Many newly discovered isoforms overlapped de novo variants in autism spectrum disorder–associated genes, including POGZ, FOXP1, and DYRK1A, suggesting that isoform diversification may contribute to neurodevelopmental risk. This resource provides a comprehensive, cell type–resolved atlas of isoform diversity in the human cortex and establishes a foundation for mechanistic studies of RNA regulation and disease vulnerability in the brain.

This repo contains all custom code and resources used to process, analyze and visualize data for this study. The raw PacBio IsoSeq and Illumina Short-Read RNA-seq data will be made available shortly on X. Additionally, processed and summarized data can be found on Zonedo and on the supplemental section of our pre-print. Lastly, cell type-resolved isoforms can be visualized on the [UCSC genome browser](https://genome.ucsc.edu/) using our [BED files](cell_type_resolved_isoforms.zip).

## Pipeline Overview

Data analysis largely consists of three mian phases: 
1. [Data Processing](#1-data-processing)
2. [Constructing and Annotating Refined Isoforms](#2-isoform-annotations)
3. [Running Downstream Analyses](#3-running-downstream-analyses)
   - Differential analyses and isoform switching
   - Predicted RNA secondary structures and Minimum Free Energy (MFE)
   - Weighted gene co-expression dynamics of isoforms
   - Neurological disease enrichment
  
<img width="9921" height="7323" alt="FIGURE_02_revised" src="https://github.com/user-attachments/assets/3a427732-cb9d-4db5-bda4-749497cf320c" />

## 1. Data Processing

- Illumina short-read RNA-sequencing data (RAPiD)
- PacBio long-read RNA-sequencing data (IsoSeq3, Minimap2, TAMA)
- Additional FANS-derived neuronal and non-neuronal nuclei data (RAPiD)

## 2. Isoform Annotations

- Characterizing HQ isoforms (SQANTI3)
- Comparing isoform annotations across cell types and brain regions (gffcompare)
- Quantifying IsoSeq data using short-read data (Salmon)
- Validating isoform full-length models (5' and 3' end support)
- Predicting isoform protein-coding potential
- Characterizing alternative splicing events
- Visualizing isoforms

## 3. Running Downstream Analyses

- Differential gene expression analysis
- Differential transcript/isoform expression analysis
- Differential transcript usage
- Isoform switching
- Weighted gene co-expression analysis
- Neurological disease enrichment analysis
