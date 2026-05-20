# IsoSeq-FANS-Human-Cortex

Authors: Andy Yang, Miguel Rodriguez de los Santos, Alexy Kozenkov, Stella Dracheva, Jack Humphrey, Michael S. Breen

Affiliation: Icahn School of Medicine at Mount Sinai, New York, NY 

Correspondance: michael.breen@mssm.edu

## Cell type-resolved long-read transcriptomics of the human cortex reveals pervasive alternative splicing and disease-relevant isoforms

Alternative splicing generates extensive transcriptomic diversity in the human brain, but the full cell type-resolved landscape of isoform variation remains unresolved due to the constraints of short-read sequencing. Here, we integrated fluorescence-activated nuclei sorting with long-read (PacBio Iso-Seq) and short-read (Illumina) RNA sequencing to generate isoform-resolved transcriptomes of five major cortical cell types, including MGE-derived GABAergic neurons, glutamatergic neurons, oligodendrocytes, astrocytes, and microglia from adult dorsolateral and orbitofrontal cortex. 

<img width="9921" height="8315" alt="Figure_01_revised" src="https://github.com/user-attachments/assets/d4c1c751-2651-4803-930f-0e88013c0281" />

This repo contains all of the code and used to process, analyze and visualize data for this study. Additionally, processed and summarized data can be found on [Zonedo](https://zenodo.org/records/20301473) and on the supplemental section of our pre-print. Lastly, cell type-resolved isoforms can be visualized on the [UCSC genome browser](https://genome.ucsc.edu/) using our [BED files](cell_type_resolved_isoforms.zip).

## Pipeline Overview

Briefly, raw long reads generated on the PacBio Sequel II platform were processsed with IsoSeq3 to generate high-quality full-length transcripts, which were aligned to the GRCh38 genome (GENCODE v38) using minimap2. Transcript models were collapsed and merged across samples within each cell type in the DLPFC and OFC using TAMA to generate nonredundant transcript annotations, followed by classification and QC filtering with SQANTI3. Transcript abundance was quantified using Salmon.

## 1. Code to produce the manuscript figure

All of the code used to produce the main figures in the manuscript can be found in code.  All of the data used in those analyses can be found on [Zonedo](https://zenodo.org/records/20301473).

## 2. Data accessibility

| Data | Description | Reference |
|--------|-------------|------|
| IsoSeq Long Reads | Raw IsoSeq long reads | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
| DLPFC Short Reads | Complementay raw Illumina short read | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
| OFC Short Reads | Complementay raw Illumina short read | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
| FANS NueN+/- Short Reads | Splice junction support | `[External Link]([https://example.com/figure1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753))` |
