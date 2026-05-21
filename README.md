# IsoSeq-FANS-Human-Cortex

Authors: Andy Yang, Miguel Rodriguez de los Santos, Alexy Kozenkov, Ramu Vadukapuram, Yasmin Hurd, Stella Dracheva, Jack Humphrey, Michael S. Breen

Affiliation: Icahn School of Medicine at Mount Sinai, New York, NY 

Correspondance: michael.breen@mssm.edu

## Cell type-resolved long-read transcriptomics of the human cortex reveals pervasive alternative splicing and disease-relevant isoforms

Alternative splicing generates extensive transcriptomic diversity in the human brain, but the full cell type-resolved landscape of isoform variation remains unresolved due to the constraints of short-read sequencing. Here, we integrated fluorescence-activated nuclei sorting with long-read (PacBio Iso-Seq) and short-read (Illumina) RNA sequencing to generate isoform-resolved transcriptomes of five major cortical cell types, including MGE-derived GABAergic neurons, glutamatergic neurons, oligodendrocytes, astrocytes, and microglia from adult dorsolateral and orbitofrontal cortex. 

<img width="2340" height="2340" alt="Graphical_Abstract" src="https://github.com/user-attachments/assets/50a82d5d-646c-4900-a063-01445ca863fb" />

This repo contains all of the code and used to process, analyze and visualize data for this study. Additionally, processed and summarized data can be found on [Zonedo](https://zenodo.org/records/20301473) and on the supplemental section of our pre-print. Lastly, cell type-resolved isoforms can be visualized on the [UCSC genome browser](https://genome.ucsc.edu/) using our [BED files](cell_type_resolved_isoforms.zip) and through our [Rshiny application](https://andyyang-isoseq.share.connect.posit.cloud/).

## Pipeline Overview

Briefly, raw long reads generated on the PacBio Sequel II platform were processsed with IsoSeq3 to generate high-quality full-length transcripts, which were aligned to the GRCh38 genome (GENCODE v38) using minimap2. Transcript models were collapsed and merged across samples within each cell type in the DLPFC and OFC using TAMA to generate nonredundant transcript annotations, followed by classification and QC filtering with SQANTI3. Transcript abundance was quantified using Salmon.

<img width="9921" height="7465" alt="PIPELINE_OVERVIEW" src="https://github.com/user-attachments/assets/39b7096d-ae02-4dd3-90a2-4036cda9442a" />


## Code to produce the manuscript figures

All of the code used to produce the main figures in the manuscript can be found in [code](https://github.com/AndyYaj/cortical-FANS-IsoSeq/tree/main/code).  All of the data used in those analyses can be found on [Zonedo](https://zenodo.org/records/20301473).

## Data accessibility

| Data | Description |
|--------|-------------|
| [IsoSeq Long Reads](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753) | Raw DLPFC and OFC IsoSeq long read (BAM) files. |
| [DLPFC Short Reads](https://www.synapse.org/#!Synapse:syn12034263) | Raw short-read RNA-seq data (cohort 1) — Human DLPFC |
| OFC Short Reads | Raw short-read RNA-seq data (cohort 2) — Human OFC (GEO accession ID in progress) | 
| [FANS NueN+/- Short Reads](https://www.synapse.org/Synapse:syn25716684/wiki/610496) | Raw short-read RNA-seq data (cohort 3) — Human BM10, BM17, BM22, BM36, BM44. | 
