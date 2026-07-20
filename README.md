# IsoSeq-FANS-Human-Cortex

Code for: **Long-read transcriptomics of purified human cortical cell types exposes glial isoform complexity and disease-relevant transcript architecture**

Authors: Yang A, Rodriguez de los Santos M, Kozlenkov A, Vadukapuram R, Hurd Y, Dracheva S, Humphrey J, Breen MS.

Affiliation: Icahn School of Medicine at Mount Sinai, New York, NY 

Preprint: https://www.biorxiv.org/content/10.1101/2025.11.25.690524v2

Interactive atlas: https://andyyang.shinyapps.io/Cortical_FANS_IsoSeq/

Archived release: https://doi.org/10.5281/zenodo.20301473

---

## Overview

Alternative splicing generates extensive transcriptomic diversity in the human brain, but the full cell type-resolved landscape of isoform variation remains unresolved due to the constraints of short-read sequencing. Here, we integrated fluorescence-activated nuclei sorting with long-read (PacBio Iso-Seq) and short-read (Illumina) RNA sequencing to generate isoform-resolved transcriptomes of five major cortical cell types, including MGE-derived GABAergic neurons, glutamatergic neurons, oligodendrocytes, astrocytes, and microglia from adult dorsolateral and orbitofrontal cortex. 

<img width="2340" height="2340" alt="Graphical_Abstract" src="https://github.com/user-attachments/assets/50a82d5d-646c-4900-a063-01445ca863fb" />

This repository contains the analysis code used to generate the cell type-resolved, full-length isoform atlas of the adult human cortex described in the manuscript above. Additionally, processed and summarized data can be found on [Zonedo](https://zenodo.org/records/20301473) and on the supplemental section of our pre-print. Lastly, cell type-resolved isoforms can be visualized on the [UCSC genome browser](https://genome.ucsc.edu/) using our [BED files](cell_type_resolved_isoforms.zip) and through our [Rshiny application](https://andyyang-isoseq.share.connect.posit.cloud/).

The repository is organized into four stages: (1) long-read processing and isoform reconstruction, (2) splice-junction validation and orthogonal isoform support, (3) differential expression and transcript usage analyses, and (4) disease-variant intersection and co-expression network analysis.

---

## Pipeline Overview

Briefly, raw long reads generated on the PacBio Sequel II platform were processsed with IsoSeq3 to generate high-quality full-length transcripts, which were aligned to the GRCh38 genome (GENCODE v38) using `minimap2`. Transcript models were collapsed and merged across samples within each cell type in the DLPFC and OFC using `TAMA` to generate nonredundant transcript annotations, followed by classification and QC filtering with `SQANTI3`. Transcript abundance was quantified using `Salmon`.

<img width="9921" height="7465" alt="PIPELINE_OVERVIEW" src="https://github.com/user-attachments/assets/39b7096d-ae02-4dd3-90a2-4036cda9442a" />


This repository contains a Snakemake workflow for processing PacBio Iso-Seq long-read transcriptome sequencing data from circular consensus sequence (CCS) reads through transcript annotation with SQANTI3. The workflow automates long-read preprocessing, genome alignment, transcript collapsing, and transcript classification using established long-read transcriptomics tools.

The pipeline is intended for reproducible analysis of PacBio Iso-Seq datasets and was developed for the analyses described in our pre-print.

---

## System requirements

**Operating systems tested:**
This workflow is designed for Linux (x86_64) systems.

Tested on:
- Rocky Linux 9.6

**Non-standard hardware**
The pipeline is intended for execution on Linux workstations or high-performance computing (HPC) clusters. No specialized hardware (e.g., GPUs) is required.

Recommended resources:

- ≥16 CPU cores
- 32–64 GB RAM

Typical memory usage:

- Most pipeline steps require <5 GB RAM
- Peak memory usage (~20–25 GB RAM) occurs during minimap2 alignment and TAMA transcript collapsing

---

## Software dependencies

The workflow uses Conda environments to manage software dependencies.

Major software components include:

- Snakemake
- IsoSeq3
- minimap2
- TAMA
- SQANTI3

Exact software versions are specified in the Conda environment files located in:

```
demo/envs/
```
---

## Prerequisites
The following software should already be installed and available in your environment before running the installation script:

- Git (v.2.46.0)
- Conda (v.26.1.1) (also: Miniconda or Anaconda)
- Snakemake (v.9.23.1)
- samtools (v.1.21)

The installation script installs the remaining workflow-specific software using Conda environments.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/AndyYaj/cortical-FANS-IsoSeq.git
cd cortical-FANS-IsoSeq
```

Install required software:

```bash
./install.sh
```

Download required reference files:

```bash
./download_references.sh
```

Typical installation time is approximately 30–60 minutes, depending on network performance and Conda environment creation.

---

## Demo
A complete end-to-end demonstration using a small example dataset is provided in:

```
demo/
```

See:

```
demo/README.md
```

for detailed instructions.

Demo input data are hosted separately on Zenodo due to GitHub file size limitations. The `demo/download_demo_data.sh` script retrieves the required files automatically.

---

## Running the pipeline on your own data

The Snakemake workflow provided in `demo/` is intended solely as a minimal reproducible example using a reduced dataset.

For analysis of full PacBio Iso-Seq datasets, the repository includes the complete analysis scripts in:

```
code/
```

These scripts implement the full workflow used for the manuscript, including processing multiple samples, donors, and cell types. Users should adapt these scripts and configuration files to their own experimental design and input data.

The demo workflow serves as a reference implementation for verifying software installation and demonstrating the major processing steps.

---

## Code

All of the code used to produce the main figures in the manuscript can be found in [`code/main_figures/`](https://github.com/AndyYaj/cortical-FANS-IsoSeq/tree/main/code/main_figures).  All of the data used in those analyses can be found on [Zonedo](https://zenodo.org/records/20301473).

---

## Data accessibility

| Data | Description | Data Processing | 
|------|-------------| ----------------------------|
| [IsoSeq Long Reads](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE330753) | Raw DLPFC and OFC IsoSeq long read (BAM) files. | [scripts](https://github.com/AndyYaj/cortical-FANS-IsoSeq/tree/main/code/data_processing) |
| [DLPFC Short Reads](https://www.synapse.org/#!Synapse:syn12034263) | Raw short-read RNA-seq data (cohort 1) — Human DLPFC | [scritps](https://github.com/CommonMindConsortium/RAPiD-nf) |
| OFC Short Reads | Raw short-read RNA-seq data (cohort 2) — Human OFC (GEO accession ID in progress) | [scripts](https://github.com/CommonMindConsortium/RAPiD-nf) |
| [FANS Short Reads](https://www.synapse.org/Synapse:syn25716684/wiki/610496) | Raw short-read RNA-seq data (cohort 3) — Human BM10, BM17, BM22, BM36, BM44. | [scripts](https://github.com/CommonMindConsortium/RAPiD-nf) |

---

## License

This repository is released under the MIT License. See `LICENSE` for full terms. This license is approved by the Open Source Initiative.

---

## Contact

Questions: Michael S. Breen (michael.breen@mssm.edu).

For issues with the code or pipeline, please open a GitHub issue.

---

## Citation

If you use this code or the associated dataset, please cite:

Yang A et al. Long-read transcriptomics of purified human cortical cell types exposes glial isoform complexity and disease-relevant transcript architecture. *bioRxiv* (2025). https://doi.org/10.1101/2025.11.25.690524
