# IsoSeq-FANS-Human-Cortex

Code for: **Long-read transcriptomics of purified human cortical cell types exposes glial isoform complexity and disease-relevant transcript architecture**

Authors: Yang A, Rodriguez de los Santos M, Kozlenkov A, Vadukapuram R, Hurd Y, Dracheva S, Humphrey J, Breen MS.

Affiliation: Icahn School of Medicine at Mount Sinai, New York, NY 

**Preprint**: https://www.biorxiv.org/content/10.1101/2025.11.25.690524v2
**Interactive atlas**: https://andyyang.shinyapps.io/Cortical_FANS_IsoSeq/
**Archived release**: https://doi.org/10.5281/zenodo.20301473

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

---

## System requirements

**Operating systems tested:**
- Linux: Ubuntu 20.04 LTS, CentOS 7
- macOS: Sonoma 14.x (analysis scripts only; long-read processing tools are Linux-only)

**Software dependencies:**

Long-read processing:
- PacBio SMRT Link v11.0
- IsoSeq v3.8
- minimap2 v2.24
- TAMA (release 2021-12-08)
- SQANTI3 v5.1.2
- GffCompare v0.12.6

Short-read processing:
- STAR v2.7.10
- Salmon v1.10.0
- Subread/featureCounts v2.0.3
- bcl2fastq2 v2.20

Analysis (R ≥ 4.3.0):
- DESeq2 v1.38, edgeR v3.40, limma v3.54
- DRIMSeq v1.26, stageR v1.20
- WGCNA v1.72
- ComplexHeatmap v2.14, pheatmap v1.0.12, ggplot2 v3.4.4

Analysis (Python ≥ 3.10):
- SUPPA2 v2.3
- CPAT v3.0.4, CPC2 v2.0, GeneMarkS-T v5.1
- ViennaRNA v2.5.1
- pandas, numpy, matplotlib v3.7, seaborn v0.12

**Non-standard hardware:**
None. All analyses run on standard Linux compute infrastructure. Long-read alignment and isoform reconstruction benefit from ≥64 GB RAM and ≥16 CPU cores.

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
