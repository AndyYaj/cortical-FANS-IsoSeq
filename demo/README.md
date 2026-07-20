# Demo

## Overview

This demo executes the complete long-read processing workflow on a reduced PacBio Iso-Seq dataset (~100 MB). The dataset consists of a subset of CCS reads from a single donor and cell type and is intended to verify successful installation and demonstrate end-to-end execution of the pipeline.

The workflow performs:

1. IsoSeq3 preprocessing
2. Genome alignment with minimap2
3. Transcript collapsing with TAMA
4. Transcript annotation using SQANTI3

The demo reproduces the complete workflow using a small dataset while requiring substantially less compute time than a full analysis.

---

## Requirements

- Linux (tested on Rocky Linux 9.6)
- Conda (Miniconda or Anaconda)
- Git

---

## Setup

Clone the repository:

```bash
git clone https://github.com/AndyYaj/cortical-FANS-IsoSeq.git
cd cortical-FANS-IsoSeq/demo
```

Install software dependencies:

```bash
./install.sh
```

Download required reference files:

```bash
./download_references.sh
```

## Download demo dataset

The demo requires the following PacBio subreads files:

```
demo/data/
├── 276_GABA_demo.subreads.bam
└── 276_GABA_demo.subreads.bam.pbi
```

These files are hosted on Zenodo and can be downloaded using:

```bash
./download_demo_data.sh
```

The download script retrieves the required files and places them in the correct directory structure.

---

## Running the demo

Execute the workflow:

```bash
snakemake --use-conda --cores 8 -p
```

---

## Expected runtime

Approximately **1 hour** using:

- 8 CPU cores
- Up to 25 GB RAM

Most workflow steps require less than 5 GB RAM. Peak memory usage occurs during minimap2 alignment and TAMA transcript collapsing.

---

## Expected outputs

| File | Description |
|------|-------------|
| `Demo_276_GABA_RulesFilter_Classification.txt` | SQANTI3 structural classification results |
| `Demo_276_GABA_RulesFilter_Filtered.gtf` | SQANTI3 reconstructed isoforms |

Reference outputs are provided in:

```
expected_outputs/
```

Successful execution should reproduce these outputs (minor differences in timestamps or log files are expected).
