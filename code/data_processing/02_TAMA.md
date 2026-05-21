[02_TAMA.md](https://github.com/user-attachments/files/28103033/02_TAMA.md)
# TAMA

## Installing TAMA
```git clone https://github.com/GenomeRIK/tama```

There is one erorr that needs to be corrected in the tama_collapse.py script. Copy and paste the code below in palce of al the "import lines":

```bash
import re
import sys
import time
from Bio import SeqIO
try:
    import StringIO
except ImportError:
    from io import StringIO
from Bio import AlignIO
import os
import argparse
```

## Setting up TAMA Environment

```bash
condact
conda create -n TAMA python=2.7.17
conda activate TAMA
conda install -n TAMA biopython
```

## 1. TAMA Collapse
Following genome alignment and coordinate sorting, transcript models are collapsed using the TAMA (Transcriptome Annotation by Modular Algorithms) pipeline. TAMA collapse merges redundant transcript alignments into non-redundant isoform models while preserving transcript structure information.

This step reduces transcript redundancy and generates high-confidence transcript annotations for downstream isoform analysis. Additionally, TAMA collapse provides the following inofmration: 
- Source information for each predicted feature
- Variation calling
Genomic poly-A detection
- Strand ambiguity

```bash
module load anaconda2/2019.10
source activate TAMA

# Define project directory
PROJECT_DIR="/path/to/IsoSeq"

# Define reference genome
REFERENCE="${PROJECT_DIR}/GRCh38.primary_assembly.genome.fa"

# Define TAMA installation directory
TAMA_DIR="${PROJECT_DIR}/tama"

# Define sample list
samples=(
    276_GABA
    276_GLU
    276_OLIG
    344_GABA
    344_OLIG
)

# Loop through samples
for sample in "${samples[@]}"; do

    echo "Running TAMA collapse for ${sample}"

    # Collapse redundant transcript models
    python ${TAMA_DIR}/tama_collapse.py \
        -s ${PROJECT_DIR}/${sample}/${sample}.flnc.sorted.sam \
        -f ${REFERENCE} \
        -p ${PROJECT_DIR}/${sample}/${sample} \
        -x no_cap \
        -i 95 \
        -z 1000

    echo "Finished TAMA collapse for ${sample}"

done
```

## 2. Transcript merging with TAMA Merge

After generating collapsed transcript models for each sample, the final step is to merge isoforms across samples to build a unified transcriptome annotation. `tama_merge.py` combines transcript models from multiple datasets while retaining consistent isoform structure definitions and reducing redundancy across conditions.

In this step, a file list (`*.txt`) is used to specify all TAMA collapse outputs that should be merged into a single annotation set. The `filelist.txt` sould contain paths to TAMA collapse outputs. Here samples, are grouped by cell type prior to merging:

```
# GABAfilelist.txt 
276_GABA.sam
344_GABA.sam
365_GABA.sam

# GLUfilelist.txt 
276_GLU.sam
344_GLU.sam
365_GLU.sam
```

Each group is merged separately so that cell type-specific transcriptomes remain distinct.

```bash
PROJECT_DIR="/path/to/IsoSeq"
OUT_DIR="${PROJECT_DIR}/TAMAoutputs"

# Run TAMA merge for GABA group
python ${PROJECT_DIR}/tama/tama_merge.py \
    -f ${PROJECT_DIR}/GABAfilelist.txt \
    -p ${OUT_DIR}/GABAmerged \
    -z 1000 \
    -a 1000 \
    -d merge_dup
```

**Note**
- Samples are merged **by cell type (within each brain region), not across all samples**
- Each `filelist.txt` should contain only samples from a single cell type (e.g., GABA datasets are processed independently to preserve biological specificity).
- `-z 1000` and `-a 1000` control merging distance thresholds
- `-d merge_dup` merges duplicate transcript models preserving isoform structure. 

## 3. Transcript Filtering and Read Support Analysis with TAMA GO

After merging transcriptomes by cell type, additional filtering is performed using the TAMA GO toolkit to improve transcript model quality and remove low-confidence or artifactual isoforms. This stage integrates read support information, removes singleton-supported models, filters poly(A)-related artifacts, and eliminates fragmented transcript models. 

**Notes**
- This workflow is performed separately for each cell type.
- Each filtering stage progressively removes low-confidence transcript models:
    - Singleton-supported transcripts
    - Poly(A)-associated artifacts
    - Fragmented models

### 3a. Read Support Assignment

First, read support levels are calculated for merged transcripts using alignment and supporting read information.

```bash
# Load environment
module load anaconda2/2019.10
source activate TAMA

# Define project directory
PROJECT_DIR="/path/to/IsoSeq"
TAMA_DIR="${PROJECT_DIR}/tama"

# Define sample group
GROUP="GABA"

# Compute read support levels
python ${TAMA_DIR}/tama_go/read_support/tama_read_support_levels.py \
    -f ${PROJECT_DIR}/TAMAfiltering/TransReadFiles/${GROUP}files.txt \
    -m ${PROJECT_DIR}/TAMAfiltering/mergedBEDfiles/${GROUP}merged_merge.txt \
    -o ${GROUP} \
    -mt tama
```

### 3b. Remove Singleton Transcript Models
```bash
python ${TAMA_DIR}/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py \
    -b ${PROJECT_DIR}/TAMAfiltering/mergedBEDfiles/${GROUP}merged.bed \
    -r ${PROJECT_DIR}/TAMAfiltering/${GROUP}_read_support.txt \
    -o ${GROUP}1 \
    -l transcript \
    -k remove_multi
```

### 3c. Remove Poly(A)-Related Artifacts
```bash
python ${TAMA_DIR}/tama_go/filter_transcript_models/tama_remove_polya_models_levels.py \
    -b ${PROJECT_DIR}/TAMAfiltering/mergedBEDfiles/${GROUP}merged.bed \
    -f ${PROJECT_DIR}/TAMAfiltering/polyAfiles/${GROUP}files.txt \
    -r ${PROJECT_DIR}/TAMAfiltering/${GROUP}_read_support.txt \
    -o ${GROUP}2 \
    -l transcript
```

### 3d. Extract Original Transcript IDs
```bash
cut -f2 ${GROUP}_read_support.txt | sed '1d' > ${GROUP}_original_ids
grep -v 'removed' ${GROUP}1_singleton_report.txt | cut -f1,2 | sed '1d' | sed 's/\t/;/g' > ${GROUP}1_original_ids
```

### 3e. Restore Transcript ID Mapping
```bash
awk 'NR==FNR{a[NR]=$0;next}{$4=a[FNR]}1' ${GROUP}1_original_ids ${GROUP}1.bed | sed 's/ /\t/g' > ${GROUP}1.origIDs.bed
```

### 3f. Remove Fragmented models
```bash
python ${TAMA_DIR}/tama_go/filter_transcript_models/tama_remove_fragment_models.py \
    -f ${GROUP}1.origIDs.bed \
    -o ${GROUP}3
```

### 3g. Final Filtering Steps
```bash
cut -f4 ${GROUP}3.bed | sed 's/^[^;]*;//' > ${GROUP}3_ids
grep -w -F -f ${GROUP}3_ids ${GROUP}2.bed > ${GROUP}.filt.bed
grep 'chr' ${GROUP}.filt.bed > ${GROUP}.filt.chr.bed
cut -f4 ${GROUP}.filt.chr.bed | sed 's/^[^;]*;//' > ${GROUP}_final_ids
grep -w -F -f ${GROUP}_final_ids ${GROUP}_read_support.txt > ${GROUP}.filt_read_support.txt
```
