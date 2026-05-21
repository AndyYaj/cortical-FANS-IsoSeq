# Creating the IsoSeq3 Environment

```bash
conda create -n IsoSeq3 python=3.7 biopython samtools isoseq3 pbccs lima pbcoretools bamtools pysam pbbam bam2fastx minimap2 -c bioconda -c conda-forge -y

conda activate IsoSeq3
```

# Iso-Seq Processing

## 1. Circular Consensus Sequence (CCS) Calling

In this stage, multiple subreads generated from the same DNA molecular are combined to produce a highly accurate consensus read. This improves sequencing accuracy and creates high-quality reads suitable for downstream transcriptome analysis. 

```bash
# Load required modules
ml anaconda3/2020.11
ml python/3.7.3

# Activate IsoSeq3 environment
source activate IsoSeq3

PROJECT_DIR="/path/to/IsoSeq"

# Run CCS calling
ccs \
    --min-rq 0.9 \
    --min-passes 3 \
    ${PROJECT_DIR}/${file_name1}/${file_name1}.subreads.bam \
    ${PROJECT_DIR}/${file_name1}/${file_name1}.ccs.${chunk_num}.bam \
    --chunk ${chunk_num}/20 \
    -j 0

echo "Finished running ${file_name1}"
```

**Notes**
- ```--min-rq 0.9``` sets the minimum read quality threshold
- ```--min-passes 3``` requires at least three passes for CCS generation.
- ```--chunk ${chunk_num}/20``` processes one chunk out of 20 total chunks.
- ```-j 0``` uses all available CPU cores.
- ```${file_name1}``` and ```${chunk_num}``` should be defined before running the script.


After CCS calling, merge CCS BAMS for each sample. This is only necessary because the chunk option was used for parallel processing of the subread bam files. 

```bash
for i in *.fofn
    do
    name=`basename $i .fofn`
    pbmerge -o ${name}.ccs.bam ${name}.fofn
done
```

## 2. FLNC

Following CCS generation, full-length non-concatemer (FLNC) reads are identified by removing sequencing primers and detecting complete transcript molecules. This step uses ```lima```, a PacBio demultiplexing and primer-removal tool, to identify reads containing both 5' and 3' Iso-Seq primers. The resulting FLNC reads represent high-confidence full-length transcript sequences.

```bash
PROJECT_DIR="/path/to/IsoSeq"

# Remove Iso-Seq primers and identify FLNC reads
lima \
    --isoseq \
    --different \
    --min-passes 1 \
    --split-bam-named \
    --dump-clips \
    --dump-removed \
    -j 0 \
    ${PROJECT_DIR}/${file_name1}/${file_name1}.ccs.bam
    ${PROJECT_DIR}/isoseq_primers.fa
    ${PROJECT_DIR}/${file_name1}/${file_name1}.fl.bam

echo "Finished running ${file_name1}"
```

After identifying FLNC reads, the next step is transcript refinement using ```isoseq3 refine```. This process trims residual primer sequences, removes artifical concatemers, and filters reads based on the prescence of poly(A) tails.

```bash
# Define project directory
PROJECT_DIR="/path/to/IsoSeq"

# Define sample list (In this example, I'm using samples from donor #276)
samples=(
    276_GABA
    276_GLU
    276_OLIG
)

# Loop through samples
for sample in "${samples[@]}"; do

    echo "Running refinement for ${sample}"

    isoseq3 refine \
        --require-polya \
        ${PROJECT_DIR}/${sample}/${sample}.fl.NEB_5p--NEB_3p.bam \
        ${PROJECT_DIR}/isoseq_primers.fa \
        ${PROJECT_DIR}/${sample}/${sample}.flnc.bam

    echo "Finished running ${sample}"

done
```

## 3. Quality Control Assessment

Following FLNC refinement, QC analysis is performed to evaluate read quality and sequencing characteristics. In this workflow, BAM files are first converted to FASTQ format using ```bam2fastq```, followed by quality assessment with ```FastQC```. 

```bash
module load fastqc

PROJECT_DIR="/path/to/IsoSeq"
mkdir -p ${PROJECT_DIR}/fastqc

samples=(
    276_GABA
    276_GLU
    276_OLIG
)

# Loop through samples
for sample in "${samples[@]}"; do

    echo "Running QC for ${sample}"

    # Convert CCS BAM to FASTQ
    bam2fastq \
        -o ${PROJECT_DIR}/${sample}/${sample}.ccs \
        ${PROJECT_DIR}/${sample}/${sample}.ccs.bam

    # Convert FLNC BAM to FASTQ
    bam2fastq \
        -o ${PROJECT_DIR}/${sample}/${sample}.flnc \
        ${PROJECT_DIR}/${sample}/${sample}.flnc.bam

    # Run FastQC on FLNC reads
    fastqc \
        ${PROJECT_DIR}/${sample}/${sample}.flnc.fastq.gz \
        -o ${PROJECT_DIR}/fastqc

    # Run FastQC on CCS reads
    fastqc \
        ${PROJECT_DIR}/${sample}/${sample}.ccs.fastq.gz \
        -o ${PROJECT_DIR}/fastqc

    echo "Finished QC for ${sample}"

done
```

## 4. Genome Alignment with Minimap2

After FLNC reads are generated and quality-controlled, transcripts are aligned to the reference genome using ```minimap2```. The ```splice``` alignment preset is optimized for long-read RNA sequencing data and enables accurate mapping across exon-intron junctions. Secondary alignments are suppressed to retain primary transcript mappings only. 

The workflow below converts FLNC BAM files to FASTQ format and aligns reads to the GRCh38 reference genome.

```bash
# Load required modules
module load minimap2

# Define project directory
PROJECT_DIR="/path/to/IsoSeq"

# Define reference genome
REFERENCE="${PROJECT_DIR}/GRCh38.primary_assembly.genome.fa"

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

    echo "Running alignment for ${sample}"

    # Convert FLNC BAM to FASTQ
    bam2fastq \
        -o ${PROJECT_DIR}/${sample}/${sample}.flnc \
        ${PROJECT_DIR}/${sample}/${sample}.flnc.bam

    # Align FLNC reads to reference genome
    minimap2 \
        -ax splice \
        -uf \
        --secondary=no \
        -C5 \
        -t 18 \
        --MD \
        ${REFERENCE} \
        ${PROJECT_DIR}/${sample}/${sample}.flnc.fastq.gz \
        > ${PROJECT_DIR}/${sample}/${sample}.flnc.sam

    echo "Finished alignment for ${sample}"

done
```

Aligned FLNC reads are next sorted using `samtools`. Coordinate-sorted SAM files are required for downstream transcript annotation and collapsing steps performed with the TAMA pipeline. Sorting ensures that alignments are ordered consistently by genomic position, which improves compatibility with downstream tools and indexing workflows.

```bash

module load samtools
PROJECT_DIR="/path/to/IsoSeq"

samples=(
    276_GABA
    276_GLU
    276_OLIG
)

for sample in "${samples[@]}"; do

    echo "Sorting SAM file for ${sample}"

    samtools sort \
        ${PROJECT_DIR}/${sample}/${sample}.flnc.sam \
        -o ${PROJECT_DIR}/${sample}/${sample}.flnc.sorted.sam

    echo "Finished sorting ${sample}"

done
```
