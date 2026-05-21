# SQANTI3

SQANTI3 is used for structural and quality annotation of long-read transcript models. It provides extensive classification of isoforms (e.g., novel, known, artifact) and is a key downstream step after transcriptome construction with TAMA.

Because SQANTI3 depends on multiple bioinformatics tools and specific Python versions, installation is performed in a dedicated conda environment with additional manual dependencies.


---

## 1. Install SQANTI3

```bash
# Download SQANTI3 source code
wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v4.2.tar.gz

# Extract archive
tar -xvf v4.2.tar.gz
cd SQANTI3-4.2

# Create conda environment from YAML file
conda env create -f SQANTI3.conda_env.yml

# Activate environment
source activate SQANTI3.env
```

## 2. Install UCSC utility: gtfToGenePred
This tool is required for GTF processing but may fail during conda installation, so it's installed manually:

```bash
# Download binary
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred \
    -P SQANTI3-4.2/utilities/

# Make executable
chmod +x SQANTI3-4.2/utilities/gtfToGenePred
```

## 3. Install cDNA_Cupcake (required dependency)
```bash
# Download Cupcake toolkit
wget https://github.com/Magdoll/cDNA_Cupcake/archive/refs/tags/v28.0.0.tar.gz

# Extract
tar -xvf v28.0.0.tar.gz
cd cDNA_Cupcake

# IMPORTANT: downgrade Python for compatibility
conda install python=3.7

# Install package
python setup.py build
python setup.py install
```

**Note**
- SQANTI3 requires a dedicated conda enviornment due to strict dependency constraints
- Manual installation of `gtfToGenePred` is often necessary due to conda compatability issues
- *Python 3.7* is required for `cDNA_Cupcake` due to NumPy compatability errors in newer version
- This environment is used for transcript classification and QC of TAMA-derived isoforms.

## 4. Preparing Data for SQANTI3 Analysis
SQANTI3 strongly prefers transcript annotations in **GTF format**, even though TAMA outputs results in BED format. To ensure compatibility, the TAMA-provided conversion script is used to generate a GTF file directly from collapsed BED annotations.

However, additional filtering is required because SQANTI3 can fail when encountering non-standard contigs (e.g., `GL000008.2`, `KI270724.1`). These are removed to retain only primary chromosome assemblies.

```bash
# Activate TAMA environment
source activate TAMA

PROJECT_DIR="/path/to/IsoSeq"

# Convert BED → GTF
python ${PROJECT_DIR}/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py \
    GABAmerged.bed \
    GABAmerged.gtf

# Filter to primary chromosomes only
grep 'chr' GABAmerged.gtf > GABAmerged_chr.gtf
```

**Notes**
- SQANTI3 may fail on alternative contigs or unplaced scaffolds, so we restrict the annotation to standard chromosomes (hence why we filter to primary chromosomes only).
- Filtering to `chr*` entries removes alternative contigs and scaffolds that can break SQANTI3 parsing.
- The resulting `*_chr.gtf` file is the recommended input for SQANTI3 analysis.

## 5. SQANTI Quality Control
After generating a high-confidence transcriptome with TAMA, SQANTI3 is used to perform structural and functional quality control of isoforms. This step classifies transcripts relative to known annotations (e.g., known, novel, fusion, antisense), evaluates splice junction support, and integrates external evidence such as short-read RNA-seq, CAGE peaks, and poly(A) motifs.

```bash
# Load environment
module load anaconda3/2020.11
source activate SQANTI3.env

PROJECT_DIR="/path/to/IsoSeq"

# Add cDNA_Cupcake to PYTHONPATH
export PYTHONPATH=$PYTHONPATH:${PROJECT_DIR}/SQANTI3-4.2/cDNA_Cupcake-28.0.0
export PYTHONPATH=$PYTHONPATH:${PROJECT_DIR}/SQANTI3-4.2/cDNA_Cupcake-28.0.0/sequence

# Run SQANTI3 QC for GABA
python /sc/arion/projects/breen_lab/IsoSeq/SQANTI3-4.2/sqanti3_qc.py \
    ${PROJECT_DIR}/TAMAfiltering/GABA.filt.chr.gtf \
    ${PROJECT_DIR}/gencode.v38.annotation.gtf \
    ${PROJECT_DIR}/GRCh38.primary_assembly.genome.fa \
    -o GABA \
    -d ${PROJECT_DIR}/TAMAfiltering/GABA_allopts_SQANTI \
    -fl ${PROJECT_DIR}/TAMAfiltering/SQANTI_optional_inputs/GABA.FLcount \
    --skipORF \
    --polyA_motif_list ${PROJECT_DIR}/TAMAfiltering/SQANTI_optional_inputs/polyAmotifs.txt \
    --cage_peak ${PROJECT_DIR}/TAMAfiltering/SQANTI_optional_inputs/refTSS_v3.3_human_coordinate.hg38.bed \
    --coverage "${PROJECT_DIR}/STAR/GABA/GABASJ/" \
    --SR_bam "${PROJECT_DIR}/STAR/GABA/GABABAM/"
```
**Notes**
- SQANTI3 classifies isoforms relative to the reference annotation (e.g., GENCODE v38).
- Short-read support (--coverage, --SR_bam) improves splice junction validation.
- CAGE peaks help validate transcription start sites (TSS).
- Poly(A) motif files improve 3′ end validation.
- This analysis is performed separately for each cell type in the DLPFC and OFC.


## 6. SQANTI3 Refinement and Filtering
After SQANTI3 QC, the final step is transcript refinement using SQANTI3 filtering rules. This stage applies evidence-based filters (e.g., splice junction support, structural categories, and classification rules) to generate a high-confidence, refined transcriptome.

The output includes corrected transcript models in both FASTA and GTF formats, along with a refined classification report.

```bash
# Load environment
module load anaconda3/2020.11
source activate SQANTI3.env

PROJECT_DIR="/path/to/IsoSeq"

# Add cDNA_Cupcake to PYTHONPATH
export PYTHONPATH=$PYTHONPATH${PROJECT_DIR}/SQANTI3-4.2/cDNA_Cupcake-28.0.0
export PYTHONPATH=$PYTHONPATH:${PROJECT_DIR}/SQANTI3-4.2/cDNA_Cupcake-28.0.0/sequence

# Run SQANTI3 filtering (refinement step)
python /sc/arion/projects/breen_lab/IsoSeq/SQANTI3-4.2/sqanti3_RulesFilter.py \
    ${PROJECT_DIR}/TAMAfiltering/GABA_allopts_SQANTI/GABA_classification.txt \
    ${PROJECT_DIR}/TAMAfiltering/GABA_allopts_SQANTI/GABA_corrected.fasta \
    ${PROJECT_DIR}/TAMAfiltering/GABA_allopts_SQANTI/GABA_corrected.gtf \
    --saturation \
    --report both
```
**Notes**
- This step applies SQANTI3’s rule-based filtering system to remove low-confidence or artifact transcripts.
- The `--saturation` flag evaluates whether transcript discovery is approaching completeness.
- Output includes:
    - Refined transcript FASTA
    - Refined transcript GTF
    - Updated classification report
