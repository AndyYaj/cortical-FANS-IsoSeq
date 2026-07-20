#!/usr/bin/env bash
set -euo pipefail

REF_DIR="references"
mkdir -p "$REF_DIR"
mkdir -p logs

echo "===================================="
echo "Downloading reference files"
echo "===================================="


# ---------------------------------------------------------
# 1. GENCODE v38 GTF
# ---------------------------------------------------------
if [[ ! -f "${REF_DIR}/gencode.v38.annotation.gtf" ]]; then
    echo "Downloading GENCODE v38 GTF..."

    wget -O "${REF_DIR}/gencode.v38.annotation.gtf.gz" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz

    gunzip -f "${REF_DIR}/gencode.v38.annotation.gtf.gz"
else
    echo "[SKIP] GTF already exists"
fi


# ---------------------------------------------------------
# 2. GRCh38 FASTA
# ---------------------------------------------------------
if [[ ! -f "${REF_DIR}/GRCh38.primary_assembly.genome.fa" ]]; then
    echo "Downloading GRCh38 genome..."

    wget -O "${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz" \
    https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz

    gunzip -f "${REF_DIR}/GRCh38.primary_assembly.genome.fa.gz"
else
    echo "[SKIP] FASTA already exists"
fi


# ---------------------------------------------------------
# 3. IsoSeq primers
# ---------------------------------------------------------
if [[ ! -f "${REF_DIR}/isoseq_primers.fa" ]]; then
    echo "Creating IsoSeq primer file..."

    cat > "${REF_DIR}/isoseq_primers.fa" <<EOF
>NEB_5p
GCAATGAAGTCGCAGGGTTGGG
>NEB_3p
GTACTCTGCGTTGATACCACTGCTT
EOF
else
    echo "[SKIP] primers exist"
fi


echo "===================================="
echo "REFERENCES READY"
echo "===================================="
