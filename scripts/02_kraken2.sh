#!/bin/bash
# 02_kraken2.sh
# Taxonomic classification of trimmed paired-end reads using Kraken2
# Tool: Kraken2 v2.17.1
# Database: k2_standard_08GB (built October 2025)
# Usage: bash scripts/02_kraken2.sh
# Run from: ~/BINF_DESKTOP/BINF6110/assignment03/

# Define directories
QC_DIR="results/fastp"
KRAKEN_DIR="results/kraken2"
DB="$HOME/kraken2"
mkdir -p $KRAKEN_DIR

# Sample IDs
SAMPLES="SRR8146968 SRR8146978 SRR8146982 SRR8146938 SRR8146935 SRR8146936"

for SRR in $SAMPLES; do
    echo "Running Kraken2 on $SRR..."

    kraken2 \
        --db $DB \
        --paired \
        --threads 4 \
        --report ${KRAKEN_DIR}/${SRR}.report \
        --output ${KRAKEN_DIR}/${SRR}.kraken \
        ${QC_DIR}/${SRR}_1_trimmed.fastq.gz \
        ${QC_DIR}/${SRR}_2_trimmed.fastq.gz

    echo "Done: $SRR"
done

echo "All samples classified. Results in ${KRAKEN_DIR}/"
