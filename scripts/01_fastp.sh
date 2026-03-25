#!/bin/bash
# 01_fastp.sh
# Quality control and adapter trimming of raw paired-end reads using fastp
# Tool: fastp v1.3.0
# Usage: bash scripts/01_fastp.sh
# Run from: ~/BINF_DESKTOP/BINF6110/assignment03/

# Define directories
DATA_DIR="data"
QC_DIR="results/fastp"
mkdir -p $QC_DIR

# Sample IDs
SAMPLES="SRR8146968 SRR8146978 SRR8146982 SRR8146938 SRR8146935 SRR8146936"

for SRR in $SAMPLES; do
    echo "Running fastp on $SRR..."

    fastp \
        --in1 ${DATA_DIR}/${SRR}_1.fastq.gz \
        --in2 ${DATA_DIR}/${SRR}_2.fastq.gz \
        --out1 ${QC_DIR}/${SRR}_1_trimmed.fastq.gz \
        --out2 ${QC_DIR}/${SRR}_2_trimmed.fastq.gz \
        --html ${QC_DIR}/${SRR}_fastp.html \
        --json ${QC_DIR}/${SRR}_fastp.json \
        --thread 4 \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --correction \
        2> ${QC_DIR}/${SRR}_fastp.log

    echo "Done: $SRR"
done

echo "All samples processed. QC reports in ${QC_DIR}/"
