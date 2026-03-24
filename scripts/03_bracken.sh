#!/bin/bash
# 03_bracken.sh
# Species-level abundance re-estimation from Kraken2 reports using Bracken
# Tool: Bracken v3.0.1
# Database: k2_standard_08GB (built October 2025)
# Read length: 300 bp (NextSeq 500, matches database300mers.kmer_distrib)
# Usage: bash scripts/03_bracken.sh
# Run from: ~/BINF_DESKTOP/BINF6110/assignment03/

# Define directories
KRAKEN_DIR="results/kraken2"
BRACKEN_DIR="results/bracken"
DB="$HOME/kraken2"
READ_LEN=300
LEVEL="S"        # Species level
THRESHOLD=10     # Minimum reads required before re-estimation
mkdir -p $BRACKEN_DIR

# Sample IDs
SAMPLES="SRR8146968 SRR8146978 SRR8146982 SRR8146938 SRR8146935 SRR8146936"

for SRR in $SAMPLES; do
    echo "Running Bracken on $SRR..."

    bracken \
        -d $DB \
        -i ${KRAKEN_DIR}/${SRR}.report \
        -o ${BRACKEN_DIR}/${SRR}.bracken \
        -w ${BRACKEN_DIR}/${SRR}_bracken.report \
        -r $READ_LEN \
        -l $LEVEL \
        -t $THRESHOLD

    echo "Done: $SRR"
done

echo "All samples processed. Combining Bracken outputs..."

# Combine all sample outputs into a single abundance table
python ~/BINF_DESKTOP/BINF6110/lectures/lecture12/Bracken-3.1/analysis_scripts/combine_bracken_outputs.py \
    --files $(ls ${BRACKEN_DIR}/*.bracken | tr '\n' ' ') \
    --names SRR8146968,SRR8146978,SRR8146982,SRR8146938,SRR8146935,SRR8146936 \
    --output ${BRACKEN_DIR}/combined_abundance.txt

echo "Combined abundance table written to ${BRACKEN_DIR}/combined_abundance.txt"
