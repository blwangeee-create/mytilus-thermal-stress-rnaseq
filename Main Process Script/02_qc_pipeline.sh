#!/bin/bash
set -euo pipefail

# ============================================================
# Phase 2: Quality Control Pipeline
# ============================================================

RAW_DIR="/Volumes/Expansion/NCBI/PRJNA934294/raw_data"
QC_DIR="/Volumes/Expansion/NCBI/PRJNA934294/qc"
TRIMMED_DIR="$RAW_DIR/trimmed"
THREADS=8

# Create directories
mkdir -p "$QC_DIR/fastqc_raw" "$QC_DIR/fastqc_trimmed" "$QC_DIR/multiqc" "$TRIMMED_DIR"

# Sample list (excluding SRR23438726 which is already done)
SAMPLES=(
    SRR23438720 SRR23438721 SRR23438722
    SRR23438723 SRR23438724 SRR23438725
    SRR23438727 SRR23438728
)

echo "=== Step 1: FastQC on raw reads ==="
for srr in "${SAMPLES[@]}"; do
    echo "Running FastQC on $srr..."
    fastqc -t 2 -o "$QC_DIR/fastqc_raw/" \
        "$RAW_DIR/${srr}_1.fastq" "$RAW_DIR/${srr}_2.fastq" &
done
wait
echo "FastQC raw complete."

echo "=== Step 2: MultiQC on raw FastQC reports ==="
multiqc -o "$QC_DIR/multiqc/raw" -n "multiqc_raw" "$QC_DIR/fastqc_raw/" --force

echo "=== Step 3: Trim Galore (adapter + quality trimming) ==="
for srr in "${SAMPLES[@]}"; do
    echo "Trimming $srr..."
    trim_galore --paired \
        --fastqc --fastqc_args "-t 2 -o $QC_DIR/fastqc_trimmed/" \
        -j 2 \
        -o "$TRIMMED_DIR/" \
        "$RAW_DIR/${srr}_1.fastq" "$RAW_DIR/${srr}_2.fastq" &
done
wait
echo "Trimming complete."

echo "=== Step 4: MultiQC on trimmed FastQC reports ==="
multiqc -o "$QC_DIR/multiqc/trimmed" -n "multiqc_trimmed" "$QC_DIR/fastqc_trimmed/" --force

echo "=== QC Pipeline Complete ==="
echo "Results in: $QC_DIR/"
