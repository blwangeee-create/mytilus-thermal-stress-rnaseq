#!/bin/bash
set -euo pipefail

SRA_BIN="/Volumes/Expansion/NCBI/PRJNA934294/sratoolkit.3.4.1-mac-arm64/bin"
RAW_DIR="/Volumes/Expansion/NCBI/PRJNA934294/raw_data"

SRRS=(
    SRR23438720 SRR23438721 SRR23438722
    SRR23438723 SRR23438724 SRR23438725
    SRR23438727 SRR23438728
)

echo "=== Step 1: Download SRA files ==="
for srr in "${SRRS[@]}"; do
    if [ -f "$RAW_DIR/${srr}/${srr}.sra" ]; then
        echo "$srr already downloaded, skipping..."
    else
        echo "Downloading $srr..."
        $SRA_BIN/prefetch $srr --max-size 50G -O "$RAW_DIR/"
    fi
done

echo "=== Step 2: Convert SRA to FASTQ ==="
for srr in "${SRRS[@]}"; do
    if [ -f "$RAW_DIR/${srr}_1.fastq" ] && [ -f "$RAW_DIR/${srr}_2.fastq" ]; then
        echo "$srr already converted, skipping..."
    else
        echo "Converting $srr to FASTQ..."
        $SRA_BIN/fasterq-dump "$RAW_DIR/${srr}/${srr}.sra" \
            --split-files --threads 4 -O "$RAW_DIR/"
    fi
done

echo "=== Download and conversion complete ==="
ls -lh "$RAW_DIR/"*.fastq | awk '{print $NF, $5}'
