#!/bin/bash
set -euo pipefail

# ============================================================
# Phase 4: STAR Alignment
# ============================================================

INDEX_DIR="/Volumes/Expansion/NCBI/PRJNA934294/index/STAR_MCOR1.1"
TRIMMED_DIR="/Volumes/Expansion/NCBI/PRJNA934294/qc/fastp"
BAM_DIR="/Volumes/Expansion/NCBI/PRJNA934294/bam"
THREADS=4

mkdir -p "$BAM_DIR"

# All 9 samples
SAMPLES=(
    SRR23438720 SRR23438721 SRR23438722  # Heat (32°C)
    SRR23438723 SRR23438724 SRR23438725  # Medium (24°C)
    SRR23438726 SRR23438727 SRR23438728  # Normal (13°C)
)

echo "=== STAR Alignment ==="
for srr in "${SAMPLES[@]}"; do
    echo "Aligning $srr..."
    STAR --genomeDir "$INDEX_DIR" \
        --readFilesIn "$TRIMMED_DIR/${srr}_clean_1.fastq" "$TRIMMED_DIR/${srr}_clean_2.fastq" \
        --outFileNamePrefix "$BAM_DIR/${srr}_" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard \
        --runThreadN $THREADS \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000

    echo "  $srr done."
done

echo "=== Alignment complete ==="
for srr in "${SAMPLES[@]}"; do
    if [ -f "$BAM_DIR/${srr}_Aligned.sortedByCoord.out.bam" ]; then
        SIZE=$(du -sh "$BAM_DIR/${srr}_Aligned.sortedByCoord.out.bam" | awk '{print $1}')
        echo "  $srr BAM: $SIZE"
    fi
done
