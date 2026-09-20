#!/bin/bash
set -euo pipefail

# ============================================================
# Phase 5: featureCounts Quantification
# ============================================================

BAM_DIR="/Volumes/Expansion/NCBI/PRJNA934294/bam"
ANNO_DIR="/Volumes/Expansion/NCBI/PRJNA934294/liftoff"
COUNTS_DIR="/Volumes/Expansion/NCBI/PRJNA934294/counts"
THREADS=8

mkdir -p "$COUNTS_DIR"

# Collect all BAM files
BAM_FILES=$(ls "$BAM_DIR/"*_Aligned.sortedByCoord.out.bam)

echo "=== featureCounts: Gene-level quantification ==="
featureCounts -T $THREADS \
    -p \
    -B \
    -C \
    -t exon \
    -g gene_id \
    -a "$ANNO_DIR/Mcoruscus_HiC_liftoff.gff3" \
    -o "$COUNTS_DIR/gene_counts.txt" \
    $BAM_FILES

echo "=== Quantification complete ==="
echo "Count matrix: $COUNTS_DIR/gene_counts.txt"
head -5 "$COUNTS_DIR/gene_counts.txt"
echo "..."
wc -l "$COUNTS_DIR/gene_counts.txt"

# Also create a simplified count matrix for R
echo "=== Creating simplified count matrix ==="
cut -f1,7- "$COUNTS_DIR/gene_counts.txt" > "$COUNTS_DIR/counts_matrix.txt"
echo "Simplified matrix: $COUNTS_DIR/counts_matrix.txt"

# Generate alignment statistics
echo "=== Alignment statistics ==="
for bam in "$BAM_DIR/"*_Aligned.sortedByCoord.out.bam; do
    SNAME=$(basename "$bam" | sed 's/_Aligned.sortedByCoord.out.bam//')
    TOTAL=$(samtools view -c "$bam" 2>/dev/null || echo "0")
    echo "  $SNAME: $TOTAL reads"
done
