#!/bin/bash
set -euo pipefail

# ============================================================
# Phase 3: STAR Genome Index Building
# ============================================================

GENOME_DIR="/Volumes/Expansion/NCBI/PRJNA934294/genome"
ANNO_DIR="/Volumes/Expansion/NCBI/PRJNA934294/liftoff"
INDEX_DIR="/Volumes/Expansion/NCBI/PRJNA934294/index/star"
THREADS=8
MEM_LIMIT=32000000000  # 32GB RAM limit

mkdir -p "$INDEX_DIR"

echo "=== Building STAR genome index ==="
echo "Genome: GCA_017311375.1 (HiC)"
echo "Annotation: Liftoff GFF3 (48,545 genes)"
echo "Threads: $THREADS"
echo "Index dir: $INDEX_DIR"

STAR --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_DIR/GCA_017311375.1_Mcoruscus_HiC_genomic.fna" \
    --sjdbGTFfile "$ANNO_DIR/Mcoruscus_HiC_liftoff.gff3" \
    --sjdbOverhang 149 \
    --runThreadN $THREADS \
    --limitGenomeGenerateRAM $MEM_LIMIT \
    --genomeSAindexNbases 12

echo "=== STAR index built ==="
ls -lh "$INDEX_DIR/"
echo "Index size: $(du -sh $INDEX_DIR | awk '{print $1}')"
