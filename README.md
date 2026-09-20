# Re-analysis of RNA-seq data from *Mytilus coruscus* under thermal stress

This repository contains the analysis code for a re-analysis of publicly available RNA-seq data from the hard-shelled mussel *Mytilus coruscus* exposed to three temperature conditions. The study integrates differential expression analysis, functional enrichment, co-expression network analysis, protein–protein interaction network analysis, qPCR validation, and protein structure prediction.

## Data source

The raw RNA-seq data were obtained from the NCBI Sequence Read Archive and were originally described in the source publication.

| Item | Value |
|---|---|
| BioProject | [PRJNA934294](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA934294) |
| SRA study | SRP422249 |
| Samples | 9 paired-end libraries (3 groups × 3 biological replicates) |
| Tissue | Gill |
| Treatment | 12 h at 13 °C (control), 24 °C, 32 °C |

| Group | Temperature | Samples |
|---|---|---|
| Control | 13 °C | SRR23438726, SRR23438727, SRR23438728 |
| T24 | 24 °C | SRR23438723, SRR23438724, SRR23438725 |
| T32 | 32 °C | SRR23438720, SRR23438721, SRR23438722 |

**Raw sequencing data are not redistributed here.** They are available from the SRA under the accessions listed above.

## Reference genome

| Item | Value |
|---|---|
| Assembly | GCA_011752425.2 (MCOR1.1) |
| Annotation | MCOR1.1 GFF/GTF |
| Functional annotation | eggNOG-mapper v2.1.12 |

## Pipeline overview

```
SRA (PRJNA934294)
  │
  ├─ fasterq-dump (SRA Toolkit v3.4.1)
  ▼
Raw FASTQ
  │
  ├─ fastp v1.1.0                              [scripts/02_qc_pipeline.sh]
  ▼
Clean reads
  │
  ├─ HISAT2 v2.2.2 → SAMtools v1.22.1          [scripts/04_hisat2_align.sh]
  ▼
Sorted BAM
  │
  ├─ featureCounts v2.1.1                      [scripts/05_featureCounts.sh]
  ▼
Count matrix
  │
  ├─ DESeq2 v1.52.0 + edgeR v3.42.4 + limma-voom v3.54.2
  ├─ GO / KEGG ORA + GSEA (clusterProfiler v4.20.0)
  ├─ WGCNA v1.72 (power = 18; ANOVA + BH for module–trait association)
  ├─ k-means clustering (K = 3)
  ├─ PPI network (STRING v12.0 + Cytoscape v3.10)
  ├─ Candidate gene identification (multi-evidence convergence)
  ├─ qPCR validation (TUBA + RPL13 dual reference, 2^-ΔΔCt)
  └─ Protein structure prediction (AlphaFold, SAVES v6.1, CABS-flex 2.0, ConSurf)
```

## Repository structure

```
├── README.md
├── LICENSE
├── sessionInfo.txt
├── scripts/                 # Upstream pipeline (shell)
│   ├── 01_download_convert.sh
│   ├── 02_qc_pipeline.sh
│   ├── run_fastp.sh
│   ├── 03_build_star_index.sh
│   ├── 04_hisat2_align.sh
│   ├── 05_featureCounts.sh
│   ├── 06_build_qc_summary.sh
│   ├── 07_build_final_qc_tables.sh
│   └── convert_gff3_for_featureCounts.py
├── analysis/                # Downstream analysis (R)
│   ├── 01_DESeq2/           # Differential expression, high-confidence DEG filter
│   ├── 02_GO_KEGG/          # GO/KEGG ORA, ORA×GSEA dual-positive test
│   ├── 03_GSEA/             # Gene set enrichment analysis
│   ├── 04_Figure/           # Volcano, MA, PCA, UpSet figures
│   ├── 02_heatmap/          # Consensus DEG heatmap
│   ├── K/                   # k-means clustering
│   ├── PPI/                 # STRING network preparation
│   ├── wgcna/               # WGCNA (original, corrected, plotting)
│   ├── primer&exon.R        # cDNA assembly and primer design
│   ├── protein structure.R  # AlphaFold-based structure analysis
│   └── Sample Distance.../  # Sample QC (distance, correlation, clustering)
└── qpcr/                    # qPCR analysis
    └── qpcr.R               # Dual-reference 2^-ΔΔCt, Kruskal-Wallis, Cliff's delta, PCA
```

## Software versions

### Upstream tools

| Tool | Version |
|---|---|
| SRA Toolkit | v3.4.1 |
| fastp | v1.1.0 |
| HISAT2 | v2.2.2 |
| SAMtools | v1.22.1 |
| featureCounts (Subread) | v2.1.1 |
| STAR | v2.7.11b (index built, alternate aligner) |

### R and packages

R version 4.6.0

| Package | Version |
|---|---|
| DESeq2 | v1.52.0 |
| edgeR | v3.42.4 |
| limma | v3.54.2 |
| clusterProfiler | v4.20.0 |
| WGCNA | v1.72 |
| ComplexUpset | v1.3.3 |
| ComplexHeatmap | v2.28.0 |
| pheatmap | v1.0.12 |
| factoextra | v1.0.7 |
| ggplot2 | v3.5.1 |
| org.Hs.eg.db | v3.23.1 |
| GO.db | v3.23.1 |

### Other software

| Software | Version |
|---|---|
| eggNOG-mapper | v2.1.12 |
| STRING | v12.0 |
| Cytoscape | v3.10 |
| PyMOL | v2.5 (Schrödinger LLC) |
| SAVES | v6.1 |
| CABS-flex | v2.0 |
| ConSurf | 2022 |

## Analysis parameters

### Differential expression

- Three algorithms: DESeq2 (Wald test), edgeR (quasi-likelihood F test), limma-voom (empirical Bayes)
- Filtering: genes with counts < 10 in more than 6 samples excluded
- Threshold: |log₂FC| > 1 and adjusted *p* < 0.05 (DESeq2: padj; edgeR: FDR; limma: adj.P.Val)
- High-confidence DEGs: intersection of at least two algorithms

### Functional enrichment

- ORA: clusterProfiler `enricher()`, foreground = consensus DEGs, background = all detected genes; BH correction, *p* < 0.05, *q* < 0.2
- GSEA: genes ranked by DESeq2 Wald statistic; minGSSize = 10, maxGSSize = 500, BH correction, padj < 0.05, `set.seed(123)`
- Dual-positive pathways: significant in both ORA (FDR < 0.05) and GSEA (FDR < 0.05, NES > 0)

### WGCNA

- Input: top 5,000 most variable genes (MAD) from the VST matrix
- Soft-thresholding power: **18** (recommended for signed networks with < 20 samples; signed R² = 0.821)
- Network type: signed; minModuleSize = 60; mergeCutHeight = 0.35
- **Module–trait association: one-way ANOVA (3 groups × 3 replicates, df = 2,6) with Benjamini-Hochberg correction**, using the three independent treatment groups as the unit of analysis rather than treating the nine samples as independent
- Hub genes: top 30 by kME
- Parameter sensitivity: power = 7 re-run for comparison

### qPCR

- Reference genes: TUBA + RPL13 (dual reference, geometric mean)
- Quantification: 2^-ΔΔCt, 13 °C as calibrator
- Replicates: 12 biological replicates per gene per temperature
- Statistics: Kruskal-Wallis test + pairwise Wilcoxon (Holm correction); effect size by Cliff's delta
- Consistency with RNA-seq: Spearman correlation

## Requirements

### System

- macOS (tested on macOS, Apple Silicon) or Linux
- 16 GB RAM minimum (WGCNA and structure analysis)
- ~50 GB disk for intermediate files (excluding raw FASTQ and BAM)

### Command-line tools

Install via conda (recommended):

```bash
conda create -n rnaseq_env -c bioconda -c conda-forge \
    sra-tools fastp fastqc multiqc hisat2 samtools subread star
```

### R packages

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
    "DESeq2", "edgeR", "limma", "clusterProfiler",
    "WGCNA", "ComplexHeatmap", "org.Hs.eg.db", "GO.db"
))

install.packages(c(
    "ggplot2", "dplyr", "tidyr", "pheatmap",
    "factoextra", "ComplexUpset", "ggrepel", "patchwork"
))
```

## Usage

The pipeline is written for a specific directory layout. Paths are hard-coded at the top of each script and must be adjusted to your environment before running.

Typical order:

```bash
# 1. Download and convert SRA to FASTQ
bash scripts/01_download_convert.sh

# 2. Quality control
bash scripts/run_fastp.sh
bash scripts/02_qc_pipeline.sh

# 3. Alignment
bash scripts/04_hisat2_align.sh

# 4. Quantification
bash scripts/05_featureCounts.sh

# 5. QC summary tables
bash scripts/06_build_qc_summary.sh
bash scripts/07_build_final_qc_tables.sh
```

```r
# 6. Differential expression
source("analysis/01_DESeq2/DESeq2_analysis.R")
source("analysis/01_DESeq2/DEG_high_confidence_filter.R")

# 7. Functional enrichment
source("analysis/02_GO_KEGG/GO_ORA_analysis.R")
source("analysis/02_GO_KEGG/KEGG_ORA_analysis.R")
source("analysis/02_GO_KEGG/dual positive test & go.R")
source("analysis/03_GSEA/GSEA_plot.R")

# 8. Co-expression and clustering
source("analysis/wgcna/wgcna_corrected.R")     # corrected version (ANOVA + power 7/18)
source("analysis/K/k-means.R")

# 9. PPI network
source("analysis/PPI/PPI.R")

# 10. Primer design
source("analysis/primer&exon.R")

# 11. Protein structure analysis
source("analysis/protein structure.R")

# 12. qPCR
source("qpcr/qpcr.R")
```

## Key analysis notes

### WGCNA module–trait association

The module–trait association was tested using one-way ANOVA with the three treatment groups (each with three biological replicates) as the unit of analysis, rather than applying a correlation-based test that assumes nine independent samples. This avoids pseudoreplication, which would overstate significance when the trait has only three distinct values.

At power = 18, a core heat-stress-responsive module (turquoise, n = 2,024) was identified with ANOVA FDR = 3.5 × 10⁻⁴, containing the chaperone genes DNAJA4, DNAJB5, HSPBP1, HSPA5, and AHSA1 among its hub genes.

### Parameter sensitivity

WGCNA was also run at power = 7 for comparison. Module partitions differed (20 vs 7 modules), but the core candidate genes were retained under both parameter settings.

## Data availability

- **Raw RNA-seq data**: NCBI SRA, BioProject [PRJNA934294](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA934294) (originally deposited by the source study)
- **Analysis results**: Zenodo, DOI [10.5281/zenodo.21890501]
- **Analysis code**: this repository

## Citation

[Manuscript in preparation]

## License

Code in this repository is released under the MIT License. See [LICENSE](LICENSE) for details.

## Contact

[Corresponding author contact to be added]
