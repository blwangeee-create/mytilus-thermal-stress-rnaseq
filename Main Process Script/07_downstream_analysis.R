#!/usr/bin/env Rscript

# ============================================================
# Phase 7: Downstream Analysis
# GO/KEGG Enrichment, Clustering, PPI prep
# ============================================================

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

results_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results"
dir.create(file.path(results_dir, "enrichment"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(results_dir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load DESeq2 results ----
rlog_file <- file.path(results_dir, "rlog_normalized_counts.csv")
rlog_mat <- as.matrix(read.csv(rlog_file, row.names = 1))

# ---- 2. GO/KEGG Enrichment for each comparison ----
# NOTE: For Mytilus coruscus, we need to map gene IDs.
# Strategy: BLAST to SwissProt or use InterProScan for GO annotation.
# Since we may not have org.db for Mytilus, we'll:
# 1) Extract gene sequences from the annotation
# 2) BLAST against SwissProt for GO mapping (or use online tools)

cat("For Mytilus coruscus GO/KEGG enrichment:\n")
cat("Option 1: Use online tools (recommended for non-model species)\n")
cat("  - DAVID: https://david.ncifcrf.gov/\n")
cat("  - KOBAS: http://kobas.cbi.pku.edu.cn/\n")
cat("  - ShinyGO: http://bioinformatics.sdstate.edu/go/\n")
cat("  - eggNOG-mapper: http://eggnog-mapper.embl.de/\n")
cat("\n")
cat("Option 2: Use InterProScan locally for GO mapping\n")
cat("  interscan.sh -i proteins.faa -f TSV -o interpro.tsv\n")

# ---- 3. Heatmap of top DEGs ----
comparisons <- c("Heat_vs_Normal", "Medium_vs_Normal", "Heat_vs_Medium")

for (comp in comparisons) {
  deg_file <- file.path(results_dir, paste0("DEGs_", comp, "_significant.csv"))
  if (!file.exists(deg_file)) next

  degs <- read.csv(deg_file, row.names = 1)

  if (nrow(degs) < 5) {
    cat("Too few DEGs for", comp, "heatmap\n")
    next
  }

  cat(sprintf("%s: %d DEGs\n", comp, nrow(degs)))

  # Top 50 genes heatmap
  top50 <- head(degs[order(degs$padj), ], 50)
  top50_genes <- rownames(top50)

  # Extract rlog values for these genes
  rlog_subset <- rlog_mat[top50_genes, , drop = FALSE]
  rlog_subset <- rlog_subset[complete.cases(rlog_subset), ]

  if (nrow(rlog_subset) < 5) next

  # Z-score normalization
  rlog_z <- t(scale(t(rlog_subset)))

  pdf(file.path(results_dir, paste0("heatmap_top50_", comp, ".pdf")),
      width = 10, height = 14)
  pheatmap(rlog_z,
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           show_rownames = FALSE,
           main = paste("Top 50 DEGs -", comp),
           annotation_col = data.frame(
             Condition = factor(c(rep("Normal",3), rep("Medium",3), rep("Heat",3))),
             row.names = colnames(rlog_mat)
           ),
           annotation_colors = list(
             Condition = c(Normal = "#2E86AB", Medium = "#F18F01", Heat = "#C73E1D")
           ))
  dev.off()
}

# ---- 4. Identify apoptosis and circadian rhythm genes ----
# Known gene families from literature
apoptosis_keywords <- c("caspase", "CASP", "bcl", "BCL2", "bax", "BAX",
                        "bak", "BAK", "bid", "BID", "bad", "BAD",
                        "p53", "tp53", "TP53", "apaf", "APAF",
                        "cytochrome c", "CYCS", "fas", "FAS", "fasl",
                        "traf", "TRAF", "iap", "IAP", "survivin",
                        "diablo", "smac", "SMAC", "xiap", "XIAP",
                        "parp", "PARP", "aif", "AIF", "endog",
                        "bid", "BID", "noxa", "NOXA", "puma", "PUMA")

circadian_keywords <- c("clock", "CLOCK", "bmal", "BMAL", "bmal1", "BMAL1",
                        "period", "per", "PER", "timeless", "tim", "TIM",
                        "cryptochrome", "cry", "CRY", "npas", "NPAS",
                        "ror", "ROR", "rev-erba", "rev-erbb",
                        "nr1d1", "nr1d2", "rev-erba", "rev-erbb",
                        "dbp", "DBP", "hlf", "HLF", "tef", "TEF",
                        "csnk1", "csnk1d", "csnk1e", "ck1",
                        "fbxl", "fbxw", "arf", "dec",
                        "timeless", "cycle", "cyc", "Clk",
                        "doubletime", "dbt")

# Search for these keywords in gene names from the annotation
gff_file <- "/Volumes/Expansion/NCBI/PRJNA934294/liftoff/Mcoruscus_HiC_liftoff.gff3"
gff_lines <- readLines(gff_file, n = 1000)  # Preview first 1000 lines
cat("\nGFF3 preview:\n")
for (line in head(gff_lines[grep("gene_biotype=protein_coding", gff_lines)], 10)) {
  cat(substr(line, 1, 200), "\n")
}

cat("\n=== Downstream analysis script complete ===\n")
cat("Key results files:\n")
cat("  - DEG tables: ", file.path(results_dir, "DEGs_*.csv\n"))
cat("  - PCA plot: ", file.path(results_dir, "PCA_plot.pdf\n"))
cat("  - Volcano plots: ", file.path(results_dir, "volcano_*.pdf\n"))
cat("  - Heatmaps: ", file.path(results_dir, "heatmap_*.pdf\n"))
