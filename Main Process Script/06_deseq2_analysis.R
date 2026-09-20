#!/usr/bin/env Rscript

# ============================================================
# Phase 6: DESeq2 Differential Expression Analysis
# Mytilus coruscus Heat Stress - 3-group comparison
# ============================================================

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(tibble)

# ---- 1. Load data ----
counts_file <- "/Volumes/Expansion/NCBI/PRJNA934294/counts/counts_matrix.txt"
counts <- read.table(counts_file, header = TRUE, row.names = 1, check.names = FALSE)

# Remove the header summary lines from featureCounts output
# First column should be gene IDs, rest are sample counts
colnames(counts) <- gsub(".*(SRR\\d+).*", "\\1", colnames(counts))

# Rename columns to meaningful sample names
sample_names <- c(
  "SRR23438728" = "Normal_1", "SRR23438727" = "Normal_2", "SRR23438726" = "Normal_3",
  "SRR23438725" = "Medium_1", "SRR23438724" = "Medium_2", "SRR23438723" = "Medium_3",
  "SRR23438722" = "Heat_1",    "SRR23438721" = "Heat_2",    "SRR23438720" = "Heat_3"
)

# Map existing column names
for (old in names(sample_names)) {
  colnames(counts) <- gsub(old, sample_names[old], colnames(counts))
}

# ---- 2. Create sample metadata ----
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c(
    rep("Normal", 3), rep("Medium", 3), rep("Heat", 3)
  ), levels = c("Normal", "Medium", "Heat")),
  replicate = factor(c(1,2,3, 1,2,3, 1,2,3))
)

cat("Sample metadata:\n")
print(coldata)
cat("\nCount matrix dimensions:", dim(counts), "\n")

# ---- 3. DESeq2 analysis ----
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition
)

# Pre-filtering: keep genes with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
cat("Genes after filtering:", nrow(dds), "\n")

# Run DESeq2
dds <- DESeq(dds)

# ---- 4. Results for three comparisons ----
results_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

comparisons <- list(
  "Heat_vs_Normal" = c("condition", "Heat", "Normal"),
  "Medium_vs_Normal" = c("condition", "Medium", "Normal"),
  "Heat_vs_Medium" = c("condition", "Heat", "Medium")
)

res_list <- list()

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  cat("\n=== ", comp_name, " ===\n")

  res <- results(dds, contrast = comp, alpha = 0.05)
  res <- res[order(res$padj), ]

  # Summary
  cat(sprintf("Total DEGs (padj < 0.05): %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
  cat(sprintf("  Up-regulated: %d\n", sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)))
  cat(sprintf("  Down-regulated: %d\n", sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)))

  # Save full results
  write.csv(as.data.frame(res), file.path(results_dir, paste0("DEGs_", comp_name, ".csv")))

  # Save significant DEGs
  sig_genes <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
  write.csv(as.data.frame(sig_genes), file.path(results_dir, paste0("DEGs_", comp_name, "_significant.csv")))

  res_list[[comp_name]] <- res
}

# ---- 5. Save normalized counts ----
rld <- rlog(dds, blind = FALSE)
rlog_mat <- assay(rld)
write.csv(rlog_mat, file.path(results_dir, "rlog_normalized_counts.csv"))

# ---- 6. PCA plot ----
pdf(file.path(results_dir, "PCA_plot.pdf"), width = 8, height = 6)
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA - Mytilus coruscus Heat Stress") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("Normal" = "#2E86AB", "Medium" = "#F18F01", "Heat" = "#C73E1D"))
print(p)
dev.off()

# ---- 7. Sample distance heatmap ----
pdf(file.path(results_dir, "sample_distance_heatmap.pdf"), width = 8, height = 7)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rld)
colnames(sampleDistMatrix) <- colnames(rld)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = "Sample distance matrix")
dev.off()

# ---- 8. Volcano plots for each comparison ----
for (comp_name in names(res_list)) {
  pdf(file.path(results_dir, paste0("volcano_", comp_name, ".pdf")), width = 8, height = 7)
  res_df <- as.data.frame(res_list[[comp_name]])

  # Add significance categories
  res_df$sig <- "NS"
  res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
  res_df$sig[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

  # Top 20 genes to label
  top_genes <- res_df[order(res_df$padj), ][1:20, ]

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "#2E86AB", "Up" = "#C73E1D", "NS" = "grey70")) +
    geom_text_repel(data = top_genes,
                    aes(label = rownames(top_genes)),
                    size = 3, max.overlaps = 15, color = "black") +
    theme_bw(base_size = 14) +
    labs(title = comp_name, x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5)
  print(p)
  dev.off()
}

cat("\n\n=== DESeq2 analysis complete ===\n")
cat("Results saved to:", results_dir, "\n")
