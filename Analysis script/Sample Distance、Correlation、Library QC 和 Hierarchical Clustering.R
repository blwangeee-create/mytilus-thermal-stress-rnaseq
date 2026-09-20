# ============================================================
# RNA-seq QC Save & Export: 保存所有结果和图表
# ============================================================

required_packages <- c("DESeq2", "pheatmap", "RColorBrewer", "ggplot2", 
                       "dplyr", "factoextra", "dendextend")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 创建输出目录
qc_outdir <- "./results/RNAseq_analysis/Figures/QC"
res_outdir <- "./results/RNAseq_analysis/Results"
dir.create(qc_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(res_outdir, showWarnings = FALSE, recursive = TRUE)

cat("Output directories:\n")
cat("  Figures:", qc_outdir, "\n")
cat("  Results:", res_outdir, "\n\n")

# 1. 读取 Counts -----------------------------------------------------------
counts_raw <- read.table("./results/counts/counts_s0.txt", header = TRUE, 
                         sep = "\t", check.names = FALSE, comment.char = "#")
counts_matrix <- as.matrix(counts_raw[, 7:ncol(counts_raw)])
rownames(counts_matrix) <- counts_raw$Geneid

# 2. 清理列名 --------------------------------------------------------------
clean_names <- gsub(".*/(SRR\\d+)\\.sorted\\.bam", "\\1", colnames(counts_matrix))
colnames(counts_matrix) <- clean_names

# 3. 样本元数据 ------------------------------------------------------------
sample_meta <- data.frame(
  SRR = c("SRR23438720", "SRR23438721", "SRR23438722",
          "SRR23438723", "SRR23438724", "SRR23438725",
          "SRR23438726", "SRR23438727", "SRR23438728"),
  SampleID = c("32-3", "32-2", "32-1", "24-3", "24-2", "24-1", "13-3", "13-2", "13-1"),
  Condition = factor(c(rep("T32", 3), rep("T24", 3), rep("Ctrl", 3)),
                     levels = c("Ctrl", "T24", "T32"))
)
sample_meta <- sample_meta[match(colnames(counts_matrix), sample_meta$SRR), ]
colnames(counts_matrix) <- sample_meta$SampleID

# 按条件排序
ord <- order(sample_meta$Condition)
counts_matrix <- counts_matrix[, ord]
sample_meta <- sample_meta[ord, ]
rownames(sample_meta) <- sample_meta$SampleID

cat("========== Sample Metadata ==========\n")
print(sample_meta)
cat("\n")

# 4. DESeq2 + 过滤 ---------------------------------------------------------
coldata <- data.frame(
  sample = sample_meta$SampleID,
  condition = sample_meta$Condition,
  row.names = sample_meta$SampleID
)
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep, ]
cat("After filtering:", nrow(dds), "genes retained.\n\n")

# 5. Library QC ------------------------------------------------------------
lib_qc <- data.frame(
  Sample = colnames(dds),
  Condition = coldata$condition,
  Total_Reads = colSums(counts(dds)),
  Detected_Genes = colSums(counts(dds) > 0),
  Detected_Genes_5reads = colSums(counts(dds) >= 5),
  Detected_Genes_10reads = colSums(counts(dds) >= 10),
  Million_Reads = round(colSums(counts(dds)) / 1e6, 2)
)

# 保存 Library QC
write.csv(lib_qc, file.path(res_outdir, "Library_QC_Summary.csv"), row.names = FALSE)
cat("Saved: Library_QC_Summary.csv\n")

# 6. VST + Ward.D2 聚类 ----------------------------------------------------
vsd <- vst(dds, blind = TRUE)
vst_mat <- assay(vsd)
sampleDists <- dist(t(vst_mat))
hc <- hclust(sampleDists, method = "ward.D2")

annotation_col <- data.frame(Condition = coldata$condition)
rownames(annotation_col) <- rownames(coldata)
ann_colors <- list(Condition = c(Ctrl = "#EBECEF", T24 = "#3070B3", T32 = "#FFCD00"))

# 7. Sample Distance --------------------------------------------------------
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

write.csv(sampleDistMatrix, file.path(res_outdir, "Sample_Distance_Matrix.csv"))
cat("Saved: Sample_Distance_Matrix.csv\n")

pdf(file.path(qc_outdir, "Fig_S1_Sample_Distance_Heatmap.pdf"), width = 8, height = 7)
pheatmap(
  sampleDistMatrix,
  cluster_rows = hc, cluster_cols = hc,
  treeheight_row = 50, treeheight_col = 30,
  cutree_rows = 3, cutree_cols = 3,
  border_color = "grey60",
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  annotation_col = annotation_col, annotation_row = annotation_col,
  annotation_colors = ann_colors,
  display_numbers = TRUE, number_color = "black", number_format = "%.2f",
  fontsize = 10,
  main = "Sample-to-Sample Distance (Euclidean on VST, Ward.D2)"
)
dev.off()
cat("Saved: Fig_S1_Sample_Distance_Heatmap.pdf\n")

# 8. Pearson Correlation ---------------------------------------------------
cor_matrix <- cor(vst_mat, method = "pearson")

write.csv(cor_matrix, file.path(res_outdir, "Sample_Correlation_Matrix.csv"))
cat("Saved: Sample_Correlation_Matrix.csv\n")

cor_colors <- colorRampPalette(c("white", "#C74375"))(100)
pdf(file.path(qc_outdir, "Fig_S2_Sample_Correlation_Heatmap.pdf"), width = 8, height = 7)
pheatmap(
  cor_matrix,
  cluster_rows = hc, cluster_cols = hc,
  treeheight_row = 50, treeheight_col = 30,
  cutree_rows = 3, cutree_cols = 3,
  border_color = "grey60",
  col = cor_colors,
  annotation_col = annotation_col, annotation_row = annotation_col,
  annotation_colors = ann_colors,
  display_numbers = TRUE, number_color = "black", number_format = "%.2f",
  fontsize = 10,
  main = "Sample Pearson Correlation (VST)"
)
dev.off()
cat("Saved: Fig_S2_Sample_Correlation_Heatmap.pdf\n")

# 9. Hierarchical Clustering -------------------------------------------------
cat("\n========== Clustering Order ==========\n")
cat("Order:", hc$order, "\n")
cat("Samples:", colnames(vsd)[hc$order], "\n\n")

# dendextend 旋转
dend <- as.dendrogram(hc)
ideal_order <- c("13-1", "13-2", "13-3", "24-1", "24-2", "24-3", "32-1", "32-2", "32-3")
present_order <- ideal_order[ideal_order %in% labels(dend)]
if (length(present_order) == length(labels(dend))) {
  dend <- rotate(dend, present_order)
  cat("Dendrogram rotated to group by condition.\n")
}

# 自动颜色映射
groups <- as.character(coldata$condition)
cutree_clusters <- cutree(hc, k = 3)
cluster_to_cond <- sapply(1:3, function(cl) {
  conds_in_cl <- groups[cutree_clusters == cl]
  names(sort(table(conds_in_cl), decreasing = TRUE))[1]
})
cond_colors <- c("Ctrl" = "#EBECEF", "T24" = "#3070B3", "T32" = "#FFCD00")
palette_ordered <- unname(cond_colors[cluster_to_cond])

cat("Cluster -> Condition:", paste(cluster_to_cond, collapse = ", "), "\n")
cat("Palette:", paste(palette_ordered, collapse = ", "), "\n\n")

pdf(file.path(qc_outdir, "Fig_1B_Hierarchical_Clustering_Dendrogram.pdf"), width = 8, height = 5)
p_dendro <- fviz_dend(
  dend, k = 3, cex = 0.9, palette = palette_ordered,
  labels_track_height = 0.3,
  xlab = "Sample", ylab = "Euclidean Distance (VST)",
  main = "Hierarchical Clustering of Samples (Ward.D2)"
) + theme(
  plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
  axis.text = element_text(size = 10)
)
print(p_dendro)
dev.off()
cat("Saved: Fig_1B_Hierarchical_Clustering_Dendrogram.pdf\n")

# 10. Library QC 图 ---------------------------------------------------------
pdf(file.path(qc_outdir, "Fig_S3_Library_Size.pdf"), width = 8, height = 5)
p1 <- ggplot(lib_qc, aes(x = Sample, y = Million_Reads, fill = Condition)) +
  geom_bar(stat = "identity", color = "grey40", width = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Ctrl" = "#EBECEF", "T24" = "#3070B3", "T32" = "#FFCD00")) +
  labs(title = "Library Size per Sample", y = "Total Reads (Million)", x = "") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)
dev.off()
cat("Saved: Fig_S3_Library_Size.pdf\n")

pdf(file.path(qc_outdir, "Fig_S4_Detected_Genes.pdf"), width = 8, height = 5)
p2 <- ggplot(lib_qc, aes(x = Sample, y = Detected_Genes_5reads, fill = Condition)) +
  geom_bar(stat = "identity", color = "grey40", width = 0.7) +
  scale_fill_manual(values = c("Ctrl" = "#EBECEF", "T24" = "#3070B3", "T32" = "#FFCD00")) +
  labs(title = "Detected Genes per Sample (>= 5 reads)", y = "Gene Count", x = "") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p2)
dev.off()
cat("Saved: Fig_S4_Detected_Genes.pdf\n")

# 11. 交叉验证 + 保存 RData ------------------------------------------------
cat("\n========== Cross-validation ==========\n")
get_group_stats <- function(mat, groups, metric_name) {
  within_vals <- c(); between_vals <- c(); n <- length(groups)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      val <- mat[i, j]
      if (groups[i] == groups[j]) { within_vals <- c(within_vals, val) }
      else { between_vals <- c(between_vals, val) }
    }
  }
  cat(sprintf("%s — Within-group: %.3f | Between-group: %.3f\n", 
              metric_name, mean(within_vals), mean(between_vals)))
}

get_group_stats(sampleDistMatrix, groups, "Distance")
get_group_stats(cor_matrix, groups, "Correlation")

# 保存 session
save(dds, vsd, sampleDistMatrix, cor_matrix, hc, lib_qc, dend, 
     file = file.path(res_outdir, "QC_analysis.RData"))
cat("\nSaved: QC_analysis.RData\n")

cat("\n========== All Files Saved ==========\n")
cat("Figures (PDF):", qc_outdir, "\n")
cat("Tables (CSV) & RData:", res_outdir, "\n")
cat("\nPlease paste the Cross-validation output above, and I will write M&M + Results.\n")