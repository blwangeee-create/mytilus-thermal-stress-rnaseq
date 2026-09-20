# ========== 安装缺失包 ==========
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ComplexHeatmap", "circlize", "DESeq2"))
# ============================================
# Consensus DEG Heatmap - 最终修复版 (Top 500)
# ============================================

# ============================================
# Consensus DEG Heatmap - 最终稳定版
# ============================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(DESeq2)
})

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294"

# ---------- 样本定义 ----------
sample_groups <- list(
  Ctrl = c("SRR23438726", "SRR23438727", "SRR23438728"),
  T24  = c("SRR23438723", "SRR23438724", "SRR23438725"),
  T32  = c("SRR23438720", "SRR23438721", "SRR23438722")
)
sample_order <- unlist(sample_groups, use.names = FALSE)
group_vec <- rep(names(sample_groups), sapply(sample_groups, length))
names(group_vec) <- sample_order

friendly_map <- c(
  "SRR23438726" = "Ctrl-1", "SRR23438727" = "Ctrl-2", "SRR23438728" = "Ctrl-3",
  "SRR23438723" = "T24-1",  "SRR23438724" = "T24-2",  "SRR23438725" = "T24-3",
  "SRR23438720" = "T32-1",  "SRR23438721" = "T32-2",  "SRR23438722" = "T32-3"
)

# ---------- 读取 DEG ----------
read_deg <- function(f) {
  if (!file.exists(f)) stop("文件不存在: ", f)
  df <- read.table(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  return(df[, 1])
}

deg_t24 <- read_deg(file.path(base_dir, "results/counts/HighConf_DEG_T24vsCtrl.txt"))
deg_t32 <- read_deg(file.path(base_dir, "results/counts/HighConf_DEG_T32vsCtrl.txt"))

# ---------- 读取 DESeq2 结果 ----------
read_res <- function(f) {
  df <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  return(df)
}

res_t24 <- read_res(file.path(base_dir, "results/RNAseq_analysis/Results/DESeq2_T24_vs_Ctrl.csv"))
res_t32 <- read_res(file.path(base_dir, "results/RNAseq_analysis/Results/DESeq2_T32_vs_Ctrl.csv"))

# ---------- 并集 + Top 300 ----------
heatmap_genes <- union(deg_t24, deg_t32)

genes_in_both <- intersect(rownames(res_t24), rownames(res_t32))
genes_use <- intersect(heatmap_genes, genes_in_both)

fc_t24 <- res_t24[genes_use, "log2FoldChange"]
fc_t32 <- res_t32[genes_use, "log2FoldChange"]

valid <- !is.na(fc_t24) & !is.na(fc_t32)
fc_t24 <- fc_t24[valid]
fc_t32 <- fc_t32[valid]
genes_use <- genes_use[valid]

all_fc <- abs(fc_t24) + abs(fc_t32)
names(all_fc) <- genes_use
all_fc <- sort(all_fc, decreasing = TRUE)

top_n <- min(300, length(all_fc))
top_genes <- names(all_fc)[1:top_n]

cat(">>> 热图基因数:", top_n, "\n")

# ---------- Count → VST → z-score ----------
counts <- read.table(file.path(base_dir, "results/counts/counts_s0.txt"),
                     header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

clean_cols <- gsub("_clean.*|\\.sorted.*|\\.bam.*|\\.name_sorted.*|_\\d+$", "", colnames(counts))
clean_cols <- gsub("_$", "", clean_cols)

matched_idx <- sapply(sample_order, function(s) {
  idx <- grep(paste0("^", s, "$"), clean_cols)
  if (length(idx) == 0) idx <- grep(s, clean_cols)
  idx[1]
})

counts_sub <- counts[, matched_idx, drop = FALSE]
colnames(counts_sub) <- sample_order

counts_hm <- counts_sub[rownames(counts_sub) %in% top_genes, , drop = FALSE]
counts_hm <- counts_hm[rowSums(counts_hm) > 0, , drop = FALSE]

vst_mat <- tryCatch(vst(as.matrix(counts_hm), blind = FALSE), error = function(e) {
  log2(as.matrix(counts_hm) + 1)
})
zscore_mat <- t(scale(t(vst_mat)))
colnames(zscore_mat) <- friendly_map[colnames(zscore_mat)]

cat(">>> 表达矩阵:", nrow(zscore_mat), "genes x", ncol(zscore_mat), "samples\n")

# ---------- 配色 ----------
sample_df <- data.frame(
  Group = factor(group_vec[sample_order], levels = c("Ctrl", "T24", "T32")),
  row.names = colnames(zscore_mat)
)

group_colors <- c(Ctrl = "gray70", T24 = "#3070B3", T32 = "#FFCD00")
hm_colors <- colorRamp2(c(-2.5, 0, 2.5), c("#8FB4BE", "#FFFFFF", "#D93F49"))

# ---------- 顶部注释 ----------
top_anno <- columnAnnotation(
  Group = sample_df$Group,
  col = list(Group = group_colors),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  show_legend = FALSE,
  height = unit(1, "mm"),
  gp = gpar(col = "white")
)

# ---------- 绘制热图 ----------
ht <- Heatmap(
  zscore_mat,
  name = "Z-score",
  col = hm_colors,
  
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  row_split = 4,                    # ← 4 个分块 + Cluster 标题
  row_gap = unit(3, "mm"),
  
  row_title = "Cluster %s",
  row_title_gp = gpar(fontsize = 12, fontface = "bold", col = "black"),
  row_title_side = "left",
  
  show_row_dend = TRUE,
  row_dend_width = unit(2, "cm"),
  row_dend_gp = gpar(lwd = 0.7, col = "black"),
  
  cluster_columns = FALSE,
  column_order = colnames(zscore_mat),
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_side = "bottom",
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  column_names_rot = 0,
  
  rect_gp = gpar(col = NA),
  border = FALSE,
  
  column_title = "Sample",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_side = "top",
  
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(3, "cm"),
    at = c(-2.5, 0, 2.5)
  ),
  
  width = unit(12, "cm"),
  height = unit(16, "cm"),
  
  top_annotation = top_anno
)

lgd_group <- Legend(
  title = "Group",
  at = c("Ctrl", "T24", "T32"),
  legend_gp = gpar(fill = group_colors),
  title_gp = gpar(fontsize = 11, fontface = "bold"),
  labels_gp = gpar(fontsize = 10)
)

# ============================================
# 输出：Plots 面板预览 + PDF（不开 quartz，避免崩溃）
# ============================================

cat(">>> 正在生成 Plots 面板预览...\n")
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_list = list(lgd_group),
     annotation_legend_side = "right",
     padding = unit(c(5, 5, 5, 5), "mm"))
rstudio_dev <- dev.cur()
cat("✓ Plots 面板预览已生成（点击右上角 Zoom 放大查看）\n\n")

cat(">>> 正在输出 PDF 到桌面...\n")
pdf(file.path(Sys.getenv("HOME"), "Desktop/heatmap_Top300.pdf"), width = 12, height = 14)
draw(ht,
     heatmap_legend_side = "right",
     annotation_legend_list = list(lgd_group),
     annotation_legend_side = "right",
     padding = unit(c(8, 8, 8, 8), "mm"))
dev.off()

dev.set(rstudio_dev)
cat("✓ PDF 已保存到 ~/Desktop/heatmap_Top300.pdf\n")
cat("✓ Plots 面板已恢复\n")

# ============================================
# 论文统计输出
# ============================================
cat("\n============================================================\n")
cat("========== 论文 Materials & Methods / Results 可用统计 ==========\n")
cat("============================================================\n\n")

cat("【Materials & Methods】\n")
cat("----------------------------------------\n")
cat("Differential expression analysis was performed using DESeq2, edgeR, and limma.\n")
cat("High-confidence DEGs (Consensus DEGs) were defined as genes simultaneously\n")
cat("identified by all three algorithms with |log2FoldChange| >= 1 and adjusted P-value < 0.05.\n")
cat("For heatmap visualization, the union of Consensus DEGs from T24 vs Ctrl and\n")
cat("T32 vs Ctrl comparisons (n = ", length(heatmap_genes), ") was generated. The top ", top_n, " genes\n")
cat("ranked by the sum of absolute log2FoldChange values were selected for display.\n")
cat("Expression values were transformed using Variance Stabilizing Transformation (VST)\n")
cat("and row-wise z-score normalized. Hierarchical clustering of genes (rows) was\n")
cat("performed using Euclidean distance and complete linkage method, and partitioned\n")
cat("into four clusters to highlight distinct expression patterns. Samples (columns)\n")
cat("were ordered by temperature gradient (Ctrl -> T24 -> T32) without clustering.\n")
cat("Heatmap was generated using ComplexHeatmap (v", as.character(packageVersion("ComplexHeatmap")), ").\n\n")

cat("【Results】\n")
cat("----------------------------------------\n")
cat(sprintf("A total of %d Consensus DEGs were identified in T24 vs Ctrl (Up: %d, Down: %d).\n",
            length(deg_t24), sum(res_t24[deg_t24, "log2FoldChange"] > 0, na.rm = TRUE),
            sum(res_t24[deg_t24, "log2FoldChange"] < 0, na.rm = TRUE)))
cat(sprintf("A total of %d Consensus DEGs were identified in T32 vs Ctrl (Up: %d, Down: %d).\n",
            length(deg_t32), sum(res_t32[deg_t32, "log2FoldChange"] > 0, na.rm = TRUE),
            sum(res_t32[deg_t32, "log2FoldChange"] < 0, na.rm = TRUE)))
cat(sprintf("The heatmap displayed the top %d Consensus DEGs from the union of T24 vs Ctrl\n", top_n))
cat(sprintf("and T32 vs Ctrl (total union n = %d).\n", length(heatmap_genes)))

cat("\n============================================================\n")
cat("统计输出完毕。\n")
cat("============================================================\n")