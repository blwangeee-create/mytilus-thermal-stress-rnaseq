rm(list = ls())

# ================= 1. 安装并加载必要包 =================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")

packages <- c("ggplot2", "ggrepel", "ggforce", "patchwork", "aplot")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) install.packages(pkg)
}

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(patchwork)
library(aplot)

# ================= 2. 加载数据 =================
search_folder <- "/Volumes/Expansion/NCBI/PRJNA934294"
found <- list.files(search_folder, pattern = "counts_s0\\.txt$", 
                    recursive = TRUE, full.names = TRUE)
if (length(found) == 0) stop("counts_s0.txt not found.")
counts_file <- found[1]

counts_raw <- read.table(counts_file, header = TRUE, row.names = 1, 
                         check.names = FALSE, comment.char = "#")
counts <- counts_raw[, 6:ncol(counts_raw)]
colnames(counts) <- gsub("\\.sorted\\.bam$", "", basename(colnames(counts)))

# ================= 3. 样本信息 =================
sample_info <- data.frame(
  row.names = c("SRR23438720","SRR23438721","SRR23438722",
                "SRR23438723","SRR23438724","SRR23438725",
                "SRR23438726","SRR23438727","SRR23438728"),
  Condition = c("Treatment_32","Treatment_32","Treatment_32",
                "Treatment_24","Treatment_24","Treatment_24",
                "Control","Control","Control")
)
sample_info$Condition <- factor(sample_info$Condition,
                                levels = c("Control", "Treatment_24", "Treatment_32"))
counts <- counts[, rownames(sample_info)]

# ================= 4. DESeq2 + VST =================
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ Condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# ================= 5. 手动计算 PCA（替代 plotPCA，更灵活） =================
mat <- assay(vsd)
# 筛选高变异基因（前5000），与 DESeq2::plotPCA 逻辑一致
rv <- rowVars(mat)
select <- order(rv, decreasing = TRUE)[seq_len(min(5000, length(rv)))]
mat_pca <- t(mat[select, ])

pca_res <- prcomp(mat_pca, scale. = TRUE, center = TRUE)

# 提取坐标
pca_df <- as.data.frame(pca_res$x[, 1:2])
colnames(pca_df) <- c("PC1", "PC2")
pca_df$Condition <- sample_info$Condition
pca_df$Sample <- c("T32-1", "T32-2", "T32-3",
                   "T24-1", "T24-2", "T24-3",
                   "Ctrl-1", "Ctrl-2", "Ctrl-3")

# 方差贡献
percentVar <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

# ================= 6. 颜色方案 =================
pal <- c(
  Control        = "#6A757E",
  Treatment_24   = "#0B8BEE",
  Treatment_32   = "#FFCD00"
)

# ================= 6.5 计算凸包（base R，无冲突）=================
hull_list <- by(pca_df, pca_df$Condition, function(df) df[chull(df$PC1, df$PC2), ])
hull_df <- do.call(rbind, hull_list)
rownames(hull_df) <- NULL

# ================= 7. 基础 PCA 图 =================
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Condition, fill = Condition)) +
  
  geom_polygon(
    data = hull_df,
    aes(x = PC1, y = PC2, fill = Condition, colour = Condition),
    alpha = 0.15,
    linewidth = 0.8
  ) +
  
  geom_point(
    aes(fill = Condition),
    shape = 21,
    size = 5,
    stroke = 0.8,
    colour = "black"
  ) +
  
  geom_text_repel(
    aes(label = Sample),
    size = 3.5,
    colour = "black",
    max.overlaps = 20,
    box.padding = 0.45,
    point.padding = 0.35,
    seed = 123
  ) +
  
  scale_fill_manual(values = pal, name = "") +
  scale_colour_manual(values = pal, name = "") +
  
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  
  theme_bw(base_size = 13) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 11),
    legend.title = element_blank(),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(colour = "black", size = 11),
    axis.line = element_line(linewidth = 0.8, colour = "black"),
    axis.ticks = element_line(linewidth = 0.8, colour = "black"),
    panel.grid = element_blank()
  )

print(p_pca)

# ================= 8. 组合图（密度图 + 箱线图 + PCA，参考代码风格） =================

# 上方：PC1 密度图
p_density <- ggplot(pca_df, aes(x = PC1, fill = Condition)) +
  geom_density(alpha = 0.6) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  )

# 右侧：PC2 箱线图
p_box <- ggplot(pca_df, aes(x = Condition, y = PC2, fill = Condition)) +
  geom_boxplot(alpha = 0.75) +
  scale_fill_manual(values = pal) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# 组合（使用 aplot，参考代码风格）
p_combine <- p_pca %>%
  insert_top(p_density, height = 0.25) %>%
  insert_right(p_box, width = 0.25)

print(p_combine)