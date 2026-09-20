# ========== 终极修复：确保 label_display 绝对不为 NA ==========

# 1. 重新安全构建所有注释列
res_df$mcor_id <- gsub("^gene-", "", res_df$gene_name)

# eggNOG 注释
res_df$gene_symbol <- combined_map$Preferred_name[match(res_df$mcor_id, combined_map$locus_tag)]
res_df$description <- combined_map$Description[match(res_df$mcor_id, combined_map$locus_tag)]

# 2. 构建 label_display：三层 fallback，最终保底为 mcor_id
res_df$label_display <- res_df$mcor_id  # 先全部设为基因ID（保底）

# 有 symbol 的覆盖
has_sym <- !is.na(res_df$gene_symbol) & res_df$gene_symbol != "" & res_df$gene_symbol != "-"
res_df$label_display[has_sym] <- res_df$gene_symbol[has_sym]

# 没有 symbol 但有 description 的覆盖
has_desc <- !has_sym & !is.na(res_df$description) & res_df$description != "" & res_df$description != "-"
res_df$label_display[has_desc] <- paste0(substr(res_df$description[has_desc], 1, 20), "...")

# 确认没有 NA
cat("label_display NA 数量：", sum(is.na(res_df$label_display)), "\n")

# 3. 重新提取 top 24
top_diff_genes <- res_df %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 12) %>%
  pull(gene_name)

# 4. 安全设置 to_label（不用 ifelse）
res_df$to_label <- ""
is_top <- res_df$gene_name %in% top_diff_genes
res_df$to_label[is_top] <- res_df$label_display[is_top]

# 再次确认
cat("to_label NA 数量：", sum(is.na(res_df$to_label)), "\n")
cat("to_label 非空数量：", sum(res_df$to_label != ""), "\n")

# 5. 查看结果
cat("\nTop 基因标签：\n")
print(res_df[is_top, c("mcor_id", "gene_symbol", "label_display", "to_label")])
# ========== 绘图（24个标签全显示） ==========

plot_labels <- res_df[res_df$to_label != "" & !is.na(res_df$to_label), ]
cat("即将绘制的标签数量：", nrow(plot_labels), "\n")

p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  
  geom_point(data = subset(res_df, group == "NS"), color = "grey80", alpha = 0.4, size = 1.7) +
  geom_point(data = subset(res_df, group == "Sig"), color = "#3070B3", alpha = 0.8, size = 1.7) +
  #geom_point(data = plot_labels, color = "#3070B3", size = 2.5, alpha = 0.9) +
  
  scale_x_log10(breaks = 10^seq(1, 5, by = 1),
                labels = c("1e+01", "1e+02", "1e+03", "1e+04", "1e+05")) +
  
  scale_y_continuous(breaks = seq(-20, 20, by = 10)) +
  coord_cartesian(ylim = c(-20, 20)) +
  
  geom_hline(yintercept = c(-1, 0, 1), 
             linetype = c("dashed", "solid", "dashed"),
             color = c("grey50", "black", "grey50"),
             linewidth = c(0.5, 0.8, 0.5)) +
  
  geom_text_repel(
    data = plot_labels,
    aes(label = to_label,
        color = ifelse(!is.na(gene_symbol) & gene_symbol != "", "black", "grey50"),
        fontface = ifelse(!is.na(gene_symbol) & gene_symbol != "", "bold", "plain")),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.15,
    segment.size = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    force = 0.3,
    force_pull = 1.2,
    direction = "both",
    seed = 42
  ) +
  
  scale_color_identity() +
  
  labs(
    x = "mean norm. expr.",
    y = expression(log[2]~fold~change),
    title = "<span style='color:#3070B3;'>T24</span> vs <span style='color:grey70;'>Ctrl</span>"
  ) +
  
  theme_few() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_markdown(size = 25, face = "bold", hjust = 0.5),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )

print(p)

# ========== 强制刷新绘图 ==========

# 1. 关闭所有图形设备
while (dev.cur() > 1) dev.off()

# 2. 清除 ggplot 缓存
rm(p)

# 3. 重新运行完整 T32 代码（上面给你的）
# ...（粘贴完整代码）...

# 4. 强制新开窗口显示
dev.new()
print(p)

# 5. 同时保存到新文件确认
png("MA_T32_test.png", width = 7, height = 5, units = "in", res = 300)
print(p)
dev.off()
cat("文件已保存到：", getwd(), "/MA_T32_test.png\n")
# ========== 验证：提取两张图的实际数据 ==========

# 先画 T24 并保存数据
rm(list = ls())  # 清空所有
load("Phase1_DEG_analysis.RData")

res_t24 <- as.data.frame(res_T24vsCtrl)
res_t24 <- res_t24[res_t24$baseMean > 10, ]
res_t24$sig <- res_t24$padj < 0.05 & abs(res_t24$log2FoldChange) > 1

# 再画 T32
res_t32 <- as.data.frame(res_T32vsCtrl)
res_t32 <- res_t32[res_t32$baseMean > 10, ]
res_t32$sig <- res_t32$padj < 0.05 & abs(res_t32$log2FoldChange) > 1
# ========== 保存两张图到不同文件，强制对比 ==========
# ========== 重新创建 combined_map ==========

# 1. 读取 eggNOG（假设路径没变）
eggnog <- read.delim(
  "/Volumes/Expansion/NCBI/PRJNA934294/annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular",
  header = FALSE, sep = "\t", comment.char = "#", check.names = FALSE
)

eggnog_cols <- c("query", "seed_ortholog", "evalue", "score", 
                 "eggNOG_OGs", "max_annot_lvl", "COG_category", 
                 "Description", "Preferred_name", "GOs", 
                 "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", 
                 "KEGG_Reaction", "KEGG_rclass", "BRITE", 
                 "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs")
colnames(eggnog) <- eggnog_cols

# 2. 读取 GFF（假设路径没变）
library(rtracklayer)
gff <- import("/Volumes/Expansion/NCBI/PRJNA934294/annotation/GCA_011752425.2_MCOR1.1_genomic.gff")
gff_df <- as.data.frame(gff)

# 3. 创建 combined_map
protein_rows <- gff_df[!is.na(gff_df$protein_id), ]
id_mapping <- unique(protein_rows[, c("protein_id", "locus_tag")])
id_mapping <- id_mapping[!is.na(id_mapping$locus_tag), ]

eggnog_map <- eggnog[, c("query", "Preferred_name")]
eggnog_map <- eggnog_map[eggnog_map$Preferred_name != "-" & !is.na(eggnog_map$Preferred_name), ]
names(eggnog_map)[1] <- "protein_id"

combined_map <- merge(id_mapping, eggnog_map, by = "protein_id", all.x = TRUE)
combined_map <- combined_map[!is.na(combined_map$Preferred_name), ]

cat("combined_map 重建完成，行数：", nrow(combined_map), "\n")
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggthemes)
library(ggtext)

# 函数：画 MA 图并保存
draw_and_save <- function(res_obj, color_hex, title_html, filename) {
  
  res_df <- as.data.frame(res_obj)
  res_df$gene_name <- rownames(res_df)
  res_df <- res_df[res_df$baseMean > 10, ]
  
  res_df$group <- "NS"
  res_df$group[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Sig"
  res_df$group <- factor(res_df$group, levels = c("NS", "Sig"))
  
  # 注释
  res_df$mcor_id <- gsub("^gene-", "", res_df$gene_name)
  res_df$gene_symbol <- combined_map$Preferred_name[match(res_df$mcor_id, combined_map$locus_tag)]
  res_df$description <- combined_map$Description[match(res_df$mcor_id, combined_map$locus_tag)]
  
  res_df$label_display <- res_df$mcor_id
  has_sym <- !is.na(res_df$gene_symbol) & res_df$gene_symbol != "" & res_df$gene_symbol != "-"
  res_df$label_display[has_sym] <- res_df$gene_symbol[has_sym]
  
  top_diff_genes <- res_df %>%
    filter(!is.na(padj) & padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    group_by(sign(log2FoldChange)) %>%
    slice_head(n = 12) %>%
    pull(gene_name)
  
  res_df$to_label <- ""
  is_top <- res_df$gene_name %in% top_diff_genes
  res_df$to_label[is_top] <- res_df$label_display[is_top]
  
  plot_labels <- res_df[res_df$to_label != "" & !is.na(res_df$to_label), ]
  
  p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(data = subset(res_df, group == "NS"), color = "grey80", alpha = 0.5, size = 1.7) +
    geom_point(data = subset(res_df, group == "Sig"), color = color_hex, alpha = 0.9, size = 1.7) +
    scale_x_log10(breaks = 10^seq(1, 5, by = 1),
                  labels = c("1e+01", "1e+02", "1e+03", "1e+04", "1e+05")) +
    scale_y_continuous(breaks = seq(-20, 20, by = 10)) +
    coord_cartesian(ylim = c(-20, 20)) +
    geom_hline(yintercept = c(-1, 0, 1), 
               linetype = c("dashed", "solid", "dashed"),
               color = c("grey50", "black", "grey50"),
               linewidth = c(0.5, 0.8, 0.5)) +
    geom_text_repel(
      data = plot_labels,
      aes(label = to_label,
          color = ifelse(!is.na(gene_symbol) & gene_symbol != "", "black", "grey50"),
          fontface = ifelse(!is.na(gene_symbol) & gene_symbol != "", "bold", "plain")),
      size = 4, max.overlaps = Inf, box.padding = 0.2, point.padding = 0.15,
      segment.size = 0.3, segment.color = "grey50", min.segment.length = 0,
      force = 0.3, force_pull = 1.2, direction = "both", seed = 42
    ) +
    scale_color_identity() +
    labs(x = "mean norm. expr.", y = expression(log[2]~fold~change), title = title_html) +
    theme_few() +
    theme(legend.position = "none",
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18, color = "black"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          plot.title = element_markdown(size = 25, face = "bold", hjust = 0.5),
          axis.ticks = element_line(linewidth = 0.8, color = "black"),
          axis.ticks.length = unit(0.3, "cm"))
  
  # 保存
  png(filename, width = 7, height = 5, units = "in", res = 300)
  print(p)
  dev.off()
  
  cat("✅ 已保存：", filename, "\n")
  cat("   显著基因数：", sum(res_df$group == "Sig"), "\n")
  cat("   标签数：", nrow(plot_labels), "\n")
  cat("   最大 |log2FC|：", round(max(abs(res_df$log2FoldChange), na.rm = TRUE), 2), "\n\n")
  
  return(p)
}

# ========== 画 T24 ==========
p24 <- draw_and_save(
  res_obj = res_T24vsCtrl,
  color_hex = "#3070B3",
  title_html = "<span style='color:#3070B3;'>T24</span> vs <span style='color:grey70;'>Ctrl</span>",
  filename = "MA_T24_FINAL.png"
)

# ========== 画 T32 ==========
p32 <- draw_and_save(
  res_obj = res_T32vsCtrl,
  color_hex = "#FFCD00",
  title_html = "<span style='color:#FFCD00;'>T32</span> vs <span style='color:grey70;'>Ctrl</span>",
  filename = "MA_T32_FINAL.png"
)

cat("两张图已保存到工作目录：", getwd(), "\n")
cat("请用图片查看器打开 MA_T24_FINAL.png 和 MA_T32_FINAL.png 对比\n")

# 对比
cat("T24 显著点数：", sum(res_t24$sig), "\n")
cat("T32 显著点数：", sum(res_t32$sig), "\n")

cat("\nT24 最大 log2FC：", max(abs(res_t24$log2FoldChange), na.rm = TRUE), "\n")
cat("T32 最大 log2FC：", max(abs(res_t32$log2FoldChange), na.rm = TRUE), "\n")

# 随机抽一个基因对比
set.seed(42)
random_gene <- sample(rownames(res_t24), 1)
cat("\n随机基因", random_gene, "：\n")
cat("T24 log2FC =", round(res_t24[random_gene, "log2FoldChange"], 3), "\n")
cat("T32 log2FC =", round(res_t32[random_gene, "log2FoldChange"], 3), "\n")
# ========== T24 预览 ==========
p_t24 <- draw_and_save(
  res_obj = res_T24vsCtrl,
  color_hex = "#3070B3",
  title_html = "<span style='color:#3070B3;'>T24</span> vs <span style='color:grey70;'>Ctrl</span>",
  filename = "MA_T24_FINAL.png"
)

# 只预览，不保存
print(p_t24)
# ========== T32 预览 ==========
p_t32 <- draw_and_save(
  res_obj = res_T32vsCtrl,
  color_hex = "#FFCD00",
  title_html = "<span style='color:#FFCD00;'>T32</span> vs <span style='color:grey70;'>Ctrl</span>",
  filename = "MA_T32_FINAL.png"
)

# 只预览，不保存
print(p_t32)
# ========== 获取论文所需数据 ==========

# T24
res_t24 <- as.data.frame(res_T24vsCtrl)
res_t24 <- res_t24[!is.na(res_t24$padj), ]
sig_t24 <- res_t24[res_t24$padj < 0.05 & abs(res_t24$log2FoldChange) > 1, ]

cat("=== T24 ===\n")
cat("总基因数：", nrow(res_t24), "\n")
cat("显著基因数：", nrow(sig_t24), "\n")
cat("上调：", sum(sig_t24$log2FoldChange > 0), "\n")
cat("下调：", sum(sig_t24$log2FoldChange < 0), "\n")
cat("Top 5 上调：\n")
print(head(sig_t24[sig_t24$log2FoldChange > 0, c("log2FoldChange", "padj")], 5))
cat("Top 5 下调：\n")
print(head(sig_t24[sig_t24$log2FoldChange < 0, c("log2FoldChange", "padj")], 5))

# T32
res_t32 <- as.data.frame(res_T32vsCtrl)
res_t32 <- res_t32[!is.na(res_t32$padj), ]
sig_t32 <- res_t32[res_t32$padj < 0.05 & abs(res_t32$log2FoldChange) > 1, ]

cat("\n=== T32 ===\n")
cat("总基因数：", nrow(res_t32), "\n")
cat("显著基因数：", nrow(sig_t32), "\n")
cat("上调：", sum(sig_t32$log2FoldChange > 0), "\n")
cat("下调：", sum(sig_t32$log2FoldChange < 0), "\n")
cat("Top 5 上调：\n")
print(head(sig_t32[sig_t32$log2FoldChange > 0, c("log2FoldChange", "padj")], 5))
cat("Top 5 下调：\n")
print(head(sig_t32[sig_t32$log2FoldChange < 0, c("log2FoldChange", "padj")], 5))







# ========== MA 图 (T32 vs T24) - 完整代码 ==========
while (dev.cur() > 1) dev.off()

library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggthemes)
library(ggtext)

# 数据准备
res_df <- as.data.frame(res_T32vsT24)
res_df$gene_name <- rownames(res_df)
res_df <- res_df[res_df$baseMean > 10, ]

# 分组
res_df$group <- "NS"
res_df$group[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1] <- "Sig"
res_df$group <- factor(res_df$group, levels = c("NS", "Sig"))

# 注释转换
res_df$mcor_id <- gsub("^gene-", "", res_df$gene_name)
res_df$gene_symbol <- combined_map$Preferred_name[match(res_df$mcor_id, combined_map$locus_tag)]
res_df$description <- combined_map$Description[match(res_df$mcor_id, combined_map$locus_tag)]

res_df$label_display <- res_df$mcor_id
has_sym <- !is.na(res_df$gene_symbol) & res_df$gene_symbol != "" & res_df$gene_symbol != "-"
res_df$label_display[has_sym] <- res_df$gene_symbol[has_sym]

# 提取 top 24
top_diff_genes <- res_df %>%
  filter(!is.na(padj) & padj < 0.05) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 12) %>%
  pull(gene_name)

res_df$to_label <- ""
is_top <- res_df$gene_name %in% top_diff_genes
res_df$to_label[is_top] <- res_df$label_display[is_top]

plot_labels <- res_df[res_df$to_label != "" & !is.na(res_df$to_label), ]

# 绘图
p <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  
  geom_point(data = subset(res_df, group == "NS"), 
             color = "grey80", alpha = 0.5, size = 1.7) +
  
  geom_point(data = subset(res_df, group == "Sig"), 
             color = "#D4247E", alpha = 0.8, size = 1.7) +
  
  scale_x_log10(
    breaks = 10^seq(1, 5, by = 1),
    labels = c("1e+01", "1e+02", "1e+03", "1e+04", "1e+05")
  ) +
  
  scale_y_continuous(breaks = seq(-20, 20, by = 10)) +
  coord_cartesian(ylim = c(-20, 20)) +
  
  geom_hline(yintercept = c(-1, 0, 1), 
             linetype = c("dashed", "solid", "dashed"),
             color = c("grey50", "black", "grey50"),
             linewidth = c(0.5, 0.8, 0.5)) +
  
  geom_text_repel(
    data = plot_labels,
    aes(label = to_label,
        color = ifelse(!is.na(gene_symbol) & gene_symbol != "", "black", "grey50"),
        fontface = ifelse(!is.na(gene_symbol) & gene_symbol != "", "bold", "plain")),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.2,
    point.padding = 0.15,
    segment.size = 0.3,
    segment.color = "grey50",
    min.segment.length = 0,
    force = 0.3,
    force_pull = 1.2,
    direction = "both",
    seed = 42
  ) +
  
  scale_color_identity() +
  
  labs(
    x = "mean norm. expr.",
    y = expression(log[2]~fold~change),
    title = "<span style='color:#FFCD00;'>T32</span> vs <span style='color:#3070B3;'>T24</span>"
  ) +
  
  theme_few() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_markdown(size = 25, face = "bold", hjust = 0.5),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.3, "cm")
  )

print(p)

cat("\nT32 vs T24 图已绘制\n")
cat("显著基因数：", sum(res_df$group == "Sig"), "\n")
cat("标签数：", nrow(plot_labels), "\n")