library(dplyr)
library(ComplexUpset)
library(ggplot2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

cat("========================================\n")
cat("UpSet 修复：严格按你的绘图代码预览 + Overlap修复\n")
cat("========================================\n\n")

# 加载变量
load("results/RNAseq_analysis/Results/Phase1_DEG_analysis.RData")

# ========================
# 你的绘图代码（原封不动）
# ========================

build_combined_matrix <- function(des_up, ed_up, lim_up, des_down, ed_down, lim_down) {
  all_genes <- unique(c(des_up, ed_up, lim_up, des_down, ed_down, lim_down))
  
  data.frame(
    gene      = all_genes,
    DESeq2    = as.integer(all_genes %in% c(des_up, des_down)),
    edgeR     = as.integer(all_genes %in% c(ed_up, ed_down)),
    limma     = as.integer(all_genes %in% c(lim_up, lim_down)),
    Direction = sapply(all_genes, function(g) {
      if (g %in% des_up || g %in% ed_up || g %in% lim_up) return("Up")
      if (g %in% des_down || g %in% ed_down || g %in% lim_down) return("Down")
      return(NA)
    })
  )
}

plot_upset_combined <- function(mat, title) {
  my_sets <- c('DESeq2', 'edgeR', 'limma')
  
  upset(
    mat,
    intersect = my_sets,
    sort_sets = FALSE,
    name      = title,
    
    stripes = upset_stripes(
      geom    = geom_segment(size = 7),
      colors  = c('grey90', 'white')
    ),
    
    matrix = intersection_matrix(
      geom    = geom_point(aes(color = in_set), size = 3.5),
      segment = geom_segment(size = 1.2)
    ) +
      scale_color_manual(
        values = c("TRUE" = "black", "FALSE" = "grey60"),
        guide  = "none"
      ),
    
    base_annotations = list(
      'Intersection size' = intersection_size(
        counts  = TRUE,
        mapping = aes(fill = Direction),
        bar_number_threshold = 1,
        text    = list(vjust = -0.5, size = 3.5)
      ) +
        scale_fill_manual(
          values = c('Up' = '#C74375', 'Down' = '#3070B3'),
          guide  = "none"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        theme(
          axis.ticks.y  = element_line(color = "black"),
          axis.line.y   = element_line(color = "black"),
          axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          axis.title.x  = element_blank(),
          panel.background = element_blank(),
          panel.grid   = element_blank()
        )
    ),
    
    set_sizes = upset_set_size(
      geom = geom_bar(
        aes(fill = Direction),
        width = 0.7
      )
    ) +
      scale_fill_manual(
        values = c('Up' = '#C74375', 'Down' = '#3070B3'),
        guide  = "none"
      ) +
      theme(
        panel.background = element_blank(),
        panel.grid       = element_blank(),
        axis.line.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.text.y      = element_blank(),
        axis.title.y     = element_blank(),
        axis.line.x      = element_line(color = "black"),
        axis.ticks.x     = element_line(color = "black"),
        axis.text.x      = element_text(size = 12, angle = 0),
        axis.title.x     = element_text(size = 12)
      ),
    
    themes = upset_modify_themes(
      list(
        overall_sizes = theme(
          panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA),
          panel.grid       = element_blank(),
          axis.line.y      = element_blank(),
          axis.ticks.y     = element_blank(),
          axis.text.y      = element_blank(),
          axis.title.y     = element_blank(),
          axis.line.x      = element_line(color = "black"),
          axis.ticks.x     = element_line(color = "black"),
          axis.text.x      = element_text(size = 12),
          axis.title.x     = element_text(size = 12)
        ),
        intersections_matrix = theme(
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid     = element_blank(),
          axis.text.x    = element_blank(),
          axis.ticks.x   = element_blank(),
          axis.title.x   = element_blank()
        ),
        'Intersection size' = theme(
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid       = element_blank(),
          axis.text.x      = element_blank(),
          axis.ticks.x     = element_blank(),
          axis.title.x     = element_blank(),
          axis.text.y      = element_text(size = 12),
          axis.title.y     = element_text(size = 12)
        )
      )
    )
  ) +
    theme(
      plot.title    = element_text(hjust = 0.5, size = 14),
      panel.spacing = unit(0, "lines")
    ) &
    theme(panel.border = element_blank())
}

# ========================
# 构建三个组矩阵 + 预览图
# ========================

cat("--- 构建矩阵并预览 ---\n\n")

mat_T24    <- build_combined_matrix(des_T24_up, ed_T24_up, lim_T24_up,
                                    des_T24_down, ed_T24_down, lim_T24_down)
mat_T32    <- build_combined_matrix(des_T32_up, ed_T32_up, lim_T32_up,
                                    des_T32_down, ed_T32_down, lim_T32_down)
mat_T32v24 <- build_combined_matrix(des_T32v24_up, ed_T32v24_up, lim_T32v24_up,
                                    des_T32v24_down, ed_T32v24_down, lim_T32v24_down)

# 验证DESeq2列不再全FALSE
cat("T24  : DESeq2 TRUE =", sum(mat_T24$DESeq2),    "| edgeR TRUE =", sum(mat_T24$edgeR),    "| limma TRUE =", sum(mat_T24$limma),    "| Total =", nrow(mat_T24), "\n")
cat("T32  : DESeq2 TRUE =", sum(mat_T32$DESeq2),    "| edgeR TRUE =", sum(mat_T32$edgeR),    "| limma TRUE =", sum(mat_T32$limma),    "| Total =", nrow(mat_T32), "\n")
cat("T32v24: DESeq2 TRUE =", sum(mat_T32v24$DESeq2), "| edgeR TRUE =", sum(mat_T32v24$edgeR), "| limma TRUE =", sum(mat_T32v24$limma), "| Total =", nrow(mat_T32v24), "\n\n")

# 预览图（只预览，不保存PDF）
cat(">>> 预览 T24 UpSet 图...\n")
print(plot_upset_combined(mat_T24, "Differentially expressed genes (T24 vs Ctrl)"))

cat(">>> 预览 T32 UpSet 图...\n")
print(plot_upset_combined(mat_T32, "Differentially expressed genes (T32 vs Ctrl)"))

cat(">>> 预览 T32vT24 UpSet 图...\n")
print(plot_upset_combined(mat_T32v24, "Differentially expressed genes (T32 vs T24)"))

# ========================
# 修复 Overlap 文件（数据文件，必须保存）
# ========================

cat("\n--- 修复 Overlap 文件 ---\n")

# T24
write.csv(data.frame(GeneID = highconf_T24$all),  "results/RNAseq_analysis/Figures/DEG/Overlap_T24_All_DEG.csv",  row.names = FALSE)
write.csv(data.frame(GeneID = highconf_T24$up),   "results/RNAseq_analysis/Figures/DEG/Overlap_T24_Up_DEG.csv",   row.names = FALSE)
write.csv(data.frame(GeneID = highconf_T24$down), "results/RNAseq_analysis/Figures/DEG/Overlap_T24_Down_DEG.csv", row.names = FALSE)

# T32
write.csv(data.frame(GeneID = highconf_T32$all),  "results/RNAseq_analysis/Figures/DEG/Overlap_T32_All_DEG.csv",  row.names = FALSE)
write.csv(data.frame(GeneID = highconf_T32$up),   "results/RNAseq_analysis/Figures/DEG/Overlap_T32_Up_DEG.csv",   row.names = FALSE)
write.csv(data.frame(GeneID = highconf_T32$down), "results/RNAseq_analysis/Figures/DEG/Overlap_T32_Down_DEG.csv", row.names = FALSE)

# T32vsT24
write.csv(data.frame(GeneID = highconf_T32v24$all),  "results/RNAseq_analysis/Figures/DEG/Overlap_T32v24_All_DEG.csv",  row.names = FALSE)
write.csv(data.frame(GeneID = highconf_T32v24$up),   "results/RNAseq_analysis/Figures/DEG/Overlap_T32v24_Up_DEG.csv",   row.names = FALSE)
write.csv(data.frame(GeneID = highconf_T32v24$down), "results/RNAseq_analysis/Figures/DEG/Overlap_T32v24_Down_DEG.csv", row.names = FALSE)

cat("T24 All :", length(highconf_T24$all), "\n")
cat("T32 All :", length(highconf_T32$all), "\n")
cat("T32v24 All:", length(highconf_T32v24$all), "\n")

cat("\n========================================\n")
cat("完成。图已预览，Overlap文件已修复。\n")
cat("请确认图中 DESeq2 列是否有黑色圆点。\n")
cat("========================================\n")