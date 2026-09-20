# ============================================================
# WGCNA 图表重绘脚本(power = 18 主分析版)
# 说明: 绘图代码完全沿用原 wgcna.R 的风格, 仅替换:
#   ① 数据源: WGCNA_corrected/power18/ (power=18 + ANOVA FDR)
#   ② 星号依据: corPvalueStudent(n=9) → ANOVA FDR
#   ③ 模块名/大小: 按 power=18 结果(7 模块)
# 输出: results/RNAseq_analysis/Figures/wgcna_power18/
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

BASE <- "/Volumes/Expansion/NCBI/PRJNA934294"
PP18 <- file.path(BASE, "results/RNAseq_analysis/Results/WGCNA_corrected/power18")
OUT  <- file.path(BASE, "results/RNAseq_analysis/Figures/wgcna_power18")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ============ 预览开关 ============
# TRUE  = 只在 RStudio Plots 面板显示,不写文件(用于调整比例)
# FALSE = 正常输出 PDF 文件
PREVIEW <- TRUE

# 输出辅助: 预览模式只 print; 否则存 PDF
out_fig <- function(plot, filename, w, h, dpi = 300) {
  if (PREVIEW) {
    print(plot)
    cat("[预览]", filename, "-> Plots 面板\n")
  } else {
    ggsave(file.path(OUT, filename), plot, width = w, height = h, dpi = dpi)
    cat("✓ 已保存:", filename, "\n")
  }
}

# ---------- 1. 读入 power=18 的模块-性状表 ----------
temp_assoc <- read.csv(file.path(PP18, "module_trait_full.csv"), stringsAsFactors = FALSE)
cat("读到模块数:", nrow(temp_assoc), "\n")
print(temp_assoc[, c("color", "n_genes", "cor_Temperature", "fdr_Temperature")])

# ---------- 2. 长表(完全沿用原 pivot 逻辑; 三性状各自的 FDR) ----------
sig_mods <- temp_assoc %>%
  dplyr::filter(color != "grey") %>%
  dplyr::select(color, n_genes, cor_Temperature, cor_T24, cor_T32,
                fdr_Temperature, fdr_T24, fdr_T32)

plot_df <- sig_mods %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Trait", values_to = "Correlation") %>%
  mutate(
    Trait = gsub("cor_", "", Trait),
    fdr = mapply(function(tr, col) {
      sig_mods[[paste0("fdr_", tr)]][match(col, sig_mods$color)]
    }, Trait, color),
    Sig = ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", ""))),
    Label = paste0(round(Correlation, 2), Sig),
    color = factor(color, levels = rev(sig_mods$color[order(sig_mods$fdr_Temperature)]))
  )

# ---------- 3. Y 轴标签(power=18 的实际模块与大小) ----------
mod_labels <- setNames(
  paste0(sig_mods$color, "\n(n=", sig_mods$n_genes, ")"),
  sig_mods$color
)
cat("\nY 轴标签:\n"); print(mod_labels)

# ---------- 4. 模块-性状热图(沿用原样式) ----------
p <- ggplot(plot_df, aes(x = Trait, y = color, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#3070B3", mid = "white", high = "#C74375",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_y_discrete(labels = mod_labels) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = NULL, title = "Module-Trait Relationships (power = 18)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 0, vjust = 1),
        panel.grid = element_blank(),
        legend.position = "right")

out_fig(p, "wgcna heatmap s.pdf", 6, 5)

# ---------- 5. 全模块热图(含 grey, 沿用原样式) ----------
plot_df_all <- temp_assoc %>%
  dplyr::select(color, n_genes, cor_Temperature, cor_T24, cor_T32,
                fdr_Temperature, fdr_T24, fdr_T32) %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Trait", values_to = "Correlation") %>%
  mutate(
    Trait = gsub("cor_", "", Trait),
    fdr = mapply(function(tr, col) {
      temp_assoc[[paste0("fdr_", tr)]][match(col, temp_assoc$color)]
    }, Trait, color),
    Sig = ifelse(fdr < 0.001, "***", ifelse(fdr < 0.01, "**", ifelse(fdr < 0.05, "*", ""))),
    Label = paste0(round(Correlation, 2), Sig),
    color = factor(color, levels = rev(temp_assoc$color[order(temp_assoc$fdr_Temperature)]))
  )

p_all <- ggplot(plot_df_all, aes(x = Trait, y = color, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = Label), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_x_discrete(expand = expansion(add = 0.3)) +
  labs(x = NULL, y = "Module", title = "Module-Trait Relationships (All Modules, power = 18)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.text = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.height = unit(1.5, "cm"))

out_fig(p_all, "wgcna heatmap.pdf", 6, 5)

# ---------- 6. 软阈值图(沿用原样式, 标出 power=18) ----------
sft <- read.csv(file.path(PP18, "soft_threshold.csv"), stringsAsFactors = FALSE)
df1 <- data.frame(Power = sft$power, R2 = sft$signedR2)

p1 <- ggplot(df1, aes(x = Power, y = R2)) +
  geom_point(size = 4, color = "#2E86AB") +
  geom_text(aes(label = Power), vjust = -1.2, size = 3.5, color = "#2E86AB") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#E94F37", linewidth = 0.8) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "#6A4C93", linewidth = 0.8) +
  geom_vline(xintercept = 18, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  annotate("text", x = 18, y = 0.85, label = "Selected: power=18",
           color = "darkgreen", hjust = 0.5, size = 3.5) +
  labs(title = "Scale Independence", x = "Soft Threshold (power)",
       y = expression("Signed R"^2)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

df2 <- data.frame(Power = sft$power, Connectivity = sft$mean_k)
p2 <- ggplot(df2, aes(x = Power, y = Connectivity)) +
  geom_point(size = 4, color = "#2E86AB") +
  geom_text(aes(label = Power), vjust = -1.2, size = 3.5, color = "#2E86AB") +
  geom_vline(xintercept = 18, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  labs(title = "Mean Connectivity", x = "Soft Threshold (power)", y = "Mean Connectivity") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

if (PREVIEW) {
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  cat("[预览] Rplot.pdf (软阈值图) -> Plots 面板\n")
} else {
  pdf(file.path(OUT, "Rplot.pdf"), width = 10, height = 4.5)
  gridExtra::grid.arrange(p1, p2, ncol = 2)
  dev.off()
  cat("✓ 已保存: Rplot.pdf (软阈值图)\n")
}

if (!PREVIEW) cat("\n全部图输出至:", OUT, "\n")

# ============================================================
# 7. ME 趋势图(Eigengene, 沿用原样式, 数据换 power=18)
# ============================================================
MEs <- readRDS(file.path(PP18, "module_eigengenes.rds"))

# 样本分组(顺序: SRR23438720-722=T32, 723-725=T24, 726-728=Ctrl13)
sample_info <- data.frame(
  Sample = rownames(MEs),
  Temperature = c(32,32,32,24,24,24,13,13,13),
  stringsAsFactors = FALSE
)

eigengene_df <- as.data.frame(MEs) %>%
  mutate(Sample = rownames(MEs)) %>%
  left_join(sample_info, by = "Sample")

eigengene_long <- eigengene_df %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene") %>%
  mutate(Module = gsub("ME", "", Module),
         Temperature = factor(Temperature, levels = c(13, 24, 32)))

# power=18 的显著模块(turquoise/blue/red)
target_modules <- c("turquoise", "blue", "red")
plot_df <- eigengene_long %>% filter(Module %in% target_modules)

p_eig <- ggplot(plot_df, aes(x = Temperature, y = Eigengene, color = Module, group = Module)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(
    values = c("turquoise" = "#40E0D0", "blue" = "#2E86AB", "red" = "#E94F37"),
    labels = c("turquoise" = "turquoise (Core heat-stress response, n=2024)",
               "blue" = "blue (Down-regulated by heat, n=1015)",
               "red" = "red (Transient at 24C, n=151)")
  ) +
  labs(x = "Temperature (°C)", y = "Module Eigengene",
       title = "Module expression trends across temperature (power = 18)",
       color = "Module") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

out_fig(p_eig, "Eigengene.pdf", 6, 4.5)

# ============================================================
# 8. 基因聚类树 + 模块颜色条(沿用原样式, 数据换 power=18)
# ============================================================
suppressPackageStartupMessages(library(WGCNA))
net18 <- readRDS(file.path(PP18, "network.rds"))
geneTree <- net18$dendrograms[[1]]
blockGenes <- net18$blockGenes[[1]]
moduleColors <- labels2colors(net18$colors)
colorsForPlot <- moduleColors[blockGenes]

cat("Tree tips:", length(geneTree$order), "| Colors:", length(colorsForPlot), "\n")

draw_tree <- function() {
  plotDendroAndColors(
    geneTree,
    colorsForPlot,
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    lwd = 0.5,
    main = "Gene clustering dendrogram and module assignment (power = 18)"
  )
}

if (PREVIEW) {
  draw_tree()
  cat("[预览] wgcna tree.pdf -> Plots 面板\n")
} else {
  pdf(file.path(OUT, "wgcna tree.pdf"), width = 9, height = 5)
  draw_tree()
  dev.off()
  cat("✓ 已保存: wgcna tree.pdf\n")
}

if (!PREVIEW) cat("\n全部图(5 张)输出至:", OUT, "\n")
