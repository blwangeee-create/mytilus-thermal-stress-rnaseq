library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 1. 读取 VST
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)
cat("VST loaded:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

# 2. 预览 VST QC 三张图（不保存文件）
par(mfrow = c(1, 2))
boxplot(vst_mat, las = 2, main = "VST Boxplot", outline = FALSE,
        col = c(rep("blue", 3), rep("green", 3), rep("red", 3)))
plotPCA(vsd, intgroup = "Condition") + ggtitle("PCA after VST")
par(mfrow = c(1, 1))

# 3. 加载 HighConf 基因
load("results/RNAseq_analysis/Results/Phase1_DEG_analysis.RData")
target_genes <- unique(c(highconf_T24$all, highconf_T32$all))
target_genes <- target_genes[target_genes %in% rownames(vst_mat)]
cat("Union HighConf genes:", length(target_genes), "\n")

# 4. 按组求平均表达
expr <- vst_mat[target_genes, ]
sample_df <- data.frame(
  sample = colnames(expr),
  group = as.character(vsd$Condition),
  stringsAsFactors = FALSE
)

# 转换为长格式，求组内均值，再转回宽格式（确保列顺序）
expr_long <- expr %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expression") %>%
  left_join(sample_df, by = "sample") %>%
  group_by(gene_id, group) %>%
  summarise(mean_expr = mean(expression), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_expr) %>%
  as.data.frame()

rownames(expr_long) <- expr_long$gene_id
expr_long$gene_id <- NULL

# 固定列顺序：Control_13, Treatment_24, Treatment_32
expr_mean <- expr_long[, c("Control_13", "Treatment_24", "Treatment_32")]
stopifnot(all(c("Control_13", "Treatment_24", "Treatment_32") %in% colnames(expr_mean)))
cat("Group means computed. Dimensions:", nrow(expr_mean), "genes x 3 groups\n")

# 5. Z-score 标准化（按基因）
expr_scaled <- t(apply(expr_mean, 1, scale))
colnames(expr_scaled) <- c("Ctrl", "T24", "T32")
expr_scaled <- expr_scaled[complete.cases(expr_scaled), ]
cat("Genes after scaling:", nrow(expr_scaled), "\n")

# 6. k-means 聚类 (K=3,4,5) 并预览趋势图
set.seed(123)
for (k in 3:5) {
  km <- kmeans(expr_scaled, centers = k, nstart = 25)
  
  # 构建簇平均趋势数据
  cluster_means <- data.frame(expr_scaled, cluster = km$cluster) %>%
    pivot_longer(cols = c("Ctrl", "T24", "T32"), names_to = "temperature", values_to = "zscore") %>%
    mutate(temperature = factor(temperature, levels = c("Ctrl", "T24", "T32"))) %>%
    group_by(cluster, temperature) %>%
    summarise(mean_zscore = mean(zscore), .groups = "drop")
  
  # 预览趋势图
  p <- ggplot(cluster_means, aes(x = temperature, y = mean_zscore, group = cluster, color = factor(cluster))) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    labs(title = paste("k-means (K =", k, ")"), subtitle = paste("n =", nrow(expr_scaled), "genes"),
         x = "Temperature", y = "Mean Z-score") +
    theme_bw() +
    scale_color_brewer(palette = "Set1", name = "Cluster")
  print(p)
  
  cat("\nK =", k, "| Cluster sizes:\n")
  print(table(km$cluster))
}
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(123)
km3 <- kmeans(expr_scaled, centers = 3, nstart = 25)

# 每个基因的折线
plot_df <- data.frame(expr_scaled, gene_id = rownames(expr_scaled), cluster = km3$cluster) %>%
  pivot_longer(cols = c("Ctrl", "T24", "T32"), names_to = "temperature", values_to = "zscore") %>%
  mutate(temperature = factor(temperature, levels = c("Ctrl", "T24", "T32")))

# 平均折线
mean_df <- plot_df %>%
  group_by(cluster, temperature) %>%
  summarise(mean_zscore = mean(zscore), .groups = "drop")

# 同一张图
p <- ggplot() +
  # 半透明细线（背景基因轨迹）
  geom_line(data = plot_df, 
            aes(x = temperature, y = zscore, group = gene_id, color = factor(cluster)), 
            alpha = 0.05, linewidth = 0.1, linetype = "solid") +
  
  # 平均折线（粗线 + 虚线区分）
  geom_line(data = mean_df, 
            aes(x = temperature, y = mean_zscore, group = factor(cluster), 
                color = factor(cluster), linetype = factor(cluster)), 
            linewidth = 0.5) +  # ← 改粗细
  
  # 节点（改成方块或去掉）
  geom_point(data = mean_df, 
             aes(x = temperature, y = mean_zscore, color = factor(cluster)), 
             size = 4, shape = NA) +  # ← 15=方块, 16=圆, 17=三角, NA=去掉
  
  scale_linetype_manual(values = c("1" = "dashed", "2" = "dashed", "3" = "dashed"),
                        name = "Cluster") +
  
  scale_x_discrete(expand = expansion(add = 0.05)) +
  
  scale_color_manual(values = c("1" = "#C74375", "2" = "#FFCD00", "3" = "#0B8BEE"),
                     labels = c("1" = "T24 峰值响应 (n=273)", 
                                "2" = "高温抑制 (n=465)", 
                                "3" = "热胁迫持续响应 (n=571)"),
                     name = "Cluster") +
  
  labs(title = "k-means Clustering (K = 3)",
       subtitle = paste("n =", nrow(expr_scaled), "genes"),
       x = "Temperature", y = "Z-score") +
  
  theme_bw() +
  theme(legend.position = "top")

print(p)


library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 1. 重新加载数据
load("./results/RNAseq_analysis/Results/Phase1_DEG_analysis.RData")
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)

# 2. 重新跑 k-means K=3
target_genes <- unique(c(highconf_T24$all, highconf_T32$all))
target_genes <- target_genes[target_genes %in% rownames(vst_mat)]

expr <- vst_mat[target_genes, ]

sample_df <- data.frame(sample = colnames(expr), group = as.character(vsd$Condition))

expr_mean <- as.data.frame(expr)
expr_mean$gene_id <- rownames(expr_mean)

expr_long <- expr_mean %>%
  pivot_longer(cols=-gene_id, names_to="sample", values_to="expression") %>%
  left_join(sample_df, by="sample") %>%
  group_by(gene_id, group) %>%
  summarise(mean_expr=mean(expression), .groups="drop")

expr_wide <- expr_long %>%
  pivot_wider(names_from=group, values_from=mean_expr) %>%
  as.data.frame()

rownames(expr_wide) <- expr_wide$gene_id
expr_wide$gene_id <- NULL
expr_wide <- expr_wide[, c("Control_13", "Treatment_24", "Treatment_32")]

expr_scaled <- t(apply(expr_wide, 1, scale))
colnames(expr_scaled) <- c("Ctrl", "T24", "T32")
expr_scaled <- expr_scaled[complete.cases(expr_scaled), ]

set.seed(123)
km3 <- kmeans(expr_scaled, centers=3, nstart=25)

cat("K=3 Cluster sizes:\n")
print(table(km3$cluster))

# 3. 保存簇基因
dir.create("results/RNAseq_analysis/Results/kmeans", showWarnings=FALSE, recursive=TRUE)

clusters <- list()
for(i in 1:3){
  genes <- rownames(expr_scaled)[km3$cluster == i]
  clusters[[i]] <- genes
  write.csv(data.frame(gene_id=genes, cluster=i), 
            paste0("results/RNAseq_analysis/Results/kmeans/K3_cluster_", i, ".csv"), 
            row.names=FALSE)
  cat("Cluster", i, ":", length(genes), "genes saved\n")
}

# 4. 读取注释
gene_annot <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors=FALSE)

# 5. 富集函数
run_enrich <- function(gene_ids, label){
  symbols <- gene_annot %>%
    filter(gene_id %in% gene_ids, Preferred_name != "", !is.na(Preferred_name)) %>%
    distinct(Preferred_name) %>%
    pull(Preferred_name)
  
  cat("\n===", label, "===\n")
  cat("Genes:", length(gene_ids), "| Symbols:", length(symbols), "\n")
  
  if(length(symbols) < 10){
    cat("Too few symbols, skipping enrichment\n")
    return(NULL)
  }
  
  # GO (BP)
  go <- enrichGO(gene = symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  
  # KEGG
  kegg <- enrichKEGG(gene = symbols, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  
  if(!is.null(go) && nrow(as.data.frame(go)) > 0){
    write.csv(as.data.frame(go), paste0("results/RNAseq_analysis/Results/kmeans/K3_", label, "_GO.csv"), row.names=FALSE)
    cat("GO top 5:\n")
    print(head(as.data.frame(go)[, c("Description", "p.adjust", "Count")], 5))
  } else {
    cat("GO: No significant enrichment\n")
  }
  
  if(!is.null(kegg) && nrow(as.data.frame(kegg)) > 0){
    write.csv(as.data.frame(kegg), paste0("results/RNAseq_analysis/Results/kmeans/K3_", label, "_KEGG.csv"), row.names=FALSE)
    cat("KEGG top 5:\n")
    print(head(as.data.frame(kegg)[, c("Description", "p.adjust", "Count")], 5))
  } else {
    cat("KEGG: No significant enrichment\n")
  }
}

# 6. 跑三个簇
for(i in 1:3){
  run_enrich(clusters[[i]], paste0("Cluster", i))
}

cat("\nAll done. Results saved to results/RNAseq_analysis/Results/kmeans/\n")




library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取本地 eggNOG KEGG 注释
kegg_term2gene <- readRDS("eggnog_annotation/kegg_term2gene.rds")
kegg_pathway <- readRDS("eggnog_annotation/kegg_term2gene_pathway.rds")

# 读取注释和簇
gene_annot <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors=FALSE)

clusters <- list()
for(i in 1:3){
  clusters[[i]] <- read.csv(paste0("results/RNAseq_analysis/Results/kmeans/K3_cluster_", i, ".csv"), stringsAsFactors=FALSE)$gene_id
}

# 富集函数（GO + 本地 KEGG）
run_enrich <- function(gene_ids, label){
  symbols <- gene_annot %>%
    filter(gene_id %in% gene_ids, Preferred_name != "", !is.na(Preferred_name)) %>%
    distinct(Preferred_name) %>%
    pull(Preferred_name)
  
  cat("\n===", label, "===\n")
  cat("Genes:", length(gene_ids), "| Symbols:", length(symbols), "\n")
  
  if(length(symbols) < 10) return(NULL)
  
  # GO (BP) - 本地数据库，不会超时
  go <- enrichGO(gene = symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                 ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  
  if(!is.null(go) && nrow(as.data.frame(go)) > 0){
    write.csv(as.data.frame(go), paste0("results/RNAseq_analysis/Results/kmeans/K3_", label, "_GO.csv"), row.names=FALSE)
    cat("GO top 5:\n")
    print(head(as.data.frame(go)[, c("Description", "p.adjust", "Count")], 5))
  } else {
    cat("GO: No significant enrichment\n")
  }
  
  # KEGG - 用本地 eggNOG 数据，避免网络超时
  # 把 symbol 映射到 kegg 的 gene
  mapped_kegg <- kegg_term2gene[kegg_term2gene$gene %in% symbols, ]
  
  if(nrow(mapped_kegg) > 0){
    # 构建 TERM2GENE
    term2gene <- mapped_kegg[, c("term", "gene")]
    # 去重
    term2gene <- distinct(term2gene)
    
    # 用 enricher 做超几何检验
    kegg_res <- enricher(gene = symbols, TERM2GENE = term2gene, pAdjustMethod = "BH", pvalueCutoff = 0.05)
    
    if(!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0){
      # 添加 pathway 名称
      kegg_df <- as.data.frame(kegg_res)
      kegg_df$Description <- kegg_pathway$pathway[match(kegg_df$ID, kegg_pathway$term)]
      write.csv(kegg_df, paste0("results/RNAseq_analysis/Results/kmeans/K3_", label, "_KEGG.csv"), row.names=FALSE)
      cat("KEGG top 5:\n")
      print(head(kegg_df[, c("Description", "p.adjust", "Count")], 5))
    } else {
      cat("KEGG: No significant enrichment\n")
    }
  } else {
    cat("KEGG: No mapped genes\n")
  }
}

# 跑三个簇
for(i in 1:3){
  run_enrich(clusters[[i]], paste0("Cluster", i))
}

cat("\nAll done.\n")

library(dplyr)
library(ggplot2)
library(ggprism)
library(gground)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 1. 读取并合并
go_list <- list()
for(i in 1:3) {
  f <- paste0("results/RNAseq_analysis/Results/kmeans/K3_Cluster", i, "_GO.csv")
  if(file.exists(f)) {
    df <- read.csv(f, stringsAsFactors = FALSE) %>%
      mutate(Cluster = paste0("Cluster", i))
    go_list[[i]] <- df
  }
}
go_all <- bind_rows(go_list) %>%
  filter(p.adjust < 0.05) %>%
  mutate(Cluster = factor(Cluster, levels = paste0("Cluster", 1:3)))

# 2. 每个簇选 Top5（不去重）
top_go <- go_all %>%
  group_by(Cluster) %>%
  arrange(p.adjust, desc(Count)) %>%
  slice_head(n = 5) %>%
  ungroup()

# 3. 生成 geneID_short
top_go <- top_go %>%
  mutate(geneID_short = sapply(strsplit(geneID, "/"), function(x) {
    paste(head(x, 3), collapse = ", ")
  }))

# 4. 关键：每个 cluster 内按显著性降序排列（最显著在最上）
# 先按 Cluster 和 p.adjust 排序，然后反转 levels（coord_flip 后最上）
top_go <- top_go %>%
  arrange(Cluster, p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# 5. 左侧标记数据
width <- 0.5
xaxis_max <- max(-log10(top_go$p.adjust)) + 1
rect_data <- top_go %>%
  group_by(Cluster) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  mutate(
    xmin = -3 * width,
    xmax = -2 * width,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

# 6. 画图
pal_cluster <- c('#C74375', '#FFCD00', '#0B8BEE')

p <- ggplot(top_go, aes(-log10(p.adjust), y = Description, fill = Cluster)) +
  geom_round_col(width = 0.6, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 4) +
  geom_text(aes(x = 0.2, label = geneID_short, colour = Cluster),
            hjust = 0, vjust = 2.8, size = 2.8, fontface = "italic") +
  geom_point(aes(x = -width, size = Count), shape = 21, fill = "white", stroke = 1.2) +
  geom_text(aes(x = -width, label = Count), size = 3) +
  scale_size_continuous(name = 'Count', range = c(4, 10)) +
  scale_fill_manual(values = pal_cluster) +
  scale_colour_manual(values = pal_cluster) +
  geom_round_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Cluster),
                  data = rect_data, radius = unit(2, 'mm'), inherit.aes = FALSE) +
  geom_text(aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2, label = Cluster),
            data = rect_data, inherit.aes = FALSE, size = 4, color = "white") +
  annotate("segment", x = 0, y = 0, xend = xaxis_max, yend = 0, linewidth = 1.5) +
  labs(y = NULL, x = expression(-log[10](italic(P)[adj]))) +
  scale_x_continuous(breaks = seq(0, ceiling(xaxis_max), 2),
                     expand = expansion(c(0, 0))) +
  theme_prism() +
  theme(axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "top")

print(p)