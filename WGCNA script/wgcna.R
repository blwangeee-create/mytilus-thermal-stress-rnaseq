# 安装 WGCNA（如果还没装）
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("WGCNA", quietly = TRUE)) BiocManager::install("WGCNA")

library(WGCNA)
library(DESeq2)
library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 1. 加载 VST
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)

cat("VST matrix:", nrow(vst_mat), "genes x", ncol(vst_mat), "samples\n")

# 2. 样本性状表
trait_df <- data.frame(
  row.names = colnames(vst_mat),
  Temperature = as.numeric(gsub(".*_(\\d+).*", "\\1", as.character(vsd$Condition))),
  Ctrl = ifelse(vsd$Condition == "Control_13", 1, 0),
  T24 = ifelse(vsd$Condition == "Treatment_24", 1, 0),
  T32 = ifelse(vsd$Condition == "Treatment_32", 1, 0)
)
cat("\nTrait table:\n")
print(trait_df)

# 3. 数据预处理：转置（WGCNA 要求行=样本，列=基因）
datExpr <- t(vst_mat)

# 过滤低表达/低变异基因（保留 top 5000 高变异）
mads <- apply(datExpr, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:min(5000, length(mads))]
datExpr <- datExpr[, keep]

cat("Genes after filtering (top 5000 by MAD):", ncol(datExpr), "\n")

# 4. 检查缺失值和异常样本
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# 5. 选 soft-thresholding power（预览图）
powers <- c(1:10)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 预览（不保存文件）
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", 
     ylab="Scale Free Topology Model Fit, signed R^2",
     main="Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
     labels=powers, cex=0.9, col="red")
abline(h=0.8, col="red", lty=2)

# 6. 构建网络（小样本保守参数）
# 看上图，R^2 首次 >0.8 的 power 是多少，填到下面
power <- 6  # 如果上图 R^2 在 power=6 时没到 0.8，改成 7 或 8

net <- blockwiseModules(datExpr, power = power,
                        TOMType = "signed", 
                        minModuleSize = 10,
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        verbose = 3)

cat("\nModules detected:", length(table(net$colors)), "\n")
cat("Module sizes:\n")
print(sort(table(net$colors), decreasing = TRUE))

# 7. 模块-性状关联
MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes
moduleTraitCor <- cor(MEs, trait_df, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("\nModule-Trait Correlation:\n")
print(round(moduleTraitCor, 2))
cat("\nP-values:\n")
print(signif(moduleTraitPvalue, 2))

# 8. 提取与 Temperature 显著相关的模块
temp_cor <- moduleTraitCor[, "Temperature"]
sig_modules <- names(temp_cor)[abs(temp_cor) > 0.5 & moduleTraitPvalue[, "Temperature"] < 0.05]
cat("\nModules significantly correlated with Temperature (|r|>0.5, p<0.05):\n")
print(sig_modules)

# 9. 保存结果
dir.create("results/RNAseq_analysis/Results/WGCNA", showWarnings=FALSE, recursive=TRUE)

# 所有模块基因
module_genes <- data.frame(gene = colnames(datExpr), module = net$colors)
write.csv(module_genes, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_genes.csv", row.names=FALSE)

# 显著模块的 Hub 基因（top 30 by kME）
for(mod in sig_modules){
  mod_label <- as.numeric(gsub("ME", "", mod))
  mod_genes <- colnames(datExpr)[net$colors == mod_label]
  
  kME <- abs(cor(datExpr[, mod_genes], MEs[, mod], use="p"))
  hub_genes <- mod_genes[order(kME, decreasing=TRUE)[1:min(30, length(mod_genes))]]
  
  write.csv(data.frame(gene=hub_genes, module=mod_label, kME=sort(kME, decreasing=TRUE)[1:length(hub_genes)]),
            paste0("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_", mod_label, ".csv"), row.names=FALSE)
  cat("Module", mod_label, ":", length(mod_genes), "genes, top Hub saved.\n")
}

# 保存网络
saveRDS(net, "results/RNAseq_analysis/Results/WGCNA/WGCNA_network.rds")
saveRDS(MEs, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_eigengenes.rds")

cat("\nWGCNA complete.\n")



library(WGCNA)
library(DESeq2)
library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 重新加载（如果内存变量还在可跳过）
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)
datExpr <- t(vst_mat)
mads <- apply(datExpr, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:min(5000, length(mads))]
datExpr <- datExpr[, keep]
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 构建网络（power = 8）
power <- 8
cat("Building network with power =", power, "...\n")

net <- blockwiseModules(datExpr, power = power,
                        TOMType = "signed", 
                        minModuleSize = 10,
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        verbose = 3)

cat("\nModules detected:", length(table(net$colors)), "\n")
cat("Module sizes:\n")
print(sort(table(net$colors), decreasing = TRUE))

# 模块-性状关联
trait_df <- data.frame(
  row.names = colnames(vst_mat),
  Temperature = as.numeric(gsub(".*_(\\d+).*", "\\1", as.character(vsd$Condition))),
  Ctrl = ifelse(vsd$Condition == "Control_13", 1, 0),
  T24 = ifelse(vsd$Condition == "Treatment_24", 1, 0),
  T32 = ifelse(vsd$Condition == "Treatment_32", 1, 0)
)

MEs <- moduleEigengenes(datExpr, net$colors)$eigengenes
moduleTraitCor <- cor(MEs, trait_df, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

cat("\nModule-Trait Correlation:\n")
print(round(moduleTraitCor, 2))
cat("\nP-values:\n")
print(signif(moduleTraitPvalue, 2))

# 与 Temperature 显著相关的模块
temp_cor <- moduleTraitCor[, "Temperature"]
sig_modules <- names(temp_cor)[abs(temp_cor) > 0.5 & moduleTraitPvalue[, "Temperature"] < 0.05]
cat("\nModules significantly correlated with Temperature (|r|>0.5, p<0.05):\n")
print(sig_modules)

# 保存
dir.create("results/RNAseq_analysis/Results/WGCNA", showWarnings=FALSE, recursive=TRUE)
module_genes <- data.frame(gene = colnames(datExpr), module = net$colors)
write.csv(module_genes, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_genes.csv", row.names=FALSE)

# Hub 基因
for(mod in sig_modules){
  mod_label <- as.numeric(gsub("ME", "", mod))
  mod_genes <- colnames(datExpr)[net$colors == mod_label]
  kME <- abs(cor(datExpr[, mod_genes], MEs[, mod], use="p"))
  hub_genes <- mod_genes[order(kME, decreasing=TRUE)[1:min(30, length(mod_genes))]]
  write.csv(data.frame(gene=hub_genes, module=mod_label, kME=sort(kME, decreasing=TRUE)[1:length(hub_genes)]),
            paste0("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_", mod_label, ".csv"), row.names=FALSE)
  cat("Module", mod_label, ":", length(mod_genes), "genes, top", length(hub_genes), "Hub saved.\n")
}

saveRDS(net, "results/RNAseq_analysis/Results/WGCNA/WGCNA_network.rds")
saveRDS(MEs, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_eigengenes.rds")
cat("\nWGCNA complete.\n")


# 1. 清除
rm(list = ls())

# 2. 先加载 WGCNA（它的 cor 会覆盖其他包的 cor）
library(WGCNA)

# 3. 再加载 DESeq2（提供 assay 函数）
library(DESeq2)

# 4. 显式指定 WGCNA 的 cor（防止后续冲突）
cor <- WGCNA::cor

# 5. 重新加载数据
setwd("/Volumes/Expansion/NCBI/PRJNA934294")
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)
cat("VST:", nrow(vst_mat), "x", ncol(vst_mat), "\n")

# 6. 准备数据
datExpr <- t(vst_mat)
mads <- apply(datExpr, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:min(5000, length(mads))]
datExpr <- datExpr[, keep]
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 7. 构建网络
power <- 8
cat("Building network with power =", power, "...\n")

net <- blockwiseModules(datExpr, power = power,
                        TOMType = "signed", 
                        minModuleSize = 10,
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        verbose = 3)

cat("\nModules detected:", length(table(net$colors)), "\n")
cat("Module sizes:\n")
print(sort(table(net$colors), decreasing = TRUE))
# 清除之前结果，重新跑
rm(list = ls())
library(WGCNA)
library(DESeq2)
cor <- WGCNA::cor

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 重新加载
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)
datExpr <- t(vst_mat)
mads <- apply(datExpr, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:min(5000, length(mads))]
datExpr <- datExpr[, keep]
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# 调整参数重新跑（关键改动）
power <- 7  # 降低 power，网络更密集，模块更少
cat("Re-building network with power =", power, "...\n")

net <- blockwiseModules(datExpr, power = power,
                        TOMType = "signed", 
                        minModuleSize = 20,      # 提高最小模块大小，过滤碎模块
                        reassignThreshold = 0, 
                        mergeCutHeight = 0.20,   # 放宽合并阈值，相似模块合并
                        numericLabels = TRUE, 
                        pamRespectsDendro = FALSE,
                        verbose = 3)

cat("\nModules detected:", length(table(net$colors)), "\n")
cat("Module sizes:\n")
print(sort(table(net$colors), decreasing = TRUE))


library(WGCNA)
library(DESeq2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")
options(stringsAsFactors = FALSE)

# 加载数据（你已有的）
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)
datExpr <- t(vst_mat)

# 保留 top 5000 高变异基因
mads <- apply(datExpr, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:5000]
datExpr <- datExpr[, keep]

# 检查
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

# ========== 关键修正参数 ==========
power <- 7
net <- blockwiseModules(datExpr, 
                        power = power,
                        TOMType = "signed",
                        minModuleSize = 60,        # ← 提高：过滤小碎片
                        mergeCutHeight = 0.35,     # ← 放宽：积极合并相似模块
                        reassignThreshold = 0,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        verbose = 3)

cat("\nModules detected:", length(table(net$colors)), "\n")
print(sort(table(net$colors), decreasing = TRUE))








library(WGCNA)
library(DESeq2)
library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 确认 net 和 datExpr 在内存中
cat("datExpr:", nrow(datExpr), "samples x", ncol(datExpr), "genes\n")

# ========== 1. 模块颜色 & MEs ==========
moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

cat("Modules (including grey):", length(unique(moduleColors)), "\n")
cat("Module sizes:\n")
print(sort(table(moduleColors), decreasing = TRUE))

# ========== 2. 样本性状表 ==========
traitData <- data.frame(
  row.names = rownames(datExpr),
  Temperature = as.numeric(gsub(".*_(\\d+).*", "\\1", as.character(vsd$Condition))),
  T24 = ifelse(vsd$Condition == "Treatment_24", 1, 0),
  T32 = ifelse(vsd$Condition == "Treatment_32", 1, 0)
)

# ========== 3. 模块-性状关联 ==========
modTraitCor <- cor(MEs, traitData, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nrow(datExpr))

# 预览热图
textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) <- dim(modTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")

# ========== 4. 提取显著模块 ==========
sig_threshold <- 0.05
temp_assoc <- data.frame(
  module = rownames(modTraitCor),
  color = gsub("ME", "", rownames(modTraitCor)),
  cor_Temperature = modTraitCor[, "Temperature"],
  p_Temperature = modTraitP[, "Temperature"],
  cor_T24 = modTraitCor[, "T24"],
  p_T24 = modTraitP[, "T24"],
  cor_T32 = modTraitCor[, "T32"],
  p_T32 = modTraitP[, "T32"]
) %>%
  filter(p_Temperature < sig_threshold | p_T32 < sig_threshold | p_T24 < sig_threshold) %>%
  arrange(p_Temperature)

cat("\n=== Significant modules (p < 0.05) ===\n")
print(temp_assoc)

# ========== 5. 提取 Hub 基因（所有模块，或只保留显著且 >100 基因的模块） ==========
dir.create("results/RNAseq_analysis/Results/WGCNA", showWarnings = FALSE, recursive = TRUE)

# 保存所有模块基因
module_gene_df <- data.frame(
  gene = colnames(datExpr),
  module = moduleColors,
  stringsAsFactors = FALSE
)
write.csv(module_gene_df, "results/RNAseq_analysis/Results/WGCNA/WGCNA_all_module_genes.csv", row.names = FALSE)

# 计算 kME（基因与模块特征向量的相关性）
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# 提取显著模块的 Hub 基因（top 30 by kME），只保留 >100 基因的模块
target_modules <- temp_assoc$color[temp_assoc$color != "grey"]
# 进一步过滤：只保留模块大小 > 100 的
mod_sizes <- table(moduleColors)
target_modules <- target_modules[target_modules %in% names(mod_sizes)[mod_sizes >= 100]]

cat("\nTarget modules for Hub extraction (size >= 100 & significant):\n")
print(target_modules)

for (mod in target_modules) {
  modGenes <- colnames(datExpr)[moduleColors == mod]
  if (length(modGenes) < 10) next
  
  kME <- abs(geneModuleMembership[modGenes, paste0("ME", mod)])
  hubGenes <- modGenes[order(kME, decreasing = TRUE)[1:min(30, length(modGenes))]]
  
  hub_df <- data.frame(
    gene = hubGenes,
    module = mod,
    kME = sort(kME, decreasing = TRUE)[1:length(hubGenes)],
    stringsAsFactors = FALSE
  )
  
  write.csv(hub_df, paste0("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_", mod, ".csv"), row.names = FALSE)
  cat("Module", mod, ":", length(modGenes), "genes, top", nrow(hub_df), "Hub saved\n")
}

# 保存网络对象
saveRDS(net, "results/RNAseq_analysis/Results/WGCNA/WGCNA_network.rds")
saveRDS(MEs, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_eigengenes.rds")
write.csv(temp_assoc, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_trait_correlation.csv", row.names = FALSE)

cat("\nWGCNA complete. Check:\n")
cat("  WGCNA_module_trait_correlation.csv\n")
cat("  WGCNA_hub_*.csv\n")



library(ggplot2)
library(gridExtra)

# 如果 sft 不在内存，重新算（10秒）
if(!exists("sft")){
  powers <- c(1:10)
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
}

# 左图：Scale independence
df1 <- data.frame(
  Power = sft$fitIndices[, 1],
  R2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
)

p1 <- ggplot(df1, aes(x = Power, y = R2)) +
  geom_point(size = 4, color = "#2E86AB") +
  geom_text(aes(label = Power), vjust = -1.2, size = 3.5, color = "#2E86AB") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "#E94F37", linewidth = 0.8) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "#6A4C93", linewidth = 0.8) +
  geom_vline(xintercept = 7, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  annotate("text", x = 7.3, y = 0.85, label = "Selected: power=7", 
           color = "darkgreen", hjust = 0, size = 3.5) +
  labs(title = "Scale Independence", 
       x = "Soft Threshold (power)", 
       y = expression("Signed R"^2)) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

# 右图：Mean connectivity
df2 <- data.frame(
  Power = sft$fitIndices[, 1],
  Connectivity = sft$fitIndices[, 5]
)

p2 <- ggplot(df2, aes(x = Power, y = Connectivity)) +
  geom_point(size = 4, color = "#2E86AB") +
  geom_text(aes(label = Power), vjust = -1.2, size = 3.5, color = "#2E86AB") +
  geom_vline(xintercept = 7, linetype = "dotted", color = "darkgreen", linewidth = 1) +
  labs(title = "Mean Connectivity", 
       x = "Soft Threshold (power)", 
       y = "Mean Connectivity") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

grid.arrange(p1, p2, ncol = 2)


library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

temp_assoc <- read.csv("results/RNAseq_analysis/Results/WGCNA/WGCNA_module_trait_correlation.csv", stringsAsFactors=FALSE)

sig_mods <- temp_assoc %>%
  dplyr::filter(color != "grey") %>%
  dplyr::select(color, cor_Temperature, cor_T24, cor_T32, p_Temperature, p_T24, p_T32)

plot_df <- sig_mods %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Trait", values_to = "Correlation") %>%
  mutate(
    Trait = gsub("cor_", "", Trait),
    pvalue = case_when(
      Trait == "Temperature" ~ sig_mods$p_Temperature[match(color, sig_mods$color)],
      Trait == "T24" ~ sig_mods$p_T24[match(color, sig_mods$color)],
      Trait == "T32" ~ sig_mods$p_T32[match(color, sig_mods$color)]
    ),
    Sig = ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")),
    Label = paste0(round(Correlation, 2), Sig),
    color = factor(color, levels = c("blue", "brown", "turquoise", "cyan"))
  )

# Y 轴标签（英文避免乱码，论文可用）
mod_labels <- c("blue" = "Core heat-stress response\n(n=1089)",
                "brown" = "T24 adaptation\n(n=497)",
                "turquoise" = "Metabolic suppression\n(n=1242)",
                "cyan" = "T24-specific transient\n(n=130)")

p <- ggplot(plot_df, aes(x = Trait, y = color, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#3070B3", mid = "white", high = "#C74375",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_y_discrete(labels = mod_labels) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = NULL, title = "Module-Trait Relationships") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 0, vjust = 1),
        panel.grid = element_blank(),
        legend.position = "right")

print(p)

library(WGCNA)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 1. 恢复 WGCNA 对象
net <- readRDS("results/RNAseq_analysis/Results/WGCNA/WGCNA_network.rds")
MEs <- readRDS("results/RNAseq_analysis/Results/WGCNA/WGCNA_module_eigengenes.rds")

# 2. 恢复 datExpr 和性状表
vsd <- readRDS("results/RNAseq_analysis/Results/vst_object.rds")
vst_mat <- assay(vsd)
datExpr <- t(vst_mat)
mads <- apply(datExpr, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:5000]
datExpr <- datExpr[, keep]
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

traitData <- data.frame(
  row.names = rownames(datExpr),
  Temperature = as.numeric(gsub(".*_(\\d+).*", "\\1", as.character(vsd$Condition))),
  T24 = ifelse(vsd$Condition == "Treatment_24", 1, 0),
  T32 = ifelse(vsd$Condition == "Treatment_32", 1, 0)
)

# 3. 重新计算所有模块的模块-性状关联（不筛选）
moduleColors <- labels2colors(net$colors)
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

modTraitCor <- cor(MEs, traitData, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, nrow(datExpr))

# 4. 保存完整的表格（20个模块）
temp_assoc_all <- data.frame(
  module = rownames(modTraitCor),
  color = gsub("ME", "", rownames(modTraitCor)),
  cor_Temperature = modTraitCor[, "Temperature"],
  p_Temperature = modTraitP[, "Temperature"],
  cor_T24 = modTraitCor[, "T24"],
  p_T24 = modTraitP[, "T24"],
  cor_T32 = modTraitCor[, "T32"],
  p_T32 = modTraitP[, "T32"]
)

write.csv(temp_assoc_all, "results/RNAseq_analysis/Results/WGCNA/WGCNA_module_trait_correlation_ALL.csv", row.names = FALSE)
cat("Saved ALL modules:", nrow(temp_assoc_all), "rows\n")

# 5. 画图：全部20个模块（含grey）
plot_df <- temp_assoc_all %>%
  dplyr::select(color, cor_Temperature, cor_T24, cor_T32, p_Temperature, p_T24, p_T32) %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "Trait", values_to = "Correlation") %>%
  mutate(
    Trait = gsub("cor_", "", Trait),
    pvalue = case_when(
      Trait == "Temperature" ~ temp_assoc_all$p_Temperature[match(color, temp_assoc_all$color)],
      Trait == "T24" ~ temp_assoc_all$p_T24[match(color, temp_assoc_all$color)],
      Trait == "T32" ~ temp_assoc_all$p_T32[match(color, temp_assoc_all$color)]
    ),
    Sig = ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", ""))),
    Label = paste0(round(Correlation, 2), Sig),
    mod_size = as.numeric(table(moduleColors)[match(color, names(table(moduleColors)))]),
    color = factor(color, levels = rev(names(sort(table(moduleColors), decreasing = TRUE))))
  )

p <- ggplot(plot_df, aes(x = Trait, y = color, fill = Correlation)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = Label), size = 3, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_x_discrete(expand = expansion(add = 0.3)) +
  labs(x = NULL, y = "Module", title = "Module-Trait Relationships (All Modules)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        axis.text = element_text(size = 9),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.key.height = unit(1.5, "cm"))

print(p)


library(WGCNA)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 如果 net 不在内存，从 RDS 读取
if(!exists("net")){
  net <- readRDS("results/RNAseq_analysis/Results/WGCNA/WGCNA_network.rds")
}

# 直接取 WGCNA 已算好的聚类树和颜色
geneTree <- net$dendrograms[[1]]
blockGenes <- net$blockGenes[[1]]
moduleColors <- labels2colors(net$colors)
colorsForPlot <- moduleColors[blockGenes]

cat("Tree tips:", length(geneTree$order), "| Colors:", length(colorsForPlot), "\n")

# 预览聚类树图
plotDendroAndColors(
  geneTree,
  colorsForPlot,
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  lwd = 0.5,              # ← 树枝粗细：默认约 1，0.3 很细，0.5 适中
  main = "Gene clustering dendrogram and module assignment"
)


library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 安全提取样本分组
sample_info <- data.frame(
  Sample = rownames(MEs),
  Condition = as.character(vsd$Condition[match(rownames(MEs), colnames(vst_mat))]),
  stringsAsFactors = FALSE
) %>%
  mutate(Temperature = as.numeric(gsub(".*_(\\d+).*", "\\1", Condition)))

eigengene_df <- as.data.frame(MEs) %>%
  mutate(Sample = rownames(MEs)) %>%
  left_join(sample_info, by = "Sample")

eigengene_long <- eigengene_df %>%
  pivot_longer(cols = starts_with("ME"), names_to = "Module", values_to = "Eigengene") %>%
  mutate(Module = gsub("ME", "", Module),
         Temperature = factor(Temperature, levels = c(13, 24, 32)))

# 4 个显著模块（对应 4 条热图）
target_modules <- c("blue", "brown", "turquoise", "cyan")
plot_df <- eigengene_long %>% filter(Module %in% target_modules)

p <- ggplot(plot_df, aes(x = Temperature, y = Eigengene, color = Module, group = Module)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(
    values = c("blue" = "#2E86AB", "brown" = "#8B4513", 
               "turquoise" = "#40E0D0", "cyan" = "#00CED1"),
    labels = c("blue" = "blue (Core heat-stress response, n=1089)",
               "brown" = "brown (T24 adaptation, n=497)",
               "turquoise" = "turquoise (Metabolic suppression, n=1242)",
               "cyan" = "cyan (T24-specific transient, n=130)")
  ) +
  labs(x = "Temperature (°C)", y = "Module Eigengene",
       title = "Module expression trends across temperature",
       color = "Module") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(p)




library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 1. 读取各来源核心基因
hub_blue <- read.csv("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_blue.csv", stringsAsFactors=FALSE)$gene
hub_brown <- read.csv("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_brown.csv", stringsAsFactors=FALSE)$gene

# PPI Hub（注意路径和文件名中的空格）
ppi_t32 <- read.csv("results/PPI/32/T32_string_edges.tsv_MCC_top20 default node.csv", stringsAsFactors=FALSE)
ppi_t24 <- read.csv("results/PPI/24/T24_string_edges.tsv_MCC_top20 default node.csv", stringsAsFactors=FALSE)

# 查看列名，确认基因名在哪一列
cat("T32 columns:", colnames(ppi_t32), "\n")
cat("T24 columns:", colnames(ppi_t24), "\n")

# 取前20（根据实际列名调整，通常是 name 或 node）
ppi_t32_genes <- ppi_t32[[1]][1:20]  # 第一列通常是基因名
ppi_t24_genes <- ppi_t24[[1]][1:20]

k3 <- read.csv("results/RNAseq_analysis/Results/kmeans/K3_cluster_3.csv", stringsAsFactors=FALSE)$gene_id
k1 <- read.csv("results/RNAseq_analysis/Results/kmeans/K3_cluster_1.csv", stringsAsFactors=FALSE)$gene_id

# 2. 取交集
core_heat <- Reduce(intersect, list(hub_blue, ppi_t32_genes, k3))
core_t24  <- Reduce(intersect, list(hub_brown, ppi_t24_genes, k1))

cat("=== Core heat-stress candidates (blue ∩ PPI_T32 ∩ Cluster3) ===\n")
print(core_heat)

cat("\n=== T24 adaptation candidates (brown ∩ PPI_T24 ∩ Cluster1) ===\n")
print(core_t24)

# 3. 保存
write.csv(data.frame(gene=core_heat, category="Core_heat_stress"), 
          "results/RNAseq_analysis/Results/core_candidates_heat.csv", row.names=FALSE)
write.csv(data.frame(gene=core_t24, category="T24_adaptation"), 
          "results/RNAseq_analysis/Results/core_candidates_T24.csv", row.names=FALSE)

library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 手动读列名，避免 # 干扰
header_line <- readLines("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", n = 1)
col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]

feat <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", 
                   stringsAsFactors = FALSE, skip = 1, header = FALSE, 
                   col.names = col_names)

cat("列名:", colnames(feat), "\n")
print(head(feat, 3))





library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 1. MCOR_ID ↔ 蛋白ID ==========
header_line <- readLines("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", n = 1)
col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]

feat <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", 
                   stringsAsFactors = FALSE, skip = 1, header = FALSE, 
                   col.names = col_names)

mrna_map <- feat %>%
  dplyr::filter(feature == "mRNA", related_accession != "") %>%
  dplyr::select(locus_tag, protein_id = related_accession) %>%
  dplyr::distinct()

cat("Feature table mapping:", nrow(mrna_map), "rows\n")
print(head(mrna_map, 5))

# ========== 2. 蛋白ID ↔ symbol (eggNOG) ==========
eggnog <- read.delim("annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular", 
                     stringsAsFactors = FALSE, comment.char = "")

# 清理列名：去掉 # 前缀，把 X.query 改回 query
names(eggnog) <- gsub("^#", "", names(eggnog))
names(eggnog) <- gsub("^X\\.", "", names(eggnog))  # 关键：X.query → query

cat("eggNOG columns:", names(eggnog)[1:5], "...\n")

symbol_map <- eggnog %>%
  dplyr::select(query, Preferred_name) %>%
  dplyr::filter(Preferred_name != "" & Preferred_name != "-") %>%
  dplyr::distinct()

cat("eggNOG symbol mapping:", nrow(symbol_map), "rows\n")
print(head(symbol_map, 5))

# ========== 3. 合并：gene-MCOR_... ↔ symbol ==========
id_map <- mrna_map %>%
  dplyr::left_join(symbol_map, by = c("protein_id" = "query")) %>%
  dplyr::filter(!is.na(Preferred_name)) %>%
  dplyr::mutate(gene_id = paste0("gene-", locus_tag)) %>%
  dplyr::select(gene_id, symbol = Preferred_name) %>%
  dplyr::distinct()

cat("\nFinal ID map (MCOR → symbol):", nrow(id_map), "genes\n")
print(head(id_map, 10))

# ========== 4. WGCNA Hub / K-means 转 symbol ==========
hub_blue <- read.csv("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_blue.csv", stringsAsFactors=FALSE)
hub_brown <- read.csv("results/RNAseq_analysis/Results/WGCNA/WGCNA_hub_brown.csv", stringsAsFactors=FALSE)
k1 <- read.csv("results/RNAseq_analysis/Results/kmeans/K3_cluster_1.csv", stringsAsFactors=FALSE)
k3 <- read.csv("results/RNAseq_analysis/Results/kmeans/K3_cluster_3.csv", stringsAsFactors=FALSE)

hub_blue_sym <- id_map$symbol[id_map$gene_id %in% hub_blue$gene]
hub_brown_sym <- id_map$symbol[id_map$gene_id %in% hub_brown$gene]
k1_sym <- id_map$symbol[id_map$gene_id %in% k1$gene_id]
k3_sym <- id_map$symbol[id_map$gene_id %in% k3$gene_id]

hub_blue_sym <- hub_blue_sym[!is.na(hub_blue_sym) & hub_blue_sym != ""]
hub_brown_sym <- hub_brown_sym[!is.na(hub_brown_sym) & hub_brown_sym != ""]
k1_sym <- k1_sym[!is.na(k1_sym) & k1_sym != ""]
k3_sym <- k3_sym[!is.na(k3_sym) & k3_sym != ""]

cat("\nWGCNA blue Hub with symbol:", length(hub_blue_sym), "\n")
cat("WGCNA brown Hub with symbol:", length(hub_brown_sym), "\n")
cat("K-means Cluster1 with symbol:", length(k1_sym), "\n")
cat("K-means Cluster3 with symbol:", length(k3_sym), "\n")

# ========== 5. PPI Hub (symbol) ==========
ppi_t32 <- read.csv("results/PPI/32/T32_string_edges.tsv_MCC_top20 default node.csv", stringsAsFactors=FALSE)
ppi_t24 <- read.csv("results/PPI/24/T24_string_edges.tsv_MCC_top20 default node.csv", stringsAsFactors=FALSE)

ppi_t32_sym <- ppi_t32$name[1:20]
ppi_t24_sym <- ppi_t24$name[1:20]

# ========== 6. 三重交集 ==========
core_heat_3way <- Reduce(intersect, list(hub_blue_sym, ppi_t32_sym, k3_sym))
core_t24_3way  <- Reduce(intersect, list(hub_brown_sym, ppi_t24_sym, k1_sym))

cat("\n=== Core heat-stress (blue ∩ PPI_T32 ∩ Cluster3) ===\n")
cat("Count:", length(core_heat_3way), "\n")
print(core_heat_3way)

cat("\n=== T24 adaptation (brown ∩ PPI_T24 ∩ Cluster1) ===\n")
cat("Count:", length(core_t24_3way), "\n")
print(core_t24_3way)

# ========== 7. 保存 ==========
if(length(core_heat_3way) > 0) {
  write.csv(data.frame(gene = core_heat_3way, category = "Core_heat_stress_3way"), 
            "results/RNAseq_analysis/Results/core_candidates_heat_3way.csv", row.names = FALSE)
}
if(length(core_t24_3way) > 0) {
  write.csv(data.frame(gene = core_t24_3way, category = "T24_adaptation_3way"), 
            "results/RNAseq_analysis/Results/core_candidates_T24_3way.csv", row.names = FALSE)
}