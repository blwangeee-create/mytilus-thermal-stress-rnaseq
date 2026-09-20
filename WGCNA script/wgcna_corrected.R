# ============================================================
# WGCNA 修正版脚本
# 修正内容:
#   ① 模块-性状显著性: corPvalueStudent(n=9, 伪重复) → ANOVA(3组×3重复, df=2,6) + BH FDR
#   ② 软阈值: 同时跑 power=7(原参数) 与 power=18(n<20 signed 网络推荐值) 做稳健性对比
# 输入: results/RNAseq_analysis/Results/vst_expression_matrix.csv (9 样本 × 20,030 基因)
# 输出: results/RNAseq_analysis/Results/WGCNA_corrected/{power7,power18}/
# 说明: 不覆盖任何原始文件; 原脚本与结果已备份至 *_backup_20260919
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  enableWGCNAThreads()
})

BASE <- "/Volumes/Expansion/NCBI/PRJNA934294"
OUT  <- file.path(BASE, "results/RNAseq_analysis/Results/WGCNA_corrected")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---------- 1. 读入 VST 矩阵 ----------
cat("=== 读取 VST 表达矩阵 ===\n")
vst <- read.csv(file.path(BASE, "results/RNAseq_analysis/Results/vst_expression_matrix.csv"),
                row.names = 1, check.names = FALSE)
cat("矩阵维度:", nrow(vst), "基因 x", ncol(vst), "样本\n")
cat("样本名:", colnames(vst), "\n")

# 样本顺序即 SRR23438720..728 = T32×3, T24×3, Ctrl13×3 (与原始脚本一致)
sample_names <- colnames(vst)
# 组别判定: 前3个=T32, 中间3个=T24, 后3个=Ctrl13 (依据原始脚本 trait 表)
grp <- factor(c(rep("T32", 3), rep("T24", 3), rep("Ctrl13", 3)),
              levels = c("Ctrl13", "T24", "T32"))
names(grp) <- sample_names
cat("\n分组:\n"); print(grp)

# ---------- 2. 转置 + top 5000 MAD ----------
datExpr_full <- t(vst)                       # 行=样本, 列=基因
mads <- apply(datExpr_full, 2, mad)
keep <- order(mads, decreasing = TRUE)[1:min(5000, length(mads))]
datExpr <- datExpr_full[, keep]
cat("\n用于建网:", nrow(datExpr), "样本 x", ncol(datExpr), "基因\n")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  grp <- grp[rownames(datExpr)]
}

# ---------- 3. 性状表(用于 Pearson 相关展示; 显著性改用 ANOVA) ----------
traitData <- data.frame(
  row.names   = rownames(datExpr),
  Temperature = as.numeric(sub("Ctrl13", "13", sub("T", "", as.character(grp)))),
  T24         = as.integer(grp == "T24"),
  T32         = as.integer(grp == "T32")
)

# ---------- 4. 核心分析函数 ----------
run_wgcna <- function(power_val, minModSize, mergeHeight, tag) {

  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("运行 WGCNA: power =", power_val, "| minModuleSize =", minModSize,
      "| mergeCutHeight =", mergeHeight, "\n")
  cat(strrep("=", 60), "\n", sep = "")

  net <- blockwiseModules(datExpr, power = power_val,
                          TOMType = "signed",
                          minModuleSize = minModSize,
                          reassignThreshold = 0,
                          mergeCutHeight = mergeHeight,
                          numericLabels = TRUE,
                          pamRespectsDendro = FALSE,
                          verbose = 3)

  moduleColors <- labels2colors(net$colors)
  n_mods <- length(unique(moduleColors))
  cat("\n检测到模块数(含 grey):", n_mods, "\n")
  sz <- sort(table(moduleColors), decreasing = TRUE)
  print(sz)

  # --- 模块特征基因 ---
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs  <- orderMEs(MEs0)

  # --- 模块-性状 Pearson 相关(展示用) ---
  modTraitCor <- cor(MEs, traitData, use = "p")

  # --- ★ 修正核心: 用 ANOVA(3组×3重复) 代替 corPvalueStudent(n=9) ---
  anova_res <- t(sapply(colnames(MEs), function(m) {
    y <- MEs[[m]]
    fit <- aov(y ~ grp)
    a <- summary(fit)[[1]]
    gm <- tapply(y, grp, mean)
    c(F_value = as.numeric(a[["F value"]][1]),
      P_value = as.numeric(a[["Pr(>F)"]][1]),
      mean_Ctrl13 = as.numeric(gm["Ctrl13"]),
      mean_T24 = as.numeric(gm["T24"]),
      mean_T32 = as.numeric(gm["T32"]))
  }))
  anova_df <- data.frame(
    module      = colnames(MEs),
    color       = gsub("^ME", "", colnames(MEs)),
    n_genes     = as.integer(sz[gsub("^ME", "", colnames(MEs))]),
    mean_Ctrl13 = round(anova_res[, "mean_Ctrl13"], 3),
    mean_T24    = round(anova_res[, "mean_T24"], 3),
    mean_T32    = round(anova_res[, "mean_T32"], 3),
    cor_Temperature = round(modTraitCor[, "Temperature"], 3),
    ANOVA_F     = round(anova_res[, "F_value"], 2),
    ANOVA_p     = signif(anova_res[, "P_value"], 3),
    stringsAsFactors = FALSE
  )
  anova_df$ANOVA_FDR <- signif(p.adjust(anova_df$ANOVA_p, "BH"), 3)
  anova_df <- anova_df[order(anova_df$ANOVA_p), ]

  # --- Hub 基因(kME top 30), 对显著模块(FDR<0.05)且 size>=100 ---
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  sig_colors <- anova_df$color[anova_df$ANOVA_FDR < 0.05]
  sig_colors <- sig_colors[sig_colors != "grey"]
  mod_sizes  <- table(moduleColors)
  target_mods <- sig_colors[sig_colors %in% names(mod_sizes)[mod_sizes >= 100]]

  cat("\n显著模块(FDR<0.05 且 size>=100):", paste(target_mods, collapse = ", "), "\n")

  hub_list <- list()
  for (mod in target_mods) {
    modGenes <- colnames(datExpr)[moduleColors == mod]
    kME <- abs(geneModuleMembership[modGenes, paste0("ME", mod)])
    hubGenes <- modGenes[order(kME, decreasing = TRUE)[1:min(30, length(modGenes))]]
    hub_df <- data.frame(gene = hubGenes, module = mod,
                         kME = sort(kME, decreasing = TRUE)[1:length(hubGenes)],
                         stringsAsFactors = FALSE)
    hub_list[[mod]] <- hub_df
    write.csv(hub_df, file.path(OUT, tag, paste0("hub_", mod, ".csv")), row.names = FALSE)
    cat("  模块", mod, ":", length(modGenes), "基因, top", nrow(hub_df), "hub 已保存\n")
  }

  # --- 保存模块基因表 ---
  write.csv(data.frame(gene = colnames(datExpr), module = moduleColors),
            file.path(OUT, tag, "module_genes.csv"), row.names = FALSE)
  write.csv(anova_df, file.path(OUT, tag, "module_trait_ANOVA.csv"), row.names = FALSE)

  # --- 软阈值诊断 ---
  sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, 12, 14, 16, 18, 20), verbose = 0)
  write.csv(data.frame(power = sft$fitIndices[, 1],
                       signedR2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
                       mean_k   = sft$fitIndices[, 5]),
            file.path(OUT, tag, "soft_threshold.csv"), row.names = FALSE)

  list(net = net, MEs = MEs, moduleColors = moduleColors,
       anova_df = anova_df, hub_list = hub_list, sizes = sz,
       sft = sft, target_mods = target_mods)
}

# ---------- 5. 跑两个 power ----------
dir.create(file.path(OUT, "power7"),  showWarnings = FALSE)
dir.create(file.path(OUT, "power18"), showWarnings = FALSE)

res7  <- run_wgcna(power_val = 7,  minModSize = 60, mergeHeight = 0.35, tag = "power7")
saveRDS(res7$net,  file.path(OUT, "power7", "network.rds"))
saveRDS(res7$MEs,  file.path(OUT, "power7", "module_eigengenes.rds"))

res18 <- run_wgcna(power_val = 18, minModSize = 60, mergeHeight = 0.35, tag = "power18")
saveRDS(res18$net, file.path(OUT, "power18", "network.rds"))
saveRDS(res18$MEs, file.path(OUT, "power18", "module_eigengenes.rds"))

# ---------- 6. 对比表 ----------
cat("\n=== 生成对比表 ===\n")

# 6a. 各 power 下的模块概况
cmp <- rbind(
  data.frame(power = 7,  n_modules = length(unique(res7$moduleColors)),
             n_genes_used = ncol(datExpr),
             grey_size = as.integer(res7$sizes["grey"]),
             n_small_modules = sum(res7$sizes < 130),
             n_sig_modules_FDR05 = sum(res7$anova_df$ANOVA_FDR < 0.05)),
  data.frame(power = 18, n_modules = length(unique(res18$moduleColors)),
             n_genes_used = ncol(datExpr),
             grey_size = as.integer(res18$sizes["grey"]),
             n_small_modules = sum(res18$sizes < 130),
             n_sig_modules_FDR05 = sum(res18$anova_df$ANOVA_FDR < 0.05))
)
write.csv(cmp, file.path(OUT, "WGCNA_power_comparison.csv"), row.names = FALSE)
cat("\n模块概况对比:\n"); print(cmp)

# 6b. 关键基因(hub) 在两 power 下的归属
key_genes_mcor <- c("DNAJB5", "DNAJA4", "HSPBP1", "HSPA5", "BAG3", "NFKBIA", "BCL2", "BAX")

# 用 MCOR ID 定位(pending: symbol 映射), 这里先按 symbol 匹配——需要注释表
anno_path <- file.path(BASE, "results/RNAseq_analysis/Results/gene_annotation_table.csv")
has_anno <- file.exists(anno_path)
if (has_anno) {
  anno <- read.csv(anno_path, stringsAsFactors = FALSE)
  cat("\n注释表列名:", paste(colnames(anno), collapse = ", "), "\n")

  # gene_id ↔ symbol 映射
  id2sym <- setNames(anno$Preferred_name, anno$gene_id)

  # 每个 hub 表转成 symbol
  hub_sym <- function(hub_list) {
    out <- list()
    for (mod in names(hub_list)) {
      ids <- hub_list[[mod]]$gene
      out[[mod]] <- unique(na.omit(id2sym[ids]))
    }
    out
  }
  hs7  <- hub_sym(res7$hub_list)
  hs18 <- hub_sym(res18$hub_list)

  key_cmp <- data.frame(gene = key_genes_mcor)
  key_cmp$power7_module  <- sapply(key_genes_mcor, function(g) {
    hit <- names(hs7)[sapply(hs7, function(x) g %in% x)]
    if (length(hit) == 0) "非hub" else paste(hit, collapse = ",")
  })
  key_cmp$power18_module <- sapply(key_genes_mcor, function(g) {
    hit <- names(hs18)[sapply(hs18, function(x) g %in% x)]
    if (length(hit) == 0) "非hub" else paste(hit, collapse = ",")
  })
  write.csv(key_cmp, file.path(OUT, "key_genes_hub_membership.csv"), row.names = FALSE)
  cat("\n关键基因 hub 归属对比:\n"); print(key_cmp)
} else {
  cat("\n[警告] 未找到 gene_annotation_table.csv, 跳过 symbol 对比\n")
}

# 6c. 模块大小对比
size_cmp <- merge(
  data.frame(color = names(res7$sizes),  size_power7  = as.integer(res7$sizes)),
  data.frame(color = names(res18$sizes), size_power18 = as.integer(res18$sizes)),
  by = "color", all = TRUE)
size_cmp[is.na(size_cmp)] <- 0
write.csv(size_cmp, file.path(OUT, "module_size_comparison.csv"), row.names = FALSE)

cat("\n=== 全部完成 ===\n")
cat("输出目录:", OUT, "\n")
cat("  power7/  : ", list.files(file.path(OUT, "power7")), "\n")
cat("  power18/ : ", list.files(file.path(OUT, "power18")), "\n")
cat("  对比表   : WGCNA_power_comparison.csv, key_genes_hub_membership.csv, module_size_comparison.csv\n")
