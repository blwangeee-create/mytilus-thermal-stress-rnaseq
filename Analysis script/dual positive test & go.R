# ============================================================
# Consensus DEG 的 GO & KEGG ORA + GSEA 交叉验证
# 物种: Mytilus coruscus (厚壳贻贝) - 非模式生物，用 eggNOG 自建 term2gene
# ============================================================

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(enrichplot)

# 设置工作目录（按你的实际路径）
setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# -----------------------------------------------------------
# 1. 读取背景基因（count matrix 中检测到的所有基因）
# -----------------------------------------------------------
counts <- read.table("results/counts/counts_s0.txt", header = TRUE, row.names = 1, sep = "\t")
bg_genes <- rownames(counts)
cat("背景基因总数:", length(bg_genes), "\n")

# -----------------------------------------------------------
# 2. 读取 eggNOG 自建的 term2gene 映射
# -----------------------------------------------------------
go_t2g <- readRDS("eggnog_annotation/go_term2gene.rds")
kegg_t2g <- readRDS("eggnog_annotation/kegg_term2gene.rds")
cat("GO term2gene 行数:", nrow(go_t2g), "\n")
cat("KEGG term2gene 行数:", nrow(kegg_t2g), "\n")

# -----------------------------------------------------------
# 3. 定义一个函数：对指定比较组跑 GO + KEGG ORA
# -----------------------------------------------------------
run_ora <- function(deg_file, comparison_name) {
  
  cat("\n========================================\n")
  cat("Processing:", comparison_name, "\n")
  cat("========================================\n")
  
  # 读取 Consensus DEG（假设每行一个基因名，或第一列是基因名）
  deg_df <- read.table(deg_file, header = TRUE, stringsAsFactors = FALSE)
  
  # 自动识别基因名列（通常是第一列）
  if (ncol(deg_df) == 1) {
    gene_list <- deg_df[[1]]
  } else {
    gene_list <- deg_df[[1]]  # 假设第一列是基因ID
  }
  gene_list <- unique(as.character(gene_list))
  cat("Consensus DEG 数量:", length(gene_list), "\n")
  
  # --- GO ORA ---
  go_enrich <- enricher(
    gene          = gene_list,
    universe      = bg_genes,
    TERM2GENE     = go_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  # --- KEGG ORA ---
  kegg_enrich <- enricher(
    gene          = gene_list,
    universe      = bg_genes,
    TERM2GENE     = kegg_t2g,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  # --- 导出结果 ---
  if (!is.null(go_enrich) && nrow(as.data.frame(go_enrich)) > 0) {
    write.csv(as.data.frame(go_enrich), 
              paste0("results/RNAseq_analysis/Results/GO_ORA_", comparison_name, ".csv"),
              row.names = FALSE)
    cat("GO 显著通路数:", nrow(as.data.frame(go_enrich)), "\n")
    
    # 画图
    p <- dotplot(go_enrich, showCategory = 20, title = paste("GO ORA", comparison_name)) + 
      theme_bw(base_size = 12)
    ggsave(paste0("results/RNAseq_analysis/Results/GO_dotplot_", comparison_name, ".pdf"),
           p, width = 10, height = 8)
  } else {
    cat("GO: 无显著富集通路\n")
  }
  
  if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
    write.csv(as.data.frame(kegg_enrich), 
              paste0("results/RNAseq_analysis/Results/KEGG_ORA_", comparison_name, ".csv"),
              row.names = FALSE)
    cat("KEGG 显著通路数:", nrow(as.data.frame(kegg_enrich)), "\n")
    
    p <- dotplot(kegg_enrich, showCategory = 20, title = paste("KEGG ORA", comparison_name)) + 
      theme_bw(base_size = 12)
    ggsave(paste0("results/RNAseq_analysis/Results/KEGG_dotplot_", comparison_name, ".pdf"),
           p, width = 10, height = 8)
  } else {
    cat("KEGG: 无显著富集通路\n")
  }
  
  # 返回结果供后续交叉验证使用
  return(list(go = go_enrich, kegg = kegg_enrich, genes = gene_list))
}

# -----------------------------------------------------------
# 4. 对三个比较组分别跑 ORA
# -----------------------------------------------------------
ora_t24 <- run_ora("results/RNAseq_analysis/Results/HighConf_DEG_T24vsCtrl.txt", "T24vsCtrl")
ora_t32 <- run_ora("results/RNAseq_analysis/Results/HighConf_DEG_T32vsCtrl.txt", "T32vsCtrl")
ora_t32v24 <- run_ora("results/RNAseq_analysis/Results/HighConf_DEG_T32vsT24.txt", "T32vsT24")

cat("\n========================================\n")
cat("ORA 分析完成，结果保存在 results/RNAseq_analysis/Results/ 目录\n")
cat("========================================\n")

# -----------------------------------------------------------
# 5. 交叉验证：提取 GSEA 结果（从 RData 中加载）
# -----------------------------------------------------------
# 你的 GSEA 结果很可能在 RNAseq_final_analysis.RData 或 Phase1_DEG_analysis.RData 里
# 先加载查看对象名称

cat("\n========================================\n")
cat("Step 5: 从 RData 中提取 GSEA 结果\n")
cat("========================================\n")

# 尝试加载（先不赋值，避免覆盖当前环境）
load("results/RNAseq_analysis/Results/RNAseq_final_analysis.RData", verbose = TRUE)

# 查看环境中所有对象，找 gseaResult 类的对象
cat("\n当前环境中的对象:\n")
print(ls())

# 自动查找 gseaResult 对象
gsea_objects <- ls()[sapply(ls(), function(x) {
  obj <- get(x)
  inherits(obj, "gseaResult")
})]

cat("\n检测到的 GSEA 结果对象:\n")
print(gsea_objects)

# -----------------------------------------------------------
# 6. 交叉验证（需要你确认 GSEA 对象名后运行）
# -----------------------------------------------------------
# 假设你的 GSEA 对象叫 gsea_go_t32 和 gsea_kegg_t32（请根据实际情况修改）

# 示例代码（取消注释并修改对象名后使用）：
#
# gsea_go_df <- as.data.frame(gsea_go_t32)      # 替换为实际的 GSEA GO 对象名
# gsea_kegg_df <- as.data.frame(gsea_kegg_t32)  # 替换为实际的 GSEA KEGG 对象名
#
# ora_go_df <- as.data.frame(ora_t32$go)
# ora_kegg_df <- as.data.frame(ora_t32$kegg)
#
# # 找双阳性通路
# overlap_go <- intersect(ora_go_df$ID, gsea_go_df$ID[gsea_go_df$p.adjust < 0.05])
# overlap_kegg <- intersect(ora_kegg_df$ID, gsea_kegg_df$ID[gsea_kegg_df$p.adjust < 0.05])
#
# cat("GO 双阳性通路:", length(overlap_go), "\n")
# cat("KEGG 双阳性通路:", length(overlap_kegg), "\n")
#
# write.csv(data.frame(Pathway_ID = overlap_go), 
#           "results/RNAseq_analysis/Results/Dual_Positive_GO_T32vsCtrl.csv", 
#           row.names = FALSE)
# write.csv(data.frame(Pathway_ID = overlap_kegg), 
#           "results/RNAseq_analysis/Results/Dual_Positive_KEGG_T32vsCtrl.csv", 
#           row.names = FALSE)
# ============================================================
# 补齐缺失的 GSEA 对象（T24 KEGG + T32vsT24 GO/KEGG）
# 基于你之前的 eggNOG 注释流程，精简版
# ============================================================

library(clusterProfiler)
library(ggplot2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# -----------------------------------------------------------
# 1. 读取注释（只读一次）
# -----------------------------------------------------------
cat("========== 1. 读取 eggNOG 注释 ==========\n")

eggnog <- read.delim(
  "annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  quote = "", comment.char = "", check.names = FALSE
)

# ID 映射
ft <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
id_map <- ft[, c("product_accession", "locus_tag")]
id_map <- id_map[id_map$product_accession != "" & !is.na(id_map$product_accession) &
                   id_map$locus_tag != "" & !is.na(id_map$locus_tag), ]
id_map <- id_map[!duplicated(id_map$product_accession), ]
id_map$gene_id <- paste0("gene-", id_map$locus_tag)

m_idx <- match(eggnog[, 1], id_map[["product_accession"]])
eggnog[["gene_id"]] <- id_map[["gene_id"]][m_idx]
eggnog <- eggnog[!is.na(eggnog$gene_id), ]
cat("成功映射基因:", nrow(eggnog), "\n")

# -----------------------------------------------------------
# 2. 构建 TERM2GENE（只构建一次）
# -----------------------------------------------------------
cat("\n========== 2. 构建 TERM2GENE ==========\n")

# KEGG
kegg_list <- strsplit(eggnog$KEGG_Pathway, ",")
kegg_df <- data.frame(
  term = unlist(kegg_list),
  gene = rep(eggnog$gene_id, sapply(kegg_list, length)),
  stringsAsFactors = FALSE
)
kegg_df <- kegg_df[kegg_df$term != "" & kegg_df$term != "-" & !is.na(kegg_df$term), ]
kegg_df <- unique(kegg_df)
kegg_df <- kegg_df[grepl("^map", kegg_df$term), ]
cat("KEGG term2gene:", nrow(kegg_df), "对\n")

# GO
go_list <- strsplit(eggnog$GOs, ",")
go_df <- data.frame(
  term = unlist(go_list),
  gene = rep(eggnog$gene_id, sapply(go_list, length)),
  stringsAsFactors = FALSE
)
go_df <- go_df[go_df$term != "" & go_df$term != "-" & !is.na(go_df$term), ]
go_df <- unique(go_df)
cat("GO term2gene:", nrow(go_df), "对\n")

# -----------------------------------------------------------
# 3. 通用 GSEA 函数
# -----------------------------------------------------------
run_gsea <- function(comparison, res_csv, pvalue_cutoff = 0.05) {
  
  cat("\n========================================\n")
  cat("Running GSEA:", comparison, "\n")
  cat("========================================\n")
  
  # 读取 DESeq2 结果
  res <- read.csv(file.path("results/RNAseq_analysis/Results", res_csv),
                  row.names = 1, stringsAsFactors = FALSE)
  res_clean <- res[!is.na(res$stat), ]
  geneList <- res_clean$stat
  names(geneList) <- rownames(res_clean)
  geneList <- sort(geneList, decreasing = TRUE)
  cat("输入基因:", length(geneList), "\n")
  
  # KEGG GSEA
  gse_kegg <- NULL
  n_kegg <- sum(names(geneList) %in% kegg_df$gene)
  cat("有 KEGG 注释:", n_kegg, "\n")
  if (n_kegg > 0) {
    gse_kegg <- GSEA(
      geneList = geneList, TERM2GENE = kegg_df,
      pvalueCutoff = pvalue_cutoff, pAdjustMethod = "BH",
      minGSSize = 10, maxGSSize = 500, verbose = FALSE, seed = 123
    )
    cat("KEGG 显著通路:", nrow(gse_kegg), "\n")
  }
  
  # GO GSEA
  gse_go <- NULL
  n_go <- sum(names(geneList) %in% go_df$gene)
  cat("有 GO 注释:", n_go, "\n")
  if (n_go > 0) {
    gse_go <- GSEA(
      geneList = geneList, TERM2GENE = go_df,
      pvalueCutoff = pvalue_cutoff, pAdjustMethod = "BH",
      minGSSize = 10, maxGSSize = 500, verbose = FALSE, seed = 123
    )
    cat("GO 显著 term:", nrow(gse_go), "\n")
  }
  
  # 保存对象（关键！）
  if (!is.null(gse_kegg)) {
    obj_name <- paste0("kegg_gsea_", gsub("_vs_", "v", comparison))
    assign(obj_name, gse_kegg, envir = .GlobalEnv)
    save(list = obj_name, 
         file = paste0("results/RNAseq_analysis/Results/", obj_name, ".RData"))
    cat(">>> 已保存:", obj_name, "\n")
  }
  
  if (!is.null(gse_go)) {
    obj_name <- paste0("go_gsea_", gsub("_vs_", "v", comparison))
    assign(obj_name, gse_go, envir = .GlobalEnv)
    save(list = obj_name,
         file = paste0("results/RNAseq_analysis/Results/", obj_name, ".RData"))
    cat(">>> 已保存:", obj_name, "\n")
  }
  
  return(list(kegg = gse_kegg, go = gse_go))
}

# -----------------------------------------------------------
# 4. 跑三个缺失的 GSEA
# -----------------------------------------------------------

# 4.1 T24 KEGG GSEA（GO 已有，只补 KEGG）
gsea_t24_kegg <- run_gsea("T24_vs_Ctrl", "DESeq2_T24_vs_Ctrl.csv")

# 4.2 T32vsT24 GO + KEGG GSEA
gsea_t32v24 <- run_gsea("T32_vs_T24", "DESeq2_T32_vs_T24.csv")

cat("\n========================================\n")
cat("全部 GSEA 补齐完成，对象已保存\n")
cat("========================================\n")




# ============================================================
# 完整交叉验证：ORA vs GSEA（三个比较组）
# ============================================================

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 加载刚保存的 T32vsCtrl 对象（如果刚跑完还在内存里，这步可跳过）
if (!exists("go_gsea_T32vCtrl")) {
  load("results/RNAseq_analysis/Results/go_gsea_T32vCtrl.RData")
  load("results/RNAseq_analysis/Results/kegg_gsea_T32vCtrl.RData")
}

# 加载 T24 和 T32vsT24 对象（如果刚跑完还在内存里，可跳过）
if (!exists("kegg_gsea_T24vCtrl")) {
  load("results/RNAseq_analysis/Results/kegg_gsea_T24vCtrl.RData")
  load("results/RNAseq_analysis/Results/go_gsea_T32vT24.RData")
  load("results/RNAseq_analysis/Results/kegg_gsea_T32vT24.RData")
}

# 确认所有对象就位
cat("=== GSEA 对象检查 ===\n")
all_gsea <- c("go_gsea_T24vCtrl", "kegg_gsea_T24vCtrl",
              "go_gsea_T32vCtrl", "kegg_gsea_T32vCtrl",
              "go_gsea_T32vT24", "kegg_gsea_T32vT24")
for (obj in all_gsea) {
  cat(obj, "存在:", exists(obj), "| 行数:", 
      ifelse(exists(obj), nrow(as.data.frame(get(obj))), "NA"), "\n")
}

# -----------------------------------------------------------
# 交叉验证函数
# -----------------------------------------------------------
cross_validate <- function(ora_obj, gsea_go_obj, gsea_kegg_obj, comparison_name) {
  
  cat("\n========================================\n")
  cat("交叉验证:", comparison_name, "\n")
  cat("========================================\n")
  
  results <- list()
  
  # --- GO ---
  if (!is.null(ora_obj$go) && !is.null(gsea_go_obj) && nrow(as.data.frame(ora_obj$go)) > 0) {
    ora_df <- as.data.frame(ora_obj$go)
    gsea_df <- as.data.frame(gsea_go_obj)
    
    ora_sig <- ora_df$ID[ora_df$p.adjust < 0.05]
    gsea_sig <- gsea_df$ID[gsea_df$p.adjust < 0.05]
    dual <- intersect(ora_sig, gsea_sig)
    
    cat("GO ORA 显著:", length(ora_sig), "| GSEA 显著:", length(gsea_sig), "| 双阳性:", length(dual), "\n")
    
    if (length(dual) > 0) {
      dual_df <- ora_df[ora_df$ID %in% dual, c("ID", "Description", "p.adjust", "Count")]
      gsea_match <- gsea_df[gsea_df$ID %in% dual, c("ID", "NES", "p.adjust")]
      names(gsea_match)[2:3] <- c("GSEA_NES", "GSEA_p.adjust")
      dual_df <- merge(dual_df, gsea_match, by = "ID")
      write.csv(dual_df, paste0("results/RNAseq_analysis/Results/Dual_Positive_GO_", comparison_name, ".csv"), row.names = FALSE)
      cat("  -> 已保存 Dual_Positive_GO_", comparison_name, ".csv\n")
    }
    results$dual_go <- dual
    results$ora_only_go <- setdiff(ora_sig, gsea_sig)
    results$gsea_only_go <- setdiff(gsea_sig, ora_sig)
  }
  
  # --- KEGG ---
  if (!is.null(ora_obj$kegg) && !is.null(gsea_kegg_obj) && nrow(as.data.frame(ora_obj$kegg)) > 0) {
    ora_df <- as.data.frame(ora_obj$kegg)
    gsea_df <- as.data.frame(gsea_kegg_obj)
    
    ora_sig <- ora_df$ID[ora_df$p.adjust < 0.05]
    gsea_sig <- gsea_df$ID[gsea_df$p.adjust < 0.05]
    dual <- intersect(ora_sig, gsea_sig)
    
    cat("KEGG ORA 显著:", length(ora_sig), "| GSEA 显著:", length(gsea_sig), "| 双阳性:", length(dual), "\n")
    
    if (length(dual) > 0) {
      dual_df <- ora_df[ora_df$ID %in% dual, c("ID", "Description", "p.adjust", "Count")]
      gsea_match <- gsea_df[gsea_df$ID %in% dual, c("ID", "NES", "p.adjust")]
      names(gsea_match)[2:3] <- c("GSEA_NES", "GSEA_p.adjust")
      dual_df <- merge(dual_df, gsea_match, by = "ID")
      write.csv(dual_df, paste0("results/RNAseq_analysis/Results/Dual_Positive_KEGG_", comparison_name, ".csv"), row.names = FALSE)
      cat("  -> 已保存 Dual_Positive_KEGG_", comparison_name, ".csv\n")
    }
    results$dual_kegg <- dual
    results$ora_only_kegg <- setdiff(ora_sig, gsea_sig)
    results$gsea_only_kegg <- setdiff(gsea_sig, ora_sig)
  }
  
  return(results)
}

# -----------------------------------------------------------
# 运行三个比较组的交叉验证
# -----------------------------------------------------------

cv_t24 <- cross_validate(ora_t24, go_gsea_T24vCtrl, kegg_gsea_T24vCtrl, "T24vsCtrl")
cv_t32 <- cross_validate(ora_t32, go_gsea_T32vCtrl, kegg_gsea_T32vCtrl, "T32vsCtrl")
cv_t32v24 <- cross_validate(ora_t32v24, go_gsea_T32vT24, kegg_gsea_T32vT24, "T32vsT24")

# -----------------------------------------------------------
# 汇总表
# -----------------------------------------------------------
cat("\n========================================\n")
cat("交叉验证汇总\n")
cat("========================================\n")

summary_all <- data.frame(
  Comparison = c("T24vsCtrl", "T32vsCtrl", "T32vsT24"),
  GO_Dual = c(length(cv_t24$dual_go), length(cv_t32$dual_go), length(cv_t32v24$dual_go)),
  GO_ORA_Only = c(length(cv_t24$ora_only_go), length(cv_t32$ora_only_go), length(cv_t32v24$ora_only_go)),
  GO_GSEA_Only = c(length(cv_t24$gsea_only_go), length(cv_t32$gsea_only_go), length(cv_t32v24$gsea_only_go)),
  KEGG_Dual = c(length(cv_t24$dual_kegg), length(cv_t32$dual_kegg), length(cv_t32v24$dual_kegg)),
  KEGG_ORA_Only = c(length(cv_t24$ora_only_kegg), length(cv_t32$ora_only_kegg), length(cv_t32v24$ora_only_kegg)),
  KEGG_GSEA_Only = c(length(cv_t24$gsea_only_kegg), length(cv_t32$gsea_only_kegg), length(cv_t32v24$gsea_only_kegg))
)

print(summary_all)
write.csv(summary_all, "results/RNAseq_analysis/Results/CrossValidation_Summary_All.csv", row.names = FALSE)

cat("\n========================================\n")
cat("全部交叉验证完成！\n")
cat("========================================\n")
cat("\n核心发现（T32vsCtrl - 你的主效应）：\n")
cat("  GO 双阳性通路:", length(cv_t32$dual_go), "\n")
cat("  KEGG 双阳性通路:", length(cv_t32$dual_kegg), "\n")
cat("  GO 仅 ORA:", length(cv_t32$ora_only_go), "| 仅 GSEA:", length(cv_t32$gsea_only_go), "\n")
cat("  KEGG 仅 ORA:", length(cv_t32$ora_only_kegg), "| 仅 GSEA:", length(cv_t32$gsea_only_kegg), "\n")
setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取 T32vsCtrl 双阳性 GO
dual_go <- read.csv("results/RNAseq_analysis/Results/Dual_Positive_GO_T32vsCtrl.csv")

# 筛选核心通路
core_go <- dual_go %>%
  filter(GSEA_NES > 0, Count >= 5) %>%
  arrange(p.adjust) %>%
  head(20)  # 取 Top 20

write.csv(core_go, "results/RNAseq_analysis/Results/Core_GO_Dual_Positive_T32vsCtrl.csv", row.names = FALSE)

cat("核心 GO 双阳性通路:", nrow(core_go), "条\n")
print(core_go[, c("Description", "Count", "p.adjust", "GSEA_NES", "GSEA_p.adjust")])
setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取 T24vsCtrl 双阳性
dual_go_t24 <- read.csv("results/RNAseq_analysis/Results/Dual_Positive_GO_T24vsCtrl.csv")

core_go_t24 <- dual_go_t24 %>%
  filter(GSEA_NES > 0, Count >= 5) %>%
  arrange(p.adjust) %>%
  head(20)

write.csv(core_go_t24, "results/RNAseq_analysis/Results/Core_GO_Dual_Positive_T24vsCtrl.csv", row.names = FALSE)

cat("T24 核心 GO 双阳性:", nrow(core_go_t24), "条\n")
print(core_go_t24[, c("Description", "Count", "p.adjust", "GSEA_NES")])
setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# T32vsT24 双阳性核心（虽然只有 12 条，但补充材料需要）
dual_go_t32v24 <- read.csv("results/RNAseq_analysis/Results/Dual_Positive_GO_T32vsT24.csv")

core_go_t32v24 <- dual_go_t32v24 %>%
  filter(GSEA_NES > 0, Count >= 3) %>%  # Count 降到 3，因为总共只有 8 个 DEG
  arrange(p.adjust)

write.csv(core_go_t32v24, "results/RNAseq_analysis/Results/Core_GO_Dual_Positive_T32vsT24.csv", row.names = FALSE)

cat("T32vsT24 双阳性:", nrow(dual_go_t32v24), "条 | 核心 (NES>0):", nrow(core_go_t32v24), "条\n")
print(core_go_t32v24[, c("Description", "Count", "p.adjust", "GSEA_NES")])

install.packages("ggrepel")
# ============================================================
# Fig 3C-D: ORA vs GSEA 散点图（T32 + T24）
# ============================================================
# 第一步：删掉内存里的旧函数（关键！）
rm(make_scatter_premium)

# 第二步：重新定义完整函数（复制下面整块）
library(ggplot2)
library(dplyr)
library(ggrepel)

make_scatter_premium <- function(ora_file, gsea_obj, title_text) {
  
  ora_df <- read.csv(ora_file)
  gsea_df <- as.data.frame(gsea_obj)
  
  merged <- merge(
    ora_df[, c("ID", "Description", "p.adjust", "Count")],
    gsea_df[, c("ID", "NES", "p.adjust")],
    by = "ID", all = TRUE,
    suffixes = c("_ORA", "_GSEA")
  )
  
  merged$p.adjust_ORA[is.na(merged$p.adjust_ORA)] <- 1
  merged$x_val <- -log10(merged$p.adjust_ORA)
  merged$NES[is.na(merged$NES)] <- 0
  merged$p.adjust_GSEA[is.na(merged$p.adjust_GSEA)] <- 1
  merged$Count[is.na(merged$Count)] <- 0
  
  ora_sig  <- merged$p.adjust_ORA < 0.05
  gsea_sig <- merged$p.adjust_GSEA < 0.05
  
  merged$category <- "Neither"
  merged$category[ora_sig & !gsea_sig] <- "ORA only"
  merged$category[!ora_sig & gsea_sig] <- "GSEA only"
  merged$category[ora_sig & gsea_sig]  <- "Dual-positive"
  merged$category <- factor(merged$category, 
                            levels = c("Dual-positive", "ORA only", "GSEA only", "Neither"))
  
  tbl <- table(merged$category)
  cat("\n===", title_text, "===\n")
  print(tbl)
  
  bg  <- merged[merged$category %in% c("Neither", "GSEA only"), ]
  fg  <- merged[merged$category %in% c("Dual-positive", "ORA only"), ]
  
  set.seed(42)
  gsea_only <- bg[bg$category == "GSEA only", ]
  if (nrow(gsea_only) > 150) gsea_only <- gsea_only[sample(nrow(gsea_only), 150), ]
  neither <- bg[bg$category == "Neither", ]
  if (nrow(neither) > 100) neither <- neither[sample(nrow(neither), 100), ]
  bg_plot <- rbind(gsea_only, neither)
  
  p <- ggplot() +
    geom_point(data = bg_plot, 
               aes(x = x_val, y = NES, color = category),
               size = 2.8, alpha = 0.38, stroke = 0,
               position = position_jitter(width = 0.1, height = 0, seed = 123)) +
    geom_point(data = fg, 
               aes(x = x_val, y = NES, color = category, size = Count),
               alpha = 0.55, stroke = 0.1) +
    scale_color_manual(values = c("Dual-positive" = "#C74375", "ORA only" = "#3070B3", 
                                  "GSEA only" = "#FFCD00", "Neither" = "#6A757E"), 
                       drop = FALSE) +
    scale_size_continuous(range = c(3.5, 9), 
                          breaks = c(5, 15, 30),
                          labels = c("5", "15", "30"),
                          name = "Gene Count") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    scale_x_continuous(limits = c(-0.5, max(merged$x_val, na.rm = TRUE) * 1.05), 
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(min(merged$NES, na.rm = TRUE) * 1.15, 
                                  max(merged$NES, na.rm = TRUE) * 1.15)) +
    labs(x = expression(-log[10]~(ORA~FDR)), y = "GSEA NES",
         color = "Category", title = title_text) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.6, "cm"))
  
  # ============================================
  # 标注部分（新版 - 用 Name 而不是 GO ID）
  # ============================================
  comp_name <- gsub(" vs ", "vs", title_text)
  dual_file <- paste0("results/RNAseq_analysis/Results/Core_GO_Dual_Positive_", comp_name, ".csv")
  
  if (file.exists(dual_file)) {
    dual_named <- read.csv(dual_file, stringsAsFactors = FALSE)
    
    # 过滤 obsolete，按 NES 降序取 Top 4
    dual_top <- dual_named %>%
      filter(GSEA_NES > 0, !grepl("^obsolete", Name, ignore.case = TRUE)) %>%
      arrange(desc(GSEA_NES)) %>%
      head(4)
    
    if (nrow(dual_top) > 0) {
      p <- p + geom_text_repel(
        data = dual_top,
        aes(x = -log10(p.adjust), y = GSEA_NES, label = Name),
        size = 3.3, max.overlaps = 8, color = "black",
        box.padding = 0.4, point.padding = 0.3,
        segment.size = 0.3, fontface = "italic"
      )
    }
  }
  
  # 统计标注
  p <- p + annotate(
    "text", x = Inf, y = -Inf,
    label = paste(names(tbl), tbl, sep = ": ", collapse = "\n"),
    hjust = 1, vjust = 0, size = 3.5, color = "grey40", fontface = "italic"
  )
  
  return(p)
}

# 第三步：重新生成图
p_t32_prem <- make_scatter_premium(
  "results/RNAseq_analysis/Results/GO_ORA_T32vsCtrl.csv",
  go_gsea_T32vCtrl,
  "T32 vs Ctrl"
)

print(p_t32_prem)
# T24 散点图
p_t24_prem <- make_scatter_premium(
  "results/RNAseq_analysis/Results/GO_ORA_T24vsCtrl.csv",
  go_gsea_T24vCtrl,
  "T24 vs Ctrl"
)

print(p_t24_prem)

library(ggplot2)
library(dplyr)
library(GO.db)
library(stringr)
library(patchwork)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 辅助函数 ==========
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(y) as.numeric(y[1])/as.numeric(y[2]))
}

wrap_name <- function(x, width = 45) {
  sapply(x, function(s) str_wrap(s, width = width))
}

get_go_info <- function(go_ids) {
  go_ids <- unique(as.character(go_ids))
  go_info <- AnnotationDbi::select(GO.db, keys = go_ids, 
                                   columns = c("TERM", "ONTOLOGY"), 
                                   keytype = "GOID")
  name_map <- setNames(go_info$TERM, go_info$GOID)
  ont_map  <- setNames(go_info$ONTOLOGY, go_info$GOID)
  list(names = name_map, ontology = ont_map)
}

# ========== T32vsCtrl (Fig 3A) ==========
go_t32 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T32vsCtrl.csv", stringsAsFactors = FALSE)
go_t32 <- go_t32 %>% filter(!grepl("^obsolete", Description, ignore.case = TRUE))

go_info <- get_go_info(go_t32$ID)
go_t32$DisplayName <- go_info$names[go_t32$ID]
go_t32$DisplayName[is.na(go_t32$DisplayName)] <- go_t32$Description[is.na(go_t32$DisplayName)]
go_t32$ONTOLOGY <- go_info$ontology[go_t32$ID]

# ===== 插入开始：过滤脊椎动物特异性 term =====
exclude_terms <- c(
  "negative regulation of T cell apoptotic process",
  "Z disc", "I band", "striated muscle dense body",
  "melanosome", "male germ cell nucleus",
  "tumor necrosis factor receptor activity",
  "death receptor activity"
)

go_t32 <- go_t32 %>% 
  filter(!DisplayName %in% exclude_terms) %>%
  filter(!grepl("T cell|striated muscle|Z disc|I band|melanosome|male germ cell", DisplayName, ignore.case = TRUE))
# ===== 插入结束 =====

# 关键：每组取 Top 5（BP 5 + CC 5 + MF 5 = 15）
go_t32_top <- go_t32 %>%
  filter(p.adjust < 0.05, !is.na(ONTOLOGY)) %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust) %>%
  mutate(
    DisplayName = wrap_name(DisplayName, width = 45),
    DisplayName = factor(DisplayName, levels = rev(DisplayName))
  )

p3a <- ggplot(go_t32_top, aes(x = parse_ratio(GeneRatio), y = DisplayName)) +
  geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.85) +
  scale_color_gradient(low = "#FFE082", high = "#E65100", name = "-log10(FDR)") +
  scale_size_continuous(range = c(3, 10), breaks = c(5, 15, 30), name = "Gene Count") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Gene Ratio", y = NULL, title = "A  T32 vs Ctrl") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#E65100"),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank())

print(p3a)

# ========== T24vsCtrl (Fig 3B) ==========
go_t24 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T24vsCtrl.csv", stringsAsFactors = FALSE)
go_t24 <- go_t24 %>% filter(!grepl("^obsolete", Description, ignore.case = TRUE))

go_info_t24 <- get_go_info(go_t24$ID)
go_t24$DisplayName <- go_info_t24$names[go_t24$ID]
go_t24$DisplayName[is.na(go_t24$DisplayName)] <- go_t24$Description[is.na(go_t24$DisplayName)]
go_t24$ONTOLOGY <- go_info_t24$ontology[go_t24$ID]

# ===== 插入开始：过滤脊椎动物特异性 term =====
go_t24 <- go_t24 %>% 
  filter(!DisplayName %in% exclude_terms) %>%
  filter(!grepl("T cell|striated muscle|Z disc|I band|melanosome|male germ cell", DisplayName, ignore.case = TRUE))
# ===== 插入结束 =====

# T24 效应弱，每组取 Top 3（BP 3 + CC 3 + MF 3 = 9）
go_t24_top <- go_t24 %>%
  filter(p.adjust < 0.05, !is.na(ONTOLOGY)) %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust) %>%
  mutate(
    DisplayName = wrap_name(DisplayName, width = 45),
    DisplayName = factor(DisplayName, levels = rev(DisplayName))
  )

p3b <- ggplot(go_t24_top, aes(x = parse_ratio(GeneRatio), y = DisplayName)) +
  geom_point(aes(size = Count, color = -log10(p.adjust)), alpha = 0.85) +
  scale_color_gradient(low = "#90CAF9", high = "#1565C0", name = "-log10(FDR)") +
  scale_size_continuous(range = c(3, 10), breaks = c(5, 15, 30), name = "Gene Count") +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Gene Ratio", y = NULL, title = "B  T24 vs Ctrl") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, color = "#1565C0"),
        axis.text.y = element_text(size = 9),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "right",
        panel.grid.minor = element_blank())

print(p3b)


devtools::install_github("dxsbiocc/gground")
library(ggplot2)
library(dplyr)
library(GO.db)
library(stringr)
library(patchwork)
library(gground)    # geom_round_col / geom_round_rect
library(ggprism)    # theme_prism

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 辅助函数（完全不变） ==========
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(y) as.numeric(y[1])/as.numeric(y[2]))
}

wrap_name <- function(x, width = 45) {
  sapply(x, function(s) str_wrap(s, width = width))
}

get_go_info <- function(go_ids) {
  go_ids <- unique(as.character(go_ids))
  go_info <- AnnotationDbi::select(GO.db, keys = go_ids, 
                                   columns = c("TERM", "ONTOLOGY"), 
                                   keytype = "GOID")
  name_map <- setNames(go_info$TERM, go_info$GOID)
  ont_map  <- setNames(go_info$ONTOLOGY, go_info$GOID)
  list(names = name_map, ontology = ont_map)
}

# ========== 颜色（你示例的配色） ==========
pal <- c("BP" = '#7bc4e2', "CC" = '#acd372', "MF" = '#fbb05b')

# ========== T32vsCtrl：数据准备（完全不变） ==========
go_t32 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T32vsCtrl.csv", stringsAsFactors = FALSE)
go_t32 <- go_t32 %>% filter(!grepl("^obsolete", Description, ignore.case = TRUE))

go_info <- get_go_info(go_t32$ID)
go_t32$DisplayName <- go_info$names[go_t32$ID]
go_t32$DisplayName[is.na(go_t32$DisplayName)] <- go_t32$Description[is.na(go_t32$DisplayName)]
go_t32$ONTOLOGY <- go_info$ontology[go_t32$ID]

exclude_terms <- c(
  "negative regulation of T cell apoptotic process",
  "Z disc", "I band", "striated muscle dense body",
  "melanosome", "male germ cell nucleus",
  "tumor necrosis factor receptor activity",
  "death receptor activity"
)

go_t32 <- go_t32 %>% 
  filter(!DisplayName %in% exclude_terms) %>%
  filter(!grepl("T cell|striated muscle|Z disc|I band|melanosome|male germ cell", DisplayName, ignore.case = TRUE))

go_t32_top <- go_t32 %>%
  filter(p.adjust < 0.05, !is.na(ONTOLOGY)) %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

# 补 geneID（如果过滤后丢了）
if (!"geneID" %in% names(go_t32_top)) {
  raw_t32 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T32vsCtrl.csv", stringsAsFactors = FALSE)
  go_t32_top$geneID <- raw_t32$geneID[match(go_t32_top$ID, raw_t32$ID)]
}
# 基因ID → 名称映射
anno <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors = FALSE)
id2name <- setNames(anno$Preferred_name, anno$gene_id)
id2name[id2name == "-" | is.na(id2name)] <- names(id2name)[id2name == "-" | is.na(id2name)]
# ========== T32：条形图数据构造 ==========
compress <- 0.6   # ← 压缩系数，0.6=缩短到60%，改这个数调长短

go_t32_top <- go_t32_top %>%
  mutate(
    index = row_number(),
    logp  = -log10(p.adjust) * compress,   # ← 乘系数
    DisplayName = factor(DisplayName, levels = DisplayName),
    # geneID 截断显示前3个
    geneID_show = sapply(strsplit(as.character(geneID), "/"), function(x) {
      names <- id2name[x]
      names[is.na(names)] <- x[is.na(names)]
      if (length(names) > 3) paste(paste(names[1:3], collapse = "/"), "...", sep = "")
      else paste(names, collapse = "/")
    })
  )

width_t32 <- 0.5
xaxis_max_t32 <- 8

rect_t32 <- go_t32_top %>%
  group_by(ONTOLOGY) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    xmin = -3 * width_t32,
    xmax = -2 * width_t32,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

# ========== T32：画图（横向圆角条形） ==========
p3a <- ggplot(go_t32_top, aes(x = logp, y = index, fill = ONTOLOGY)) +
  geom_round_col(aes(y = DisplayName), width = 0.6, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = DisplayName), hjust = 0, size = 3.5) +
  geom_text(aes(x = 0.1, label = geneID_show, colour = ONTOLOGY),
            hjust = 0, vjust = 2.6, size = 2.2, fontface = 'italic', show.legend = FALSE) +
  geom_point(aes(x = -width_t32, size = Count), shape = 21) +
  geom_text(aes(x = -width_t32, label = Count), size = 2.5) +
  scale_size_continuous(name = 'Count', range = c(3, 8)) +
  geom_round_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
                  data = rect_t32, radius = unit(2, 'mm'), inherit.aes = FALSE) +
  geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
            data = rect_t32, inherit.aes = FALSE, angle = 90, size = 3, fontface = "bold") +
  geom_segment(aes(x = 0, y = 0, xend = xaxis_max_t32 * compress, yend = 0),   # ← 乘系数
               linewidth = 1.5, inherit.aes = FALSE) +
  labs(y = NULL, title = "A  T32 vs Ctrl") +
  scale_fill_manual(name = 'Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, ceiling(xaxis_max_t32), 2), expand = expansion(c(0, 0))) +
  theme_prism() +
  theme(axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

print(p3a)

# ========== T24vsCtrl：数据准备（完全不变） ==========
go_t24 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T24vsCtrl.csv", stringsAsFactors = FALSE)
go_t24 <- go_t24 %>% filter(!grepl("^obsolete", Description, ignore.case = TRUE))

go_info_t24 <- get_go_info(go_t24$ID)
go_t24$DisplayName <- go_info_t24$names[go_t24$ID]
go_t24$DisplayName[is.na(go_t24$DisplayName)] <- go_t24$Description[is.na(go_t24$DisplayName)]
go_t24$ONTOLOGY <- go_info_t24$ontology[go_t24$ID]

go_t24 <- go_t24 %>% 
  filter(!DisplayName %in% exclude_terms) %>%
  filter(!grepl("T cell|striated muscle|Z disc|I band|melanosome|male germ cell", DisplayName, ignore.case = TRUE))

go_t24_top <- go_t24 %>%
  filter(p.adjust < 0.05, !is.na(ONTOLOGY)) %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

if (!"geneID" %in% names(go_t24_top)) {
  raw_t24 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T24vsCtrl.csv", stringsAsFactors = FALSE)
  go_t24_top$geneID <- raw_t24$geneID[match(go_t24_top$ID, raw_t24$ID)]
}
# 基因ID → 名称映射
anno <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors = FALSE)
id2name <- setNames(anno$Preferred_name, anno$gene_id)
id2name[id2name == "-" | is.na(id2name)] <- names(id2name)[id2name == "-" | is.na(id2name)]
# ========== T24：条形图数据构造 ==========
go_t24_top <- go_t24_top %>%
  mutate(
    index = row_number(),
    logp  = -log10(p.adjust),
    DisplayName = factor(DisplayName, levels = DisplayName),
    geneID_show = sapply(strsplit(as.character(geneID), "/"), function(x) {
      names <- id2name[x]
      names[is.na(names)] <- x[is.na(names)]
      if (length(names) > 3) paste(paste(names[1:3], collapse = "/"), "...", sep = "")
      else paste(names, collapse = "/")
    })
  )

width_t24 <- 0.5
xaxis_max_t24 <- 8

rect_t24 <- go_t24_top %>%
  group_by(ONTOLOGY) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    xmin = -3 * width_t24,
    xmax = -2 * width_t24,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

# ========== T24：画图 ==========
p3b <- ggplot(go_t24_top, aes(x = logp, y = index, fill = ONTOLOGY)) +
  geom_round_col(aes(y = DisplayName), width = 0.6, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = DisplayName), hjust = 0, size = 3.5) +
  geom_text(aes(x = 0.1, label = geneID_show, colour = ONTOLOGY),
            hjust = 0, vjust = 2.6, size = 2.2, fontface = 'italic', show.legend = FALSE) +
  geom_point(aes(x = -width_t24, size = Count), shape = 21) +
  geom_text(aes(x = -width_t24, label = Count), size = 2.5) +
  scale_size_continuous(name = 'Count', range = c(3, 8)) +
  geom_round_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
                  data = rect_t24, radius = unit(2, 'mm'), inherit.aes = FALSE) +
  geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
            data = rect_t24, inherit.aes = FALSE, angle = 90, size = 3, fontface = "bold") +
  geom_segment(aes(x = 0, y = 0, xend = xaxis_max_t24, yend = 0),
               linewidth = 1.5, inherit.aes = FALSE) +
  labs(y = NULL, title = "B  T24 vs Ctrl") +
  scale_fill_manual(name = 'Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(
    breaks = seq(0, ceiling(xaxis_max_t32 * compress), 2),
    labels = seq(0, ceiling(xaxis_max_t32 * compress), 2) / compress,   # ← 标签显示真实值
    expand = expansion(c(0, 0))
  ) +
  theme_prism() +
  theme(axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

print(p3b)

# 组合保存（取消注释）
# fig3ab <- p3a + p3b + plot_layout(guides = "collect")
# ggsave("results/RNAseq_analysis/Figures/heatmap、kegg、go，/Fig3AB_GO_ORA_RoundBar.pdf",
#        fig3ab, width = 14, height = 9, dpi = 300)



# ============================================================
# 重建 Pathway 级别 KEGG term2gene + 重跑 ORA + 交叉验证
# ============================================================
library(clusterProfiler)
library(dplyr)
library(ggplot2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# -----------------------------------------------------------
# 1. 从 eggNOG 重建 Pathway 级别的 kegg_term2gene
# -----------------------------------------------------------
cat("========== 1. 重建 Pathway 级别 KEGG term2gene ==========\n")

eggnog <- read.delim(
  "annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular",
  header = TRUE, sep = "\t", stringsAsFactors = FALSE,
  quote = "", comment.char = "", check.names = FALSE
)

# ID 映射（CAC → gene-MCOR，同你 GSEA 代码）
ft <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
id_map <- ft[, c("product_accession", "locus_tag")]
id_map <- id_map[id_map$product_accession != "" & !is.na(id_map$product_accession) &
                   id_map$locus_tag != "" & !is.na(id_map$locus_tag), ]
id_map <- id_map[!duplicated(id_map$product_accession), ]
id_map$gene_id <- paste0("gene-", id_map$locus_tag)

m_idx <- match(eggnog[, 1], id_map[["product_accession"]])
eggnog[["gene_id"]] <- id_map[["gene_id"]][m_idx]
eggnog <- eggnog[!is.na(eggnog$gene_id), ]

# 构建 Pathway 级别 TERM2GENE（只取 map 前缀）
kegg_list <- strsplit(eggnog$KEGG_Pathway, ",")
kegg_pathway_df <- data.frame(
  term = unlist(kegg_list),
  gene = rep(eggnog$gene_id, sapply(kegg_list, length)),
  stringsAsFactors = FALSE
)
kegg_pathway_df <- kegg_pathway_df[
  kegg_pathway_df$term != "" & 
    kegg_pathway_df$term != "-" & 
    !is.na(kegg_pathway_df$term), 
]
kegg_pathway_df <- unique(kegg_pathway_df)
kegg_pathway_df <- kegg_pathway_df[grepl("^map", kegg_pathway_df$term), ]

cat("Pathway 级别 KEGG term2gene:", nrow(kegg_pathway_df), "对\n")
cat("唯一 Pathway 数:", length(unique(kegg_pathway_df$term)), "\n")

# 保存新文件（不覆盖旧的，以防万一）
saveRDS(kegg_pathway_df, "eggnog_annotation/kegg_term2gene_pathway.rds")
cat(">>> 已保存: eggnog_annotation/kegg_term2gene_pathway.rds\n")

# -----------------------------------------------------------
# 2. 用 Pathway 级别重跑 KEGG ORA（三个比较组）
# -----------------------------------------------------------
cat("\n========== 2. 重跑 Pathway 级别 KEGG ORA ==========\n")

bg_genes <- rownames(read.table("results/counts/counts_s0.txt", header = TRUE, row.names = 1, sep = "\t"))

run_kegg_ora_pathway <- function(deg_file, comparison_name) {
  
  cat("\n---", comparison_name, "---\n")
  
  deg_df <- read.table(deg_file, header = TRUE, stringsAsFactors = FALSE)
  gene_list <- unique(as.character(deg_df[[1]]))
  cat("Consensus DEG:", length(gene_list), "\n")
  
  kegg_enrich <- enricher(
    gene          = gene_list,
    universe      = bg_genes,
    TERM2GENE     = kegg_pathway_df,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  if (!is.null(kegg_enrich) && nrow(as.data.frame(kegg_enrich)) > 0) {
    df <- as.data.frame(kegg_enrich)
    write.csv(df, paste0("results/RNAseq_analysis/Results/KEGG_ORA_Pathway_", comparison_name, ".csv"), row.names = FALSE)
    cat("显著 Pathway:", nrow(df), "\n")
    cat("ID:", paste(df$ID, collapse = ", "), "\n")
    return(kegg_enrich)
  } else {
    cat("KEGG Pathway ORA: 无显著富集\n")
    return(NULL)
  }
}

ora_kegg_t24  <- run_kegg_ora_pathway("results/RNAseq_analysis/Results/HighConf_DEG_T24vsCtrl.txt", "T24vsCtrl")
ora_kegg_t32  <- run_kegg_ora_pathway("results/RNAseq_analysis/Results/HighConf_DEG_T32vsCtrl.txt", "T32vsCtrl")
ora_kegg_t32t24 <- run_kegg_ora_pathway("results/RNAseq_analysis/Results/HighConf_DEG_T32vsT24.txt", "T32vsT24")

# -----------------------------------------------------------
# 3. KEGG ORA vs GSEA 交叉验证（Pathway 级别）
# -----------------------------------------------------------
cat("\n========== 3. KEGG Pathway 级别交叉验证 ==========\n")

# 加载 GSEA 对象（Pathway 级别，你之前跑的）
load("results/RNAseq_analysis/Results/kegg_gsea_T24vCtrl.RData")
load("results/RNAseq_analysis/Results/kegg_gsea_T32vCtrl.RData")
load("results/RNAseq_analysis/Results/kegg_gsea_T32vT24.RData")

cv_kegg <- function(ora_obj, gsea_obj, comparison_name) {
  
  if (is.null(ora_obj)) {
    cat(comparison_name, ": ORA 无结果，跳过\n")
    return(list(dual = character(), ora_only = character(), gsea_only = character()))
  }
  
  ora_df <- as.data.frame(ora_obj)
  gsea_df <- as.data.frame(gsea_obj)
  
  ora_sig <- ora_df$ID[ora_df$p.adjust < 0.05]
  gsea_sig <- gsea_df$ID[gsea_df$p.adjust < 0.05]
  dual <- intersect(ora_sig, gsea_sig)
  
  cat("\n", comparison_name, ":\n", sep = "")
  cat("  ORA 显著 Pathway:", length(ora_sig), "\n")
  cat("  GSEA 显著 Pathway:", length(gsea_sig), "\n")
  cat("  双阳性:", length(dual), "\n")
  
  if (length(dual) > 0) {
    dual_df <- ora_df[ora_df$ID %in% dual, c("ID", "Description", "p.adjust", "Count")]
    gsea_match <- gsea_df[gsea_df$ID %in% dual, c("ID", "NES", "p.adjust")]
    names(gsea_match)[2:3] <- c("GSEA_NES", "GSEA_p.adjust")
    dual_df <- merge(dual_df, gsea_match, by = "ID")
    write.csv(dual_df, paste0("results/RNAseq_analysis/Results/Dual_Positive_KEGG_Pathway_", comparison_name, ".csv"), row.names = FALSE)
    cat("  -> 已保存双阳性结果\n")
    print(dual_df)
  }
  
  return(list(
    dual = dual,
    ora_only = setdiff(ora_sig, gsea_sig),
    gsea_only = setdiff(gsea_sig, ora_sig)
  ))
}

cv_t24_kegg  <- cv_kegg(ora_kegg_t24, kegg_gsea_T24vCtrl, "T24vsCtrl")
cv_t32_kegg  <- cv_kegg(ora_kegg_t32, kegg_gsea_T32vCtrl, "T32vsCtrl")
cv_t32t24_kegg <- cv_kegg(ora_kegg_t32t24, kegg_gsea_T32vT24, "T32vsT24")

cat("\n========== 汇总 ==========\n")
cat("T24vsCtrl  KEGG 双阳性:", length(cv_t24_kegg$dual), "\n")
cat("T32vsCtrl  KEGG 双阳性:", length(cv_t32_kegg$dual), "\n")
cat("T32vsT24 KEGG 双阳性:", length(cv_t32t24_kegg$dual), "\n")













# ============================================================
# KEGG ORA 美化 — keggGet 逐个查询版（已验证网络通）
# 原则：只 preview，不 ggsave
# ============================================================
library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)
library(KEGGREST)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 辅助函数 ==========
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(y) as.numeric(y[1])/as.numeric(y[2]))
}

wrap_name <- function(x, width = 50) {
  sapply(x, function(s) str_wrap(s, width = width))
}

# ========== KEGG 名称查询（带重试） ==========
kegg_get_with_retry <- function(query_id, max_retry = 3) {
  for (i in 1:max_retry) {
    info <- tryCatch(keggGet(query_id), error = function(e) NULL)
    if (!is.null(info) && length(info) > 0) return(info)
    Sys.sleep(0.3)
  }
  return(NULL)
}

get_kegg_name <- function(map_id) {
  clean <- gsub("^(map|ko|path:ko|path:map|ko:)", "", map_id)
  
  # 先试 map 前缀
  info <- kegg_get_with_retry(paste0("map", clean))
  if (is.null(info)) {
    info <- kegg_get_with_retry(paste0("ko", clean))
  }
  
  if (!is.null(info) && !is.null(info[[1]]$NAME)) {
    name <- info[[1]]$NAME[1]
    name <- gsub("\\s*\\[.*\\]$", "", name)
    return(name)
  }
  return(map_id)  # 失败返回原ID
}

# ========== 通用画图函数 ==========
plot_kegg_ora <- function(csv_path, title, color_low, color_high, top_n = 15) {
  
  if (!file.exists(csv_path)) {
    cat("文件不存在:", csv_path, "\n")
    return(NULL)
  }
  
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  if (nrow(df) == 0) {
    cat(title, ": 无显著通路\n")
    return(NULL)
  }
  
  # 逐个查名称
  cat("\n===", title, "— 查询通路名称 ===\n")
  df$DisplayName <- sapply(df$ID, get_kegg_name)
  
  # 如果查询失败（还是编号），用 Description 兜底
  is_failed <- df$DisplayName == df$ID
  if (any(is_failed) && "Description" %in% names(df)) {
    df$DisplayName[is_failed] <- ifelse(
      df$Description[is_failed] != "" & !is.na(df$Description[is_failed]),
      df$Description[is_failed],
      df$DisplayName[is_failed]
    )
  }
  
  df$DisplayName <- wrap_name(df$DisplayName, 55)
  
  n_show <- min(top_n, nrow(df))
  df_top <- df %>%
    arrange(p.adjust) %>%
    head(n_show) %>%
    mutate(
      DisplayName = factor(DisplayName, levels = rev(DisplayName)),
      logp = -log10(p.adjust),
      gene_ratio = parse_ratio(GeneRatio)
    )
  
  cat(title, ": 共", nrow(df), "条显著通路，显示 Top", n_show, "\n")
  
  p <- ggplot(df_top, aes(x = gene_ratio, y = DisplayName)) +
    geom_point(aes(size = Count, color = logp), alpha = 0.9) +
    scale_color_gradient(low = color_low, high = color_high, name = "-log10(FDR)") +
    scale_size_continuous(range = c(3, 10), name = "Gene Count") +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, color = color_high),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# ========== 三个比较组 ==========
p_t32 <- plot_kegg_ora(
  csv_path = "results/RNAseq_analysis/Results/KEGG_ORA_T32vsCtrl.csv",
  title = "KEGG ORA — T32 vs Ctrl",
  color_low = "#FFE082", color_high = "#E65100", top_n = 15
)
if (!is.null(p_t32)) print(p_t32)

p_t24 <- plot_kegg_ora(
  csv_path = "results/RNAseq_analysis/Results/KEGG_ORA_T24vsCtrl.csv",
  title = "KEGG ORA — T24 vs Ctrl",
  color_low = "#90CAF9", color_high = "#1565C0", top_n = 10
)
if (!is.null(p_t24)) print(p_t24)

p_t32t24 <- plot_kegg_ora(
  csv_path = "results/RNAseq_analysis/Results/KEGG_ORA_T32vsT24.csv",
  title = "KEGG ORA — T32 vs T24",
  color_low = "#A5D6A7", color_high = "#2E7D32", top_n = 50
)
if (!is.null(p_t32t24)) print(p_t32t24)

# 组合预览（T32 + T24 并排）
if (!is.null(p_t32) && !is.null(p_t24)) {
  p_combined <- p_t32 + p_t24 + plot_layout(guides = "collect")
  print(p_combined)
}

cat("\n>>> KEGG ORA 预览完成\n")


# ============================================================
# 完整自包含版：GO + KEGG 交叉验证散点图（6个）
# 复制整块，直接粘贴运行
# ============================================================
library(ggplot2)
library(dplyr)
library(ggrepel)
library(KEGGREST)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ---------- 0. 名称缓存（必须最先定义）----------
name_cache <- new.env()

# ---------- 1. 手动保底映射表（已补全所有 unknown）----------
manual_go_map <- list(
  "GO:0042026" = "protein refolding",
  "GO:0031929" = "TOR signaling",
  "GO:0038202" = "TORC1 signaling",
  "GO:0032781" = "positive regulation of ATP-dependent activity",
  "GO:0043462" = "regulation of ATP-dependent activity",
  "GO:0046889" = "positive regulation of lipid biosynthetic process",
  "GO:1901800" = "positive regulation of proteasomal protein catabolic process",
  "GO:0032436" = "positive regulation of proteasomal ubiquitin-dependent protein catabolic process",
  "GO:0061077" = "obsolete chaperone-mediated protein folding",
  "GO:0051082" = "obsolete unfolded protein binding",
  "GO:0043620" = "obsolete regulation of DNA-templated transcription in response to stress",
  "GO:0043618" = "obsolete regulation of transcription from RNA polymerase II promoter in response to stress",
  "GO:1903052" = "obsolete positive regulation of proteolysis involved in protein catabolic process",
  "GO:0006457" = "protein folding",
  "GO:1903599" = "positive regulation of autophagy of mitochondrion",
  "GO:0006986" = "response to unfolded protein",
  "GO:0031072" = "heat shock protein binding",
  "GO:0060590" = "ATPase regulator activity"
)

# ---------- 2. GO 名称查询（修复版：GOTERM + @Term）----------
get_go_name <- function(go_id) {
  if (exists(go_id, envir = name_cache)) return(get(go_id, envir = name_cache))
  
  name <- NULL
  
  # 方法1: GO.db::GOTERM
  tryCatch({
    if (requireNamespace("GO.db", quietly = TRUE)) {
      terms <- AnnotationDbi::mget(go_id, GO.db::GOTERM, ifnotfound = NA)
      term <- terms[[1]]
      if (!is.na(term)) name <- term@Term
    }
  }, error = function(e) {})
  
  # 方法2: 手动保底
  if (is.null(name) && go_id %in% names(manual_go_map)) {
    name <- manual_go_map[[go_id]]
  }
  
  # 方法3: 在线查询 EBI QuickGO
  if (is.null(name)) {
    tryCatch({
      url <- paste0("https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/", go_id)
      r <- httr::GET(url, httr::add_headers(Accept = "application/json"), httr::timeout(5))
      if (r$status_code == 200) {
        d <- jsonlite::fromJSON(httr::content(r, "text"))
        if ("results" %in% names(d) && nrow(d$results) > 0) {
          name <- d$results$name[1]
          if (d$results$isObsolete[1]) name <- paste0("obsolete ", name)
        }
      }
    }, error = function(e) {})
  }
  
  if (is.null(name)) name <- paste0(go_id, " (unknown)")
  
  assign(go_id, name, envir = name_cache)
  return(name)
}

# ---------- 3. KEGG 名称查询 ----------
get_kegg_name <- function(map_id) {
  if (exists(map_id, envir = name_cache)) return(get(map_id, envir = name_cache))
  clean <- gsub("^map", "", map_id)
  info <- tryCatch(keggGet(paste0("map", clean)), error = function(e) NULL)
  name <- map_id
  if (!is.null(info) && length(info) > 0 && !is.null(info[[1]]$NAME)) {
    name <- gsub("\\s*\\[.*\\]$", "", info[[1]]$NAME[1])
  }
  assign(map_id, name, envir = name_cache)
  return(name)
}

get_name <- function(id) {
  if (grepl("^GO:", id)) return(get_go_name(id))
  if (grepl("^map", id)) return(get_kegg_name(id))
  return(id)
}

# ---------- 4. 排除列表 ----------
EXCLUDE_KEYWORDS <- c(
  "leukemia", "myeloid", "collecting duct", "pertussis", "viral myocarditis",
  "osteoclast", "Epstein-Barr", "EBV", "measles", "lung cancer", "HTLV",
  "Toxoplasmosis", "toxoplasmosis", "Legionellosis", "legionellosis",
  "Th17 cell differentiation", "endometrial cancer", "cancer", "carcinoma",
  "lymphoma", "virus infection", "viral", "small cell", "chronic myeloid",
  "acute myeloid", "bladder cancer", "glioma", "melanoma", "tuberculosis",
  "malaria", "amoebiasis", "Chagas", "African trypanosomiasis", "Leishmaniasis",
  "Shigellosis", "Salmonella infection", "Pathogenic Escherichia", "Helicobacter",
  "Vibrio cholerae", "Staphylococcus aureus", "hepatitis", "influenza",
  "papillomavirus", "cytomegalovirus", "herpes", "HIV", "COVID-19", "coronavirus",
  "Kaposi", "Epstein", "HTLV", "measles", "pertussis", "autoimmune thyroid",
  "systemic lupus", "rheumatoid arthritis", "inflammatory bowel", "asthma",
  "allograft rejection", "graft-versus-host", "primary immunodeficiency",
  "type I diabetes", "type II diabetes", "maturity onset diabetes", "Cushing",
  "alcoholic liver", "non-alcoholic fatty liver", "insulin resistance",
  "adipocytokine", "thyroid cancer", "prostate cancer", "pancreatic cancer",
  "renal cell carcinoma", "colorectal cancer", "gastric cancer", "hepatocellular",
  "breast cancer", "basal cell carcinoma", "PD-L1", "cocaine", "amphetamine",
  "morphine", "nicotine", "alcoholism", "Alzheimer", "Parkinson", "Huntington",
  "amyotrophic lateral sclerosis", "prion", "cocaine addiction", "long-term potentiation",
  "synaptic", "neurotrophin", "GABAergic", "dopaminergic", "serotonergic",
  "cholinergic", "glutamatergic", "olfactory", "taste", "phototransduction",
  "salivary secretion", "gastric acid", "pancreatic secretion", "bile secretion",
  "vitamin digestion", "mineral absorption", "prolactin", "thyroid hormone",
  "ovarian steroidogenesis", "progesterone", "estrogen signaling", "GnRH",
  "oxytocin", "glucagon", "renin", "aldosterone", "relaxin", "cortisol",
  "parathyroid hormone", "adipocytokine", "melanogenesis", "GnRH secretion",
  "cardiomyopathy", "arrhythmogenic", "dilated cardiomyopathy", "hypertrophic",
  "diabetic cardiomyopathy", "viral myocarditis", "atherosclerosis", "fluid shear",
  "longevity regulating", "ferroptosis", "necroptosis", "cellular senescence",
  "oocyte meiosis", "meiosis", "spermatogenesis", "adherens junction",
  "tight junction", "gap junction", "focal adhesion", "ECM-receptor",
  "cell adhesion molecules", "axon guidance", "osteoclast", "dorso-ventral",
  "Wnt signaling", "Notch", "Hedgehog", "TGF-beta", "Hippo", "VEGF",
  "apelin", "FoxO", "NF-kappa B", "chemokine", "cytokine-cytokine receptor",
  "T cell receptor", "B cell receptor", "Fc epsilon", "Fc gamma",
  "natural killer cell", "hematopoietic cell lineage", "leukocyte transendothelial",
  "intestinal immune network", "complement and coagulation", "platelet activation",
  "neutrophil extracellular trap", "Toll-like receptor", "NOD-like receptor",
  "RIG-I-like receptor", "cytosolic DNA-sensing", "C-type lectin",
  "plant-pathogen interaction", "striated muscle"
)

should_exclude <- function(name) {
  any(sapply(EXCLUDE_KEYWORDS, function(k) grepl(k, name, ignore.case = TRUE)))
}

# ---------- 5. 主画图函数 ----------
make_cv_scatter <- function(ora_csv, gsea_input, title_text, 
                            highlight_ids = NULL, highlight_names = NULL) {
  
  ora_df <- read.csv(ora_csv, stringsAsFactors = FALSE)
  
  if (is.character(gsea_input) && file.exists(gsea_input)) {
    env_tmp <- new.env()
    load(gsea_input, envir = env_tmp)
    obj_name <- ls(env_tmp)[1]
    gsea_obj <- get(obj_name, envir = env_tmp)
  } else {
    gsea_obj <- gsea_input
  }
  gsea_df <- as.data.frame(gsea_obj)
  
  merged <- merge(
    ora_df[, c("ID", "Description", "p.adjust", "Count")],
    gsea_df[, c("ID", "NES", "p.adjust")],
    by = "ID", all = TRUE,
    suffixes = c("_ORA", "_GSEA")
  )
  
  merged$p.adjust_ORA[is.na(merged$p.adjust_ORA)] <- 1
  merged$x_val <- -log10(merged$p.adjust_ORA)
  merged$NES[is.na(merged$NES)] <- 0
  merged$p.adjust_GSEA[is.na(merged$p.adjust_GSEA)] <- 1
  merged$Count[is.na(merged$Count)] <- 0
  
  ora_sig  <- merged$p.adjust_ORA < 0.05
  gsea_sig <- merged$p.adjust_GSEA < 0.05
  
  merged$category <- "Neither"
  merged$category[ora_sig & !gsea_sig] <- "ORA only"
  merged$category[!ora_sig & gsea_sig] <- "GSEA only"
  merged$category[ora_sig & gsea_sig]  <- "Dual-positive"
  merged$category <- factor(merged$category,
                            levels = c("Dual-positive", "ORA only", "GSEA only", "Neither"))
  
  tbl <- table(merged$category)
  cat("\n===", title_text, "===\n")
  print(tbl)
  
  bg <- merged[merged$category %in% c("Neither", "GSEA only"), ]
  fg <- merged[merged$category %in% c("Dual-positive", "ORA only"), ]
  
  set.seed(42)
  gsea_only <- bg[bg$category == "GSEA only", ]
  if (nrow(gsea_only) > 150) gsea_only <- gsea_only[sample(nrow(gsea_only), 150), ]
  neither <- bg[bg$category == "Neither", ]
  if (nrow(neither) > 100) neither <- neither[sample(nrow(neither), 100), ]
  bg_plot <- rbind(gsea_only, neither)
  
  p <- ggplot() +
    geom_point(data = bg_plot,
               aes(x = x_val, y = NES, color = category),
               size = 2.8, alpha = 0.38, stroke = 0,
               position = position_jitter(width = 0.1, height = 0, seed = 123)) +
    geom_point(data = fg,
               aes(x = x_val, y = NES, color = category, size = Count),
               alpha = 0.55, stroke = 0.1) +
    scale_color_manual(values = c("Dual-positive" = "#C74375", "ORA only" = "#3070B3",
                                  "GSEA only" = "#FFCD00", "Neither" = "#6A757E"),
                       drop = FALSE) +
    scale_size_continuous(range = c(3.5, 9), name = "Gene Count") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
    coord_cartesian(xlim = c(-0.5, max(merged$x_val, na.rm = TRUE) * 1.05),
                    ylim = c(min(merged$NES, na.rm = TRUE) * 1.15,
                             max(merged$NES, na.rm = TRUE) * 1.15)) +
    labs(x = expression(-log[10]~(ORA~FDR)), y = "GSEA NES",
         color = "Category", title = title_text) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          panel.grid.minor = element_blank(),
          legend.key.size = unit(0.6, "cm"))
  
  # 标注
  dual_df <- merged[merged$category == "Dual-positive", ]
  
  if (nrow(dual_df) > 0) {
    dual_df$PathwayName <- sapply(dual_df$ID, get_name)
    
    # 过滤 obsolete
    dual_df <- dual_df[!grepl("^obsolete", dual_df$PathwayName, ignore.case = TRUE), ]
    
    cat("\n双阳性通路 Top 10:\n")
    print(dual_df[order(dual_df$NES, decreasing = TRUE), c("ID", "PathwayName", "NES")] %>% head(10))
    
    if (!is.null(highlight_ids)) {
      for (i in seq_along(highlight_ids)) {
        hid <- highlight_ids[i]
        hname <- ifelse(is.null(highlight_names), hid, highlight_names[i])
        if (hid %in% dual_df$ID) {
          cat("✅", hname, "在双阳性中，NES =", dual_df$NES[dual_df$ID == hid], "\n")
        } else if (hid %in% merged$ID) {
          cat("⚠️ ", hname, "不在双阳性中，状态:", as.character(merged$category[merged$ID == hid]), "\n")
        } else {
          cat("❌", hname, "不在当前数据集中\n")
        }
      }
    }
    
    label_ids <- character()
    if (!is.null(highlight_ids)) {
      hit <- highlight_ids[highlight_ids %in% dual_df$ID]
      label_ids <- c(label_ids, hit)
    }
    
    remaining <- dual_df[!dual_df$ID %in% label_ids, ] %>% arrange(desc(NES))
    remaining$exclude <- sapply(remaining$PathwayName, should_exclude)
    remaining_clean <- remaining[!remaining$exclude, ]
    
    n_needed <- max(0, 5 - length(label_ids))
    if (n_needed > 0 && nrow(remaining_clean) > 0) {
      label_ids <- c(label_ids, head(remaining_clean$ID, n_needed))
    }
    
    if (length(label_ids) < 5 && nrow(remaining) > 0) {
      backup <- remaining[!remaining$ID %in% label_ids, ]
      n_backup <- min(5 - length(label_ids), nrow(backup))
      label_ids <- c(label_ids, head(backup$ID, n_backup))
    }
    
    label_df <- dual_df[dual_df$ID %in% label_ids, ]
    label_df <- label_df[order(label_df$NES, decreasing = TRUE), ]
    
    cat("\n最终标注:\n")
    print(label_df[, c("ID", "PathwayName", "NES")])
    
    if (nrow(label_df) > 0) {
      p <- p + geom_text_repel(
        data = label_df,
        aes(x = x_val, y = NES, label = PathwayName),
        size = 3.3, color = "black",
        max.overlaps = Inf, force = 8, force_pull = 0.5,
        box.padding = 0.5, point.padding = 0.4,
        segment.size = 0.3, fontface = "italic",
        min.segment.length = 0
      )
    }
  }
  
  p <- p + annotate(
    "text", x = Inf, y = -Inf,
    label = paste(names(tbl), tbl, sep = ": ", collapse = "\n"),
    hjust = 1, vjust = 0, size = 3.5, color = "grey40", fontface = "italic"
  )
  
  return(p)
}

# ============================================================
# 运行 6 个散点图
# ============================================================

cat("========== GO 交叉验证 ==========\n")

load("results/RNAseq_analysis/Results/go_gsea_T24vCtrl.RData")
p_t24_go <- make_cv_scatter(
  "results/RNAseq_analysis/Results/GO_ORA_T24vsCtrl.csv",
  go_gsea_T24vCtrl, "GO Cross-Validation — T24 vs Ctrl"
)

load("results/RNAseq_analysis/Results/go_gsea_T32vCtrl.RData")
p_t32_go <- make_cv_scatter(
  "results/RNAseq_analysis/Results/GO_ORA_T32vsCtrl.csv",
  go_gsea_T32vCtrl, "GO Cross-Validation — T32 vs Ctrl"
)

load("results/RNAseq_analysis/Results/go_gsea_T32vT24.RData")
p_t32t24_go <- make_cv_scatter(
  "results/RNAseq_analysis/Results/GO_ORA_T32vsT24.csv",
  go_gsea_T32vT24, "GO Cross-Validation — T32 vs T24"
)

cat("\n========== KEGG 交叉验证 ==========\n")

p_t24_k <- make_cv_scatter(
  "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T24vsCtrl.csv",
  "results/RNAseq_analysis/Results/kegg_gsea_T24vCtrl.RData",
  "KEGG Cross-Validation — T24 vs Ctrl",
  highlight_ids = c("map04141", "map04210", "map03008"),
  highlight_names = c("Protein processing in ER", "Apoptosis", "Ribosome biogenesis")
)

p_t32_k <- make_cv_scatter(
  "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T32vsCtrl.csv",
  "results/RNAseq_analysis/Results/kegg_gsea_T32vCtrl.RData",
  "KEGG Cross-Validation — T32 vs Ctrl",
  highlight_ids = c("map04141", "map04210", "map03008"),
  highlight_names = c("Protein processing in ER", "Apoptosis", "Ribosome biogenesis")
)

p_t32t24_k <- make_cv_scatter(
  "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T32vsT24.csv",
  "results/RNAseq_analysis/Results/kegg_gsea_T32vT24.RData",
  "KEGG Cross-Validation — T32 vs T24",
  highlight_ids = c("map04141", "map04210", "map03008"),
  highlight_names = c("Protein processing in ER", "Apoptosis", "Ribosome biogenesis")
)

print(p_t24_go)
print(p_t32_go)
print(p_t32t24_go)
print(p_t24_k)
print(p_t32_k)
print(p_t32t24_k)

cat("\n>>> 6 个散点图全部完成\n")



# ============================================================
# 完整自包含版：KEGG 桑基图 + 气泡图（独立预览）
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(KEGGREST)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ---------- 安装/加载 ggsankey（只需一次）----------
if (!requireNamespace("ggsankey", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("davidsjoberg/ggsankey")
}
library(ggsankey)

# ============================================================
# 1. 排除列表 + 判断函数
# ============================================================
EXCLUDE_KEYWORDS <- c(
  "leukemia", "myeloid", "collecting duct", "pertussis", "viral myocarditis",
  "osteoclast", "Epstein-Barr", "EBV", "measles", "lung cancer", "HTLV",
  "Toxoplasmosis", "toxoplasmosis", "Legionellosis", "legionellosis",
  "Th17 cell differentiation", "endometrial cancer", "cancer", "carcinoma",
  "lymphoma", "virus infection", "viral", "small cell", "chronic myeloid",
  "acute myeloid", "bladder cancer", "glioma", "melanoma", "tuberculosis",
  "malaria", "amoebiasis", "Chagas", "African trypanosomiasis", "Leishmaniasis",
  "Shigellosis", "Salmonella infection", "Pathogenic Escherichia", "Helicobacter",
  "Vibrio cholerae", "Staphylococcus aureus", "hepatitis", "influenza",
  "papillomavirus", "cytomegalovirus", "herpes", "HIV", "COVID-19", "coronavirus",
  "Kaposi", "Epstein", "HTLV", "measles", "pertussis", "autoimmune thyroid",
  "systemic lupus", "rheumatoid arthritis", "inflammatory bowel", "asthma",
  "allograft rejection", "graft-versus-host", "primary immunodeficiency",
  "type I diabetes", "type II diabetes", "maturity onset diabetes", "Cushing",
  "alcoholic liver", "non-alcoholic fatty liver", "insulin resistance",
  "adipocytokine", "thyroid cancer", "prostate cancer", "pancreatic cancer",
  "renal cell carcinoma", "colorectal cancer", "gastric cancer", "hepatocellular",
  "breast cancer", "basal cell carcinoma", "PD-L1", "cocaine", "amphetamine",
  "morphine", "nicotine", "alcoholism", "Alzheimer", "Parkinson", "Huntington",
  "amyotrophic lateral sclerosis", "prion", "cocaine addiction", "long-term potentiation",
  "synaptic", "neurotrophin", "GABAergic", "dopaminergic", "serotonergic",
  "cholinergic", "glutamatergic", "olfactory", "taste", "phototransduction",
  "salivary secretion", "gastric acid", "pancreatic secretion", "bile secretion",
  "vitamin digestion", "mineral absorption", "prolactin", "thyroid hormone",
  "ovarian steroidogenesis", "progesterone", "estrogen signaling", "GnRH",
  "oxytocin", "glucagon", "renin", "aldosterone", "relaxin", "cortisol",
  "parathyroid hormone", "adipocytokine", "melanogenesis", "GnRH secretion",
  "cardiomyopathy", "arrhythmogenic", "dilated cardiomyopathy", "hypertrophic",
  "diabetic cardiomyopathy", "viral myocarditis", "atherosclerosis", "fluid shear",
  "longevity regulating", "ferroptosis", "necroptosis", "cellular senescence",
  "oocyte meiosis", "meiosis", "spermatogenesis", "adherens junction",
  "tight junction", "gap junction", "focal adhesion", "ECM-receptor",
  "cell adhesion molecules", "axon guidance", "osteoclast", "dorso-ventral",
  "Wnt signaling", "Notch", "Hedgehog", "TGF-beta", "Hippo", "VEGF",
  "apelin", "FoxO", "NF-kappa B", "chemokine", "cytokine-cytokine receptor",
  "T cell receptor", "B cell receptor", "Fc epsilon", "Fc gamma",
  "natural killer cell", "hematopoietic cell lineage", "leukocyte transendothelial",
  "intestinal immune network", "complement and coagulation", "platelet activation",
  "neutrophil extracellular trap", "Toll-like receptor", "NOD-like receptor",
  "RIG-I-like receptor", "cytosolic DNA-sensing", "C-type lectin",
  "plant-pathogen interaction", "striated muscle"
)

should_exclude <- function(name) {
  any(sapply(EXCLUDE_KEYWORDS, function(k) grepl(k, name, ignore.case = TRUE)))
}

# ============================================================
# 2. KEGG 名称查询
# ============================================================
kegg_name_cache <- new.env()

get_kegg_name <- function(map_id) {
  if (exists(map_id, envir = kegg_name_cache)) return(get(map_id, envir = kegg_name_cache))
  clean <- gsub("^map", "", map_id)
  info <- tryCatch(keggGet(paste0("map", clean)), error = function(e) NULL)
  name <- map_id
  if (!is.null(info) && length(info) > 0 && !is.null(info[[1]]$NAME)) {
    name <- gsub("\\s*\\[.*\\]$", "", info[[1]]$NAME[1])
  }
  assign(map_id, name, envir = kegg_name_cache)
  return(name)
}

# ============================================================
# 3. 基因名映射（gene-MCOR_xxx → 基因符号/描述）
# ============================================================
build_gene_name_map <- function() {
  cat(">>> 正在构建基因名映射表...\n")
  
  ft <- read.delim(
    "annotation/GCA_011752425.2_MCOR1.1_feature_table.txt",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
  )
  
  id_map <- ft[, c("product_accession", "locus_tag")]
  id_map <- id_map[
    id_map$product_accession != "" & !is.na(id_map$product_accession) &
      id_map$locus_tag != "" & !is.na(id_map$locus_tag),
  ]
  id_map <- id_map[!duplicated(id_map$product_accession), ]
  
  id_map$gene_id_dash  <- paste0("gene-", id_map$locus_tag)
  id_map$gene_id_under <- paste0("gene_", id_map$locus_tag)
  
  eggnog <- read.delim(
    "annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular",
    header = TRUE, sep = "\t", stringsAsFactors = FALSE,
    quote = "", comment.char = "", check.names = FALSE
  )
  
  m <- match(eggnog[, 1], id_map$product_accession)
  
  pref_col <- grep("Preferred", names(eggnog), value = TRUE, ignore.case = TRUE)
  desc_col <- grep("Description", names(eggnog), value = TRUE, ignore.case = TRUE)
  
  name_map <- list()
  
  for (i in 1:nrow(eggnog)) {
    if (is.na(m[i])) next
    gid_dash  <- id_map$gene_id_dash[m[i]]
    gid_under <- id_map$gene_id_under[m[i]]
    
    display_name <- NULL
    if (length(pref_col) > 0) {
      pref <- eggnog[[pref_col[1]]][i]
      if (!is.na(pref) && pref != "" && pref != "-") display_name <- pref
    }
    if (is.null(display_name) && length(desc_col) > 0) {
      desc <- eggnog[[desc_col[1]]][i]
      if (!is.na(desc) && desc != "" && desc != "-") display_name <- substr(desc, 1, 35)
    }
    if (is.null(display_name)) display_name <- gid_dash
    
    name_map[[gid_dash]]  <- display_name
    name_map[[gid_under]] <- display_name
  }
  
  cat(">>> 映射表完成，共", length(unique(name_map)), "条唯一名称\n")
  return(name_map)
}

gene_name_map <- build_gene_name_map()

get_gene_display_name <- function(gene_id) {
  if (gene_id %in% names(gene_name_map)) return(gene_name_map[[gene_id]])
  return(gene_id)
}

# ============================================================
# 4. 桑基图函数
# ============================================================
plot_kegg_sankey <- function(
    ora_csv,
    deg_file,
    title,
    top_n = 10,
    max_genes_per_pathway = 10
) {
  
  kegg <- read.csv(ora_csv, stringsAsFactors = FALSE)
  
  if (!"geneID" %in% names(kegg)) {
    stop("CSV 里没有 geneID 列！")
  }
  
  kegg$PathwayName <- sapply(kegg$ID, get_kegg_name)
  
  kegg$exclude <- sapply(kegg$PathwayName, should_exclude)
  kegg_clean <- kegg[!kegg$exclude, ]
  
  n_excluded <- sum(kegg$exclude)
  cat(">>> [桑基图] 排除", n_excluded, "条无关通路，剩余", nrow(kegg_clean), "条\n")
  
  if (nrow(kegg_clean) == 0) stop("过滤后无通路可用")
  
  kegg_clean <- kegg_clean[order(kegg_clean$p.adjust), ]
  kegg_top <- head(kegg_clean, top_n)
  
  cat(">>> [桑基图] 最终展示 Top", nrow(kegg_top), "条通路\n")
  
  kegg_long <- kegg_top %>%
    separate_rows(geneID, sep = "/") %>%
    filter(geneID != "" & !is.na(geneID))
  
  deg <- read.csv(deg_file, row.names = 1, stringsAsFactors = FALSE)
  deg$geneID <- rownames(deg)
  
  fc_col <- ifelse("log2FoldChange" %in% names(deg), "log2FoldChange", "logFC")
  
  sankey_data <- merge(
    kegg_long,
    deg[, c("geneID", fc_col)],
    by = "geneID",
    all.x = TRUE
  )
  names(sankey_data)[names(sankey_data) == fc_col] <- "logFC"
  
  sankey_data <- sankey_data %>%
    group_by(PathwayName) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = max_genes_per_pathway) %>%
    ungroup()
  
  sankey_data$GeneDisplay <- sapply(sankey_data$geneID, get_gene_display_name)
  
  cat("[桑基图]", title, ": 展示", length(unique(sankey_data$GeneDisplay)), "个基因\n")
  
  df_long <- sankey_data %>%
    make_long(GeneDisplay, PathwayName)
  
  cols_source <- c("#3070B3", "#eeeeee", "#FFCD00")
  limit_val <- max(abs(sankey_data$logFC), na.rm = TRUE)
  limit_range <- c(-limit_val, limit_val)
  
  pathway_means <- sankey_data %>%
    group_by(PathwayName) %>%
    summarise(mean_logFC = mean(logFC, na.rm = TRUE), .groups = "drop")
  
  col_func <- colorRampPalette(cols_source)
  n_col <- 100
  
  pathway_colors <- pathway_means %>%
    mutate(
      idx = as.numeric(cut(mean_logFC, breaks = n_col, labels = FALSE)),
      idx = pmin(pmax(idx, 1), n_col),
      color = col_func(n_col)[idx]
    ) %>%
    with(setNames(color, PathwayName))
  
  gene_idx <- as.numeric(cut(sankey_data$logFC, breaks = n_col, labels = FALSE))
  gene_idx <- pmin(pmax(gene_idx, 1), n_col)
  gene_colors <- setNames(col_func(n_col)[gene_idx], sankey_data$GeneDisplay)
  
  all_colors <- c(gene_colors, pathway_colors)
  df_long$fill_color <- all_colors[as.character(df_long$node)]
  
  pathway_order <- rev(kegg_top$PathwayName)
  gene_nodes <- unique(df_long$node[df_long$x == "GeneDisplay"])
  
  df_long <- df_long %>%
    mutate(
      node = factor(node, levels = c(gene_nodes, pathway_order)),
      next_node = factor(next_node, levels = levels(node))
    )
  
  p <- ggplot(
    df_long,
    aes(x = x, next_x = next_x, node = node, next_node = next_node,
        label = node, fill = fill_color)
  ) +
    geom_sankey(flow.alpha = 0.3, smooth = 8, width = 0.08, color = NA) +
    geom_sankey_text(aes(hjust = ifelse(x == "GeneDisplay", 1, 0)),
                     vjust = 0.5, size = 2.5, color = "black") +
    geom_point(data = sankey_data, aes(x = 1, y = 1, color = logFC),
               alpha = 0, inherit.aes = FALSE) +
    scale_fill_identity() +
    scale_color_gradientn(colors = cols_source, limits = limit_range,
                          name = "log2 Fold Change") +
    labs(title = title) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  
  return(p)
}

# ============================================================
# 5. 气泡图函数（通路顺序与桑基图严格一致）
# ============================================================
plot_kegg_dotplot <- function(ora_csv, title, top_n = 10) {
  
  kegg <- read.csv(ora_csv, stringsAsFactors = FALSE)
  kegg$PathwayName <- sapply(kegg$ID, get_kegg_name)
  
  kegg$exclude <- sapply(kegg$PathwayName, should_exclude)
  kegg_clean <- kegg[!kegg$exclude, ]
  
  n_excluded <- sum(kegg$exclude)
  cat(">>> [气泡图] 排除", n_excluded, "条无关通路，剩余", nrow(kegg_clean), "条\n")
  
  kegg_clean <- kegg_clean[order(kegg_clean$p.adjust), ]
  kegg_top <- head(kegg_clean, top_n)
  
  kegg_top$GeneRatio_num <- sapply(kegg_top$GeneRatio, function(x) {
    nums <- as.numeric(unlist(strsplit(x, "/")))
    nums[1] / nums[2]
  })
  
  kegg_top$PathwayName <- factor(kegg_top$PathwayName,
                                 levels = rev(kegg_top$PathwayName))
  
  p <- ggplot(kegg_top, aes(x = GeneRatio_num, y = PathwayName, fill = p.adjust)) +
    geom_point(aes(size = Count), shape = 21, color = "black", stroke = 0.3) +
    scale_size(range = c(3, 10), name = "Gene Count") +
    scale_fill_gradientn(
      colors = c("#FFCD00", "#eeeeee", "#3070B3"),
      values = scales::rescale(c(0, 0.5, 1)),
      name = "FDR"
    ) +
    labs(x = "Gene Ratio", y = NULL, title = title) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 12)
    )
  
  return(p)
}

# ============================================================
# 6. 运行：三个比较组，桑基图和气泡图分别独立预览
# ============================================================

cat("\n========== T32 vs Ctrl ==========\n")
p_sankey_t32 <- plot_kegg_sankey(
  ora_csv = "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T32vsCtrl.csv",
  deg_file = "results/RNAseq_analysis/Results/DESeq2_T32_vs_Ctrl.csv",
  title = "KEGG Sankey — T32 vs Ctrl",
  top_n = 10,
  max_genes_per_pathway = 10
)
p_dot_t32 <- plot_kegg_dotplot(
  ora_csv = "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T32vsCtrl.csv",
  title = "KEGG Dotplot — T32 vs Ctrl",
  top_n = 10
)
print(p_sankey_t32)
print(p_dot_t32)

cat("\n========== T24 vs Ctrl ==========\n")
p_sankey_t24 <- plot_kegg_sankey(
  ora_csv = "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T24vsCtrl.csv",
  deg_file = "results/RNAseq_analysis/Results/DESeq2_T24_vs_Ctrl.csv",
  title = "KEGG Sankey — T24 vs Ctrl",
  top_n = 10,
  max_genes_per_pathway = 10
)
p_dot_t24 <- plot_kegg_dotplot(
  ora_csv = "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T24vsCtrl.csv",
  title = "KEGG Dotplot — T24 vs Ctrl",
  top_n = 10
)
print(p_sankey_t24)
print(p_dot_t24)

cat("\n========== T32 vs T24 ==========\n")
p_sankey_t32t24 <- plot_kegg_sankey(
  ora_csv = "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T32vsT24.csv",
  deg_file = "results/RNAseq_analysis/Results/DESeq2_T32_vs_T24.csv",
  title = "KEGG Sankey — T32 vs T24",
  top_n = 10,
  max_genes_per_pathway = 10
)
p_dot_t32t24 <- plot_kegg_dotplot(
  ora_csv = "results/RNAseq_analysis/Results/KEGG_ORA_Pathway_T32vsT24.csv",
  title = "KEGG Dotplot — T32 vs T24",
  top_n = 10
)
print(p_sankey_t32t24)
print(p_dot_t32t24)

cat("\n>>> 全部预览完成\n")

library(ggplot2)
library(dplyr)
library(GO.db)
library(stringr)
library(gground)
library(ggprism)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ---------- 必须重新定义的函数和变量 ----------
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(y) as.numeric(y[1])/as.numeric(y[2]))
}

get_go_info <- function(go_ids) {
  go_ids <- unique(as.character(go_ids))
  go_info <- AnnotationDbi::select(GO.db, keys = go_ids, 
                                   columns = c("TERM", "ONTOLOGY"), 
                                   keytype = "GOID")
  name_map <- setNames(go_info$TERM, go_info$GOID)
  ont_map  <- setNames(go_info$ONTOLOGY, go_info$GOID)
  list(names = name_map, ontology = ont_map)
}

pal <- c("BP" = '#7bc4e2', "CC" = '#acd372', "MF" = '#fbb05b')

exclude_terms <- c(
  "negative regulation of T cell apoptotic process",
  "Z disc", "I band", "striated muscle dense body",
  "melanosome", "male germ cell nucleus",
  "tumor necrosis factor receptor activity",
  "death receptor activity"
)

# 基因ID → 名称映射
anno <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors = FALSE)
id2name <- setNames(anno$Preferred_name, anno$gene_id)
id2name[id2name == "-" | is.na(id2name)] <- names(id2name)[id2name == "-" | is.na(id2name)]

# ========== T32vsT24：数据准备 ==========
go_t32t24 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T32vsT24.csv", stringsAsFactors = FALSE)
go_t32t24 <- go_t32t24 %>% filter(!grepl("^obsolete", Description, ignore.case = TRUE))

go_info_t32t24 <- get_go_info(go_t32t24$ID)
go_t32t24$DisplayName <- go_info_t32t24$names[go_t32t24$ID]
go_t32t24$DisplayName[is.na(go_t32t24$DisplayName)] <- go_t32t24$Description[is.na(go_t32t24$DisplayName)]
go_t32t24$ONTOLOGY <- go_info_t32t24$ontology[go_t32t24$ID]

go_t32t24 <- go_t32t24 %>% 
  filter(!DisplayName %in% exclude_terms) %>%
  filter(!grepl("T cell|striated muscle|Z disc|I band|melanosome|male germ cell", DisplayName, ignore.case = TRUE))

go_t32t24_top <- go_t32t24 %>%
  filter(p.adjust < 0.05, !is.na(ONTOLOGY)) %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(ONTOLOGY, p.adjust)

if (!"geneID" %in% names(go_t32t24_top)) {
  raw_t32t24 <- read.csv("results/RNAseq_analysis/Results/GO_ORA_T32vsT24.csv", stringsAsFactors = FALSE)
  go_t32t24_top$geneID <- raw_t32t24$geneID[match(go_t32t24_top$ID, raw_t32t24$ID)]
}

# ========== T32vsT24：条形图数据构造 ==========
go_t32t24_top <- go_t32t24_top %>%
  mutate(
    index = row_number(),
    logp  = -log10(p.adjust),
    DisplayName = factor(DisplayName, levels = DisplayName),
    geneID_show = sapply(strsplit(as.character(geneID), "/"), function(x) {
      names <- id2name[x]
      names[is.na(names)] <- x[is.na(names)]
      if (length(names) > 3) paste(paste(names[1:3], collapse = "/"), "...", sep = "")
      else paste(names, collapse = "/")
    })
  )

width_t32t24 <- 0.5
xaxis_max_t32t24 <- 8

rect_t32t24 <- go_t32t24_top %>%
  group_by(ONTOLOGY) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    xmin = -3 * width_t32t24,
    xmax = -2 * width_t32t24,
    ymax = cumsum(n),
    ymin = lag(ymax, default = 0) + 0.6,
    ymax = ymax + 0.4
  )

# ========== T32vsT24：画图 ==========
p3c <- ggplot(go_t32t24_top, aes(x = logp, y = index, fill = ONTOLOGY)) +
  geom_round_col(aes(y = DisplayName), width = 0.6, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = DisplayName), hjust = 0, size = 3.5) +
  geom_text(aes(x = 0.1, label = geneID_show, colour = ONTOLOGY),
            hjust = 0, vjust = 2.6, size = 2.2, fontface = 'italic', show.legend = FALSE) +
  geom_point(aes(x = -width_t32t24, size = Count), shape = 21) +
  geom_text(aes(x = -width_t32t24, label = Count), size = 2.5) +
  scale_size_continuous(name = 'Count', range = c(3, 8)) +
  geom_round_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
                  data = rect_t32t24, radius = unit(2, 'mm'), inherit.aes = FALSE) +
  geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
            data = rect_t32t24, inherit.aes = FALSE, angle = 90, size = 3, fontface = "bold") +
  geom_segment(aes(x = 0, y = 0, xend = xaxis_max_t32t24, yend = 0),
               linewidth = 1.5, inherit.aes = FALSE) +
  labs(y = NULL, title = "C  T32 vs T24") +
  scale_fill_manual(name = 'Category', values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, ceiling(xaxis_max_t32t24), 2), expand = expansion(c(0, 0))) +
  theme_prism() +
  theme(axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5))

print(p3c)