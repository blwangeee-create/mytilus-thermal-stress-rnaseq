# ============================================================
# 三向交集验证(修正版)
# 用新的 hub 表(power7 / power18)重跑交集,验证候选基因名单
# 原交集逻辑见 wgcna.R 第 737-853 行(symbol 版)
# ============================================================

suppressPackageStartupMessages({ library(WGCNA) })

BASE <- "/Volumes/Expansion/NCBI/PRJNA934294"
CORR <- file.path(BASE, "results/RNAseq_analysis/Results/WGCNA_corrected")

# ---------- 1. ID 映射(MCOR ↔ symbol), 与原脚本同法 ----------
header_line <- readLines(file.path(BASE, "annotation/GCA_011752425.2_MCOR1.1_feature_table.txt"), n = 1)
col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]
feat <- read.delim(file.path(BASE, "annotation/GCA_011752425.2_MCOR1.1_feature_table.txt"),
                   stringsAsFactors = FALSE, skip = 1, header = FALSE, col.names = col_names)

mrna_map <- unique(feat[feat$feature == "mRNA" & feat$related_accession != "",
                        c("locus_tag", "related_accession")])
colnames(mrna_map) <- c("locus_tag", "protein_id")

eggnog <- read.delim(file.path(BASE, "annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular"),
                     stringsAsFactors = FALSE, comment.char = "")
names(eggnog) <- gsub("^#", "", names(eggnog))
names(eggnog) <- gsub("^X\\.", "", names(eggnog))

symbol_map <- unique(eggnog[eggnog$Preferred_name != "" & eggnog$Preferred_name != "-",
                            c("query", "Preferred_name")])
colnames(symbol_map) <- c("protein_id", "symbol")

id_map <- merge(mrna_map, symbol_map, by = "protein_id")
id_map$gene_id <- paste0("gene-", id_map$locus_tag)
id_map <- unique(id_map[, c("gene_id", "symbol")])
cat("ID 映射表:", nrow(id_map), "条\n")

# ---------- 2. 交集函数 ----------
run_intersect <- function(power_tag) {
  hub_dir <- file.path(CORR, power_tag)
  hub_files <- list.files(hub_dir, pattern = "^hub_.*\\.csv$", full.names = TRUE)

  hub_sym <- unique(unlist(lapply(hub_files, function(f) {
    ids <- read.csv(f, stringsAsFactors = FALSE)$gene
    s <- id_map$symbol[id_map$gene_id %in% ids]
    s[!is.na(s) & s != ""]
  })))

  # PPI MCC top20(T32 与 T24, 用 name 列)
  ppi_t32 <- read.csv(file.path(BASE, "results/PPI/32/T32_string_edges.tsv_MCC_top20 default node.csv"),
                      stringsAsFactors = FALSE)$name[1:20]
  ppi_t24 <- read.csv(file.path(BASE, "results/PPI/24/T24_string_edges.tsv_MCC_top20 default node.csv"),
                      stringsAsFactors = FALSE)$name[1:20]

  # k-means 簇(转 symbol)
  k3 <- read.csv(file.path(BASE, "results/RNAseq_analysis/Results/kmeans/K3_cluster_3.csv"),
                 stringsAsFactors = FALSE)$gene_id
  k1 <- read.csv(file.path(BASE, "results/RNAseq_analysis/Results/kmeans/K3_cluster_1.csv"),
                 stringsAsFactors = FALSE)$gene_id
  k3_sym <- id_map$symbol[id_map$gene_id %in% k3]; k3_sym <- k3_sym[!is.na(k3_sym) & k3_sym != ""]
  k1_sym <- id_map$symbol[id_map$gene_id %in% k1]; k1_sym <- k1_sym[!is.na(k1_sym) & k1_sym != ""]

  core_heat <- Reduce(intersect, list(hub_sym, ppi_t32, k3_sym))
  core_t24  <- Reduce(intersect, list(hub_sym, ppi_t24, k1_sym))

  list(hub_n = length(hub_sym), core_heat = core_heat, core_t24 = core_t24)
}

cat("\n=== 三向交集验证 (Hub ∩ PPI-T32 ∩ k-means Cluster3) ===\n")
r7  <- run_intersect("power7")
r18 <- run_intersect("power18")

cat(sprintf("\npower=7  : hub 基因数 %d -> 交集 %d 个: %s\n",
            r7$hub_n, length(r7$core_heat), paste(r7$core_heat, collapse = ", ")))
cat(sprintf("power=18 : hub 基因数 %d -> 交集 %d 个: %s\n",
            r18$hub_n, length(r18$core_heat), paste(r18$core_heat, collapse = ", ")))

# 原结果对照
old_core <- read.csv(file.path(BASE, "results/RNAseq_analysis/Results/core_candidates_heat_3way.csv"),
                     stringsAsFactors = FALSE)$gene
cat(sprintf("\n原结果(未修正) : %s\n", paste(old_core, collapse = ", ")))

# 交集稳定性
cat("\n=== 稳定性判定 ===\n")
cat("原结果三基因:", paste(sort(old_core), collapse = ", "), "\n")
cat("power7  三基因:", paste(sort(r7$core_heat), collapse = ", "), "\n")
cat("power18 三基因:", paste(sort(r18$core_heat), collapse = ", "), "\n")
cat("power7 与原结果一致:", identical(sort(old_core), sort(r7$core_heat)), "\n")
cat("power18 与原结果一致:", identical(sort(old_core), sort(r18$core_heat)), "\n")

# ---------- 3. 保存验证结果 ----------
ver <- data.frame(
  power = c("original", "power7", "power18"),
  n_hub_genes = c(NA, r7$hub_n, r18$hub_n),
  core_heat_n = c(length(old_core), length(r7$core_heat), length(r18$core_heat)),
  core_heat_genes = c(paste(old_core, collapse = ","),
                      paste(r7$core_heat, collapse = ","),
                      paste(r18$core_heat, collapse = ",")),
  stringsAsFactors = FALSE)
write.csv(ver, file.path(CORR, "core_candidates_verification.csv"), row.names = FALSE)
cat("\n已保存:", file.path(CORR, "core_candidates_verification.csv"), "\n")
