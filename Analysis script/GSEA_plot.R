# ============================================
# GSEA Analysis - T24 vs Ctrl (修复版，含 ID 映射)
# ============================================
if (!require("ggridges", quietly = TRUE)) install.packages("ggridges")
library(ggridges)
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(ggplot2)
})

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294"

# ============================================
# 修改这两个变量，其他不变
# ============================================

comparison     <- "T32_vs_T24"
res_csv        <- "DESeq2_T32_vs_T24.csv"
pvalue_cutoff  <- 0.05
padj_method    <- "BH"

# ------------------ 1. 读取 eggNOG 注释 ------------------
cat("========== 1. 读取 eggNOG 注释 ==========\n")

anno_dir <- file.path(base_dir, "annotation")
anno_files <- list.files(anno_dir, pattern = "eggNOG.*annotations.*tabular|Galaxy4.*eggNOG.*tabular", full.names = TRUE)
if (length(anno_files) == 0) stop("未找到 eggNOG 注释文件")

eggnog <- read.delim(anno_files[1], header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE, quote = "", comment.char = "", check.names = FALSE)
cat("注释维度:", nrow(eggnog), "x", ncol(eggnog), "\n")

# ------------------ 1.5 关键修复：ID 映射 ------------------
cat("\n========== 1.5 建立 ID 映射 (CAC → gene-MCOR) ==========\n")

ft <- read.delim(file.path(base_dir, "annotation/GCA_011752425.2_MCOR1.1_feature_table.txt"),
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# 取 product_accession 和 locus_tag，去重（一个基因多行 feature）
id_map <- ft[, c("product_accession", "locus_tag")]
id_map <- id_map[id_map$product_accession != "" & !is.na(id_map$product_accession) &
                   id_map$locus_tag != "" & !is.na(id_map$locus_tag), ]
id_map <- id_map[!duplicated(id_map$product_accession), ]  # 一个 CAC 只对应一个 MCOR

# 加 gene- 前缀，匹配 DESeq2 行名
id_map$gene_id <- paste0("gene-", id_map$locus_tag)

cat("映射对数:", nrow(id_map), "\n")
cat("示例:\n")
print(head(id_map, 5))

# 映射：eggNOG 第一列 (CAC...) → gene-MCOR...
# 用列索引 [1] 代替名字，彻底避开命名问题
m_idx <- match(eggnog[, 1], id_map[["product_accession"]])
eggnog[["gene_id"]] <- id_map[["gene_id"]][m_idx]

# 验证
mapped_n <- sum(!is.na(eggnog[["gene_id"]]))
cat("成功映射:", mapped_n, "/", nrow(eggnog), "\n\n")

# ------------------ 2. 构建 KEGG TERM2GENE ------------------
cat("========== 2. 构建 KEGG Pathway 基因集 ==========\n")

kegg_df <- data.frame(term = character(), gene = character(), stringsAsFactors = FALSE)
if ("KEGG_Pathway" %in% colnames(eggnog)) {
  kegg_list <- strsplit(eggnog$KEGG_Pathway, ",")
  kegg_df <- data.frame(
    term = unlist(kegg_list),
    gene = rep(eggnog$gene_id, sapply(kegg_list, length)),
    stringsAsFactors = FALSE
  )
  kegg_df <- kegg_df[kegg_df$term != "" & kegg_df$term != "-" & !is.na(kegg_df$term) & !is.na(kegg_df$gene), ]
  kegg_df <- unique(kegg_df)
  # 只保留 map 前缀的通路，删除 ko 前缀的
  kegg_df <- kegg_df[grepl("^map", kegg_df$term), ]
  cat("KEGG 通路-基因对:", nrow(kegg_df), "| 唯一通路:", length(unique(kegg_df$term)), "\n")
} else {
  cat("警告: 无 KEGG_Pathway 列\n")
}

# ------------------ 3. 构建 GO TERM2GENE + TERM2NAME ------------------
cat("\n========== 3. 构建 GO 基因集 ==========\n")

go_df <- data.frame(term = character(), gene = character(), stringsAsFactors = FALSE)
go_name_map <- data.frame(term = character(), name = character(), stringsAsFactors = FALSE)

if ("GOs" %in% colnames(eggnog)) {
  go_list <- strsplit(eggnog$GOs, ",")
  go_df <- data.frame(
    term = unlist(go_list),
    gene = rep(eggnog$gene_id, sapply(go_list, length)),
    stringsAsFactors = FALSE
  )
  go_df <- go_df[go_df$term != "" & go_df$term != "-" & !is.na(go_df$term) & !is.na(go_df$gene), ]
  go_df <- unique(go_df)
  cat("GO term-基因对:", nrow(go_df), "| 唯一 GO:", length(unique(go_df$term)), "\n")
  
  # ===== 新增：构建 GO ID → 名称映射表 =====
  unique_go_ids <- unique(go_df$term)
  if (length(unique_go_ids) > 0) {
    if (!require("GO.db", quietly = TRUE)) BiocManager::install("GO.db")
    library(GO.db)
    
    go_desc <- AnnotationDbi::select(GO.db, keys = unique_go_ids,
                                     columns = c("GOID", "TERM"),
                                     keytype = "GOID")
    go_desc <- go_desc[!duplicated(go_desc$GOID), ]
    
    go_name_map <- data.frame(
      term = go_desc$GOID,
      name = go_desc$TERM,
      stringsAsFactors = FALSE
    )
    # GO.db 中没有的 ID，用 ID 本身作为名称
    missing_ids <- setdiff(unique_go_ids, go_name_map$term)
    if (length(missing_ids) > 0) {
      go_name_map <- rbind(go_name_map, 
                           data.frame(term = missing_ids, name = missing_ids, stringsAsFactors = FALSE))
    }
    cat("GO 名称映射表:", nrow(go_name_map), "条\n")
  }
} else {
  cat("警告: 无 GOs 列\n")
}

# ------------------ 4. 读取 DESeq2 完整结果 ------------------
cat("\n========== 4. 读取 DESeq2 结果:", comparison, "==========\n")

res_path <- file.path(base_dir, "results/RNAseq_analysis/Results", res_csv)
res <- read.csv(res_path, row.names = 1, stringsAsFactors = FALSE)
cat("结果维度:", nrow(res), "x", ncol(res), "\n")

if (!"stat" %in% colnames(res)) stop("缺少 'stat' 列")

# ------------------ 5. 构建排序基因列表 ------------------
cat("\n========== 5. 构建 GSEA 输入 (按 stat 排序) ==========\n")

res_clean <- res[!is.na(res$stat), ]
geneList <- res_clean$stat
names(geneList) <- rownames(res_clean)
geneList <- sort(geneList, decreasing = TRUE)

cat("输入基因总数:", length(geneList), "\n")
cat("Top 5:\n"); print(head(geneList, 5))

# ------------------ 6. KEGG GSEA ------------------
cat("\n========== 6. KEGG GSEA ==========\n")

gse_kegg <- NULL
if (nrow(kegg_df) > 0) {
  n_annot <- sum(names(geneList) %in% kegg_df$gene)
  cat("有 KEGG 注释的基因:", n_annot, "/", length(geneList), "\n")
  
  if (n_annot > 0) {
    gse_kegg <- GSEA(
      geneList      = geneList,
      TERM2GENE     = kegg_df,
      pvalueCutoff  = pvalue_cutoff,
      pAdjustMethod = padj_method,
      minGSSize     = 10,
      maxGSSize     = 500,
      verbose       = FALSE,
      seed          = 123
    )
    
    cat("显著 KEGG 通路 (padj<0.05):", nrow(gse_kegg), "\n")
    if (nrow(gse_kegg) > 0) {
      kegg_top <- as.data.frame(gse_kegg)[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")]
      cat("\nTop 10:\n")
      print(kegg_top[1:min(10, nrow(kegg_top)), ])
    }
  }
}

# ------------------ 7. GO GSEA（不带 TERM2NAME）------------------
cat("\n========== 7. GO GSEA ==========\n")

gse_go <- NULL
if (nrow(go_df) > 0) {
  n_annot_go <- sum(names(geneList) %in% go_df$gene)
  cat("有 GO 注释的基因:", n_annot_go, "/", length(geneList), "\n")
  
  if (n_annot_go > 0) {
    # 关键：不传入 TERM2NAME，避免 clusterProfiler bug
    gse_go <- GSEA(
      geneList      = geneList,
      TERM2GENE     = go_df,
      pvalueCutoff  = pvalue_cutoff,
      pAdjustMethod = padj_method,
      minGSSize     = 10,
      maxGSSize     = 500,
      verbose       = FALSE,
      seed          = 123
    )
    
    cat("显著 GO term (padj<0.05):", nrow(gse_go), "\n")
    if (nrow(gse_go) > 0) {
      go_top <- as.data.frame(gse_go)[, c("ID", "Description", "setSize", "NES", "pvalue", "p.adjust")]
      cat("\nTop 10 (转换前):\n")
      print(go_top[1:min(10, nrow(go_top)), ])
    }
  }
}

# ------------------ 7.5 转换通路/GO 名称 ------------------
# ------------------ 7.5 转换通路/GO 名称（最终版）------------------
cat("\n========== 7.5 转换通路/GO 名称 ==========\n")

# --- KEGG 名称转换 ---
if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {
  # 方法1：尝试用 enrichplot 的 setReadable（如果有 OrgDb）
  # 方法2：直接用 keggGet 逐个查询，但增加超时和错误处理
  # 方法3：用本地缓存文件
  
  # 这里用最简单的方法：直接从 KEGG API 下载 map 通路列表
  # KEGG 提供了一个干净的参考通路列表
  tryCatch({
    # 下载参考通路名称（只包含 map 前缀的通用名称）
    kegg_pathway_url <- "https://rest.kegg.jp/list/pathway"
    kegg_pathway_txt <- readLines(kegg_pathway_url, warn = FALSE)
    
    # 解析：每行格式为 "path:map00010\tGlycolysis / Gluconeogenesis"
    kegg_pathway_df <- do.call(rbind, lapply(kegg_pathway_txt, function(line) {
      parts <- strsplit(line, "\t")[[1]]
      if (length(parts) >= 2) {
        id <- gsub("path:", "", parts[1])
        name <- parts[2]
        # 去掉物种后缀
        name <- gsub(" - .*$", "", name)
        data.frame(ID = id, name = name, stringsAsFactors = FALSE)
      } else NULL
    }))
    
    # 匹配
    m_idx <- match(gse_kegg$ID, kegg_pathway_df$ID)
    matched_names <- kegg_pathway_df$name[m_idx]
    
    # 应用转换
    gse_kegg@result$Description <- ifelse(is.na(matched_names), gse_kegg$ID, matched_names)
    
    cat("KEGG 名称转换完成:", sum(!is.na(matched_names)), "/", nrow(gse_kegg), "\n")
  }, error = function(e) {
    cat("KEGG 在线查询失败，使用原始 ID:", conditionMessage(e), "\n")
  })
}
# 在 KEGG 名称转换的最后，增加强制替换
kegg_forced_fix <- c(
  "Human T-cell leukemia virus 1 infection" = "HTLV-I infection",
  "Human cytomegalovirus infection" = "Cytomegalovirus infection",
  "Human papillomavirus infection" = "Papillomavirus infection",
  "Hepatitis B" = "Hepatitis B",
  "Hepatitis C" = "Hepatitis C",
  "Influenza A" = "Influenza A",
  "Epstein-Barr virus infection" = "EBV infection",
  "Kaposi sarcoma-associated herpesvirus infection" = "KSHV infection",
  "Measles" = "Measles",
  "Amoebiasis" = "Amoebiasis",
  "Malaria" = "Malaria",
  "Leishmaniasis" = "Leishmaniasis",
  "Chagas disease" = "Chagas disease",
  "Toxoplasmosis" = "Toxoplasmosis",
  "Small cell lung cancer" = "Small cell lung cancer",
  "Protein processing in endoplasmic reticulum" = "Protein processing in ER",
  "Antigen processing and presentation" = "Antigen processing and presentation",
  "Longevity regulating pathway" = "Longevity regulating pathway",
  "TNF signaling pathway" = "TNF signaling pathway",
  "Insulin resistance" = "Insulin resistance",
  "MAPK signaling pathway" = "MAPK signaling pathway",
  "Apoptosis" = "Apoptosis",
  "Toll and Imd signaling pathway" = "Toll and Imd signaling pathway",
  "NF-kappa B signaling pathway" = "NF-kappa B signaling pathway",
  "Cell cycle" = "Cell cycle",
  "Ribosome" = "Ribosome",
  "DNA replication" = "DNA replication",
  "Proteasome" = "Proteasome",
  "Ribosome biogenesis in eukaryotes" = "Ribosome biogenesis"
)

for (old_name in names(kegg_forced_fix)) {
  idx <- which(gse_kegg@result$Description == old_name)
  if (length(idx) > 0) {
    gse_kegg@result$Description[idx] <- kegg_forced_fix[old_name]
  }
}
# --- GO 名称转换 ---
if (!is.null(gse_go) && nrow(gse_go) > 0) {
  # 方法：用 GO.db + 在线 bioontology API 作为兜底
  go_ids <- gse_go$ID
  
  # 第一层：GO.db
  if (!require("GO.db", quietly = TRUE)) BiocManager::install("GO.db")
  library(GO.db)
  
  go_desc <- AnnotationDbi::select(GO.db, keys = go_ids,
                                   columns = c("GOID", "TERM"),
                                   keytype = "GOID")
  go_desc <- go_desc[!duplicated(go_desc$GOID), ]
  
  m_idx <- match(go_ids, go_desc$GOID)
  matched_go <- go_desc$TERM[m_idx]
  
  # 第二层：对于 GO.db 中没有的，用在线 API 查询
  still_na <- is.na(matched_go)
  if (any(still_na)) {
    cat("尝试在线查询", sum(still_na), "个缺失的 GO 名称...\n")
    
    missing_ids <- go_ids[still_na]
    
    # 使用 bioontology API（不需要 API key，免费）
    for (go_id in missing_ids) {
      tryCatch({
        url <- paste0("http://www.ebi.ac.uk/ols/api/ontologies/go/terms?iri=http://purl.obolibrary.org/obo/", gsub(":", "_", go_id))
        response <- httr::GET(url, httr::timeout(5))
        if (httr::status_code(response) == 200) {
          content <- httr::content(response, "parsed")
          if (!is.null(content$`_embedded`$terms[[1]]$label)) {
            name <- content$`_embedded`$terms[[1]]$label
            matched_go[go_ids == go_id] <- name
          }
        }
      }, error = function(e) {})
    }
  }
  
  # 应用转换
  gse_go@result$Description <- ifelse(is.na(matched_go), go_ids, matched_go)
  
  cat("GO 名称转换完成:", sum(!is.na(matched_go)), "/", length(go_ids), "\n")
  
  # 检查剩余未匹配的
  remaining <- sum(grepl("^GO:", gse_go@result$Description))
  if (remaining > 0) {
    cat("警告: 仍有", remaining, "个 GO term 无法解析\n")
  }
}
# ------------------ 8. 可视化（仅预览）------------------
cat("\n========== 8. 预览图 ==========\n")

if (!require("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")
library(enrichplot)

# 公共颜色定义（KEGG 和 GO 共用：Up = 金黄 #FFCD00, Down = 蓝 #0B8BEE）
go_cols <- c("Up" = "#FFCD00", "Down" = "#0B8BEE")

if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {

  # ============================================
  # KEGG GSEA 曲线图 - 修复表格颜色
  # ============================================
  
  library(ggplot2)
  library(patchwork)
  
  # 定义通路
  # ============================================
  # 固定三条目标通路（按名称查找，适用于所有对比组）
  # ============================================
  
  target_names <- c(
    "Protein processing in ER",
    "Apoptosis", 
    "Ribosome biogenesis"
  )
  
  # 按名称查找 geneSetID
  find_geneSetID <- function(gsea_obj, target_name) {
    exact <- which(gsea_obj@result$Description == target_name)
    if (length(exact) > 0) return(exact[1])
    fuzzy <- grep(target_name, gsea_obj@result$Description, ignore.case = TRUE)
    if (length(fuzzy) > 0) return(fuzzy[1])
    warning("未找到通路: ", target_name)
    return(NA)
  }
  
  my_geneSetID <- sapply(target_names, function(name) find_geneSetID(gse_kegg, name))
  my_geneSetID <- my_geneSetID[!is.na(my_geneSetID)]
  
  # 保存到全局环境
  assign("my_geneSetID", my_geneSetID, envir = .GlobalEnv)
  
  cat("固定三条通路:\n")
  for (i in seq_along(my_geneSetID)) {
    cat("  [", my_geneSetID[i], "] ", gse_kegg@result$Description[my_geneSetID[i]], 
        " (NES = ", round(gse_kegg@result$NES[my_geneSetID[i]], 3), ")\n", sep = "")
  }
  pathway_info <- data.frame(
    geneSetID = my_geneSetID,
    Description = gse_kegg@result$Description[my_geneSetID],
    NES = gse_kegg@result$NES[my_geneSetID],
    padj = gse_kegg@result$p.adjust[my_geneSetID],
    setSize = gse_kegg@result$setSize[my_geneSetID],
    stringsAsFactors = FALSE
  )
  
  # 动态颜色映射：只给找到的通路分配颜色
  color_map <- c(
    "Protein processing in ER" = "#FFCD00",
    "Apoptosis" = "#C74375",
    "Ribosome biogenesis" = "#2980B9"
  )
  pathway_info$peak_color <- color_map[pathway_info$Description]
  
  titanium_gray <- "#A0A0A0"
  
  format_p <- function(p) {
    sapply(p, function(x) {
      if (x < 0.001) sprintf("%.2e", x)
      else if (x < 0.01) sprintf("%.3f", x)
      else sprintf("%.4f", x)
    })
  }
  
  # ============================================
  # 创建渐变色线段
  # ============================================
  
  create_gradient_segments <- function(dd, peak_color, n_segments = 100) {
    
    peak_idx <- which.max(abs(dd$runningScore))
    peak_x <- dd$x[peak_idx]
    
    max_dist <- max(abs(dd$x - peak_x))
    dist_ratio <- abs(dd$x - peak_x) / max_dist
    
    col_ramp <- colorRampPalette(c(peak_color, titanium_gray))(n_segments + 1)
    
    color_idx <- round(dist_ratio * n_segments) + 1
    color_idx <- pmin(pmax(color_idx, 1), n_segments + 1)
    
    segments <- data.frame(
      x = head(dd$x, -1),
      xend = tail(dd$x, -1),
      y = head(dd$runningScore, -1),
      yend = tail(dd$runningScore, -1),
      color = head(col_ramp[color_idx], -1),
      stringsAsFactors = FALSE
    )
    
    return(segments)
  }
  
  # ============================================
  # 提取所有通路数据
  # ============================================
  
  all_segments <- list()
  all_barcodes <- list()
  all_ridges <- list()
  
  for (i in 1:nrow(pathway_info)) {
    
    gs_id <- pathway_info$geneSetID[i]
    dd <- enrichplot:::gsInfo(gse_kegg, geneSetID = gs_id)
    
    seg <- create_gradient_segments(dd, pathway_info$peak_color[i], n_segments = 80)
    seg$pathway <- pathway_info$Description[i]
    all_segments[[i]] <- seg
    
    bar <- data.frame(x = dd$x[dd$position == 1], pathway = pathway_info$Description[i])
    all_barcodes[[i]] <- bar
    
    ridge <- dd[, c("x", "runningScore", "geneList")]
    ridge$pathway <- pathway_info$Description[i]
    all_ridges[[i]] <- ridge
  }
  
  seg_df <- do.call(rbind, all_segments)
  bar_df <- do.call(rbind, all_barcodes)
  ridge_df <- do.call(rbind, all_ridges)
  
  # ============================================
  # 上部分：三条曲线 + 图例表格
  # ============================================
  
  # 先画曲线
  p_curves_base <- ggplot(seg_df) +
    geom_segment(
      aes(x = x, xend = xend, y = y, yend = yend),
      color = seg_df$color,
      linewidth = 1.3, lineend = "round"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    labs(
      title = paste("KEGG GSEA -", comparison),
      x = NULL,
      y = "Enrichment Score"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.y = element_text(size = 11),
      axis.text = element_text(size = 10),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(seg_df$x)))
  
  # 添加彩色图例（用 geom_point 模拟）
  # 添加彩色图例（用 geom_point 模拟）
  legend_data <- data.frame(
    x = rep(NA, nrow(pathway_info)),
    y = rep(NA, nrow(pathway_info)),
    pathway = pathway_info$Description,
    color = pathway_info$peak_color
  )
  
  p_curves <- p_curves_base +
    geom_point(
      data = legend_data,
      aes(x = x, y = y, color = pathway),
      size = 3, show.legend = TRUE
    ) +
    scale_color_manual(
      name = NULL,
      values = setNames(pathway_info$peak_color, pathway_info$Description)
    ) +
    theme(
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.text = element_text(size = 10),
      legend.key = element_blank()
    )
  
  # ============================================
  # 中部分：三个条形码
  # ============================================
  
  bar_df$y_offset <- match(bar_df$pathway, pathway_info$Description)
  bar_height <- 0.35
  bar_df$y_bottom <- bar_df$y_offset - bar_height
  bar_df$y_top <- bar_df$y_offset + bar_height
  bar_df$bar_color <- pathway_info$peak_color[match(bar_df$pathway, pathway_info$Description)]
  
  p_barcodes <- ggplot(bar_df, aes(x = x)) +
    geom_segment(
      aes(xend = x, y = y_bottom, yend = y_top),
      color = bar_df$bar_color,
      linewidth = 0.5
    ) +
    scale_y_continuous(
      breaks = 1:nrow(pathway_info),
      labels = pathway_info$Description,
      expand = c(0, 0),
      limits = c(0.5, 3.5)
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(bar_df$x))) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 10, color = pathway_info$peak_color),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 5, 0, 5)
    )
  
  # ============================================
  # 下部分：三个山脊
  # ============================================
  
  ridge_df$pathway_num <- match(ridge_df$pathway, pathway_info$Description)
  ridge_colors <- sapply(pathway_info$peak_color, function(col) {
    colorRampPalette(c(col, "#F5F5F5"))(5)[3]
  })
  ridge_df$ridge_fill <- ridge_colors[match(ridge_df$pathway, pathway_info$Description)]
  
  p_ridges <- ggplot(ridge_df, aes(x = x, y = geneList)) +
    geom_area(
      aes(group = pathway),
      fill = ridge_df$ridge_fill,
      color = titanium_gray,
      alpha = 0.5,
      linewidth = 0.3
    ) +
    facet_grid(pathway ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(ridge_df$x))) +
    labs(x = "Gene Rank", y = "Ranked List Metric") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "white"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # ============================================
  # 图注表格（用 ggplot2 画，带颜色）
  # ============================================
  
  # 创建表格用的数据
  table_plot_data <- data.frame(
    y = nrow(pathway_info):1,
    label = pathway_info$Description,
    NES = sprintf("%.3f", pathway_info$NES),
    padj = format_p(pathway_info$padj),
    size = pathway_info$setSize,
    color = pathway_info$peak_color,
    stringsAsFactors = FALSE
  )
  
  table_plot <- ggplot(table_plot_data) +
    # 通路名称（彩色）
    geom_text(aes(x = 0.1, y = y, label = label, color = color), 
              hjust = 0, size = 3.5, fontface = "bold") +
    # NES
    geom_text(aes(x = 0.5, y = y, label = paste0("NES = ", NES)), 
              hjust = 0, size = 3.2, color = "black") +
    # padj
    geom_text(aes(x = 0.7, y = y, label = paste0("padj = ", padj)), 
              hjust = 0, size = 3.2, color = "black") +
    # size
    geom_text(aes(x = 0.9, y = y, label = paste0("n = ", size)), 
              hjust = 0, size = 3.2, color = "black") +
    scale_color_identity() +
    scale_x_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 3.5), expand = c(0, 0)) +
    labs(title = "Significant Pathways") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # ============================================
  # 组合
  # ============================================
  
  p_combined <- p_curves / p_barcodes / p_ridges / table_plot +
    plot_layout(heights = c(3.5, 0.8, 3, 1))
  
  print(p_combined)
  
  cat(">>> KEGG ridge plot...\n")
  library(ggridges)
  
  ridge_dat <- as.data.frame(gse_kegg)
  ridge_dat <- ridge_dat[order(ridge_dat$p.adjust), ]
  
  # 剔除指定通路
  exclude_pathways <- c("Chronic myeloid leukemia", "Collecting duct acid secretion", "Pertussis", "Viral myocarditis", "Osteoclast differentiation", "EBV infection", "Measles", "Small cell lung cancer", "HTLV-I infection", "Toxoplasmosis", "Th17 cell differentiation", "Legionellosis")
  ridge_dat <- ridge_dat[!ridge_dat$Description %in% exclude_pathways, ]
  
  # 去重：按 Description 保留 p.adjust 最显著的（即第一个出现的）
  ridge_dat <- ridge_dat[!duplicated(ridge_dat$Description), ]
  
  # 取 top 20
  ridge_dat <- head(ridge_dat, 20)
  
  # 按 NES 排序：负值在上，正值在下
  ridge_dat <- ridge_dat[order(ridge_dat$NES, decreasing = FALSE), ]
  
  # 构建每个通路的核心基因排名分布
  gene_ranks <- setNames(seq_along(geneList), names(geneList))
  ridge_long <- do.call(rbind, lapply(seq_len(nrow(ridge_dat)), function(i) {
    gs <- ridge_dat$ID[i]
    genes <- gse_kegg@geneSets[[gs]]
    core_genes <- genes[genes %in% names(gene_ranks)]
    ranks <- gene_ranks[core_genes]
    ranks <- ranks[!is.na(ranks)]
    data.frame(
      rank = ranks,
      Description = ridge_dat$Description[i],
      NES = ridge_dat$NES[i],
      stringsAsFactors = FALSE
    )
  }))
  
  ridge_long$Description <- factor(ridge_long$Description, levels = ridge_dat$Description)
  
  p_ridge <- ggplot(ridge_long, aes(x = rank, y = Description, fill = NES)) +
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#0B8BEE", mid = "#F5F5F5", high = "#FFCD00", midpoint = 0) +
    labs(x = "Gene Rank", y = NULL, title = paste("KEGG Ridge Plot -", comparison)) +
    theme_ridges() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(
        color = ifelse(ridge_dat$NES > 0, "#FFCD00", "#0B8BEE"),
        size = 12
      )
    )
  print(p_ridge)
  
  cat(">>> KEGG lollipop plot with gradient stems...\n")
  
  kegg_dat <- as.data.frame(gse_kegg)
  kegg_dat <- kegg_dat[order(kegg_dat$p.adjust), ]
  
  # 剔除指定通路
  exclude_pathways <- c("Chronic myeloid leukemia", "Collecting duct acid secretion", "Pertussis", "Viral myocarditis", "Osteoclast differentiation", "EBV infection", "Measles", "Small cell lung cancer", "HTLV-I infection", "Toxoplasmosis", "Th17 cell differentiation", "Legionellosis")
  kegg_dat <- kegg_dat[!kegg_dat$Description %in% exclude_pathways, ]
  
  # 去重：按 Description 保留 p.adjust 最显著的（即第一个出现的）
  kegg_dat <- kegg_dat[!duplicated(kegg_dat$Description), ]
  
  # 取 top 20
  kegg_dat <- head(kegg_dat, 20)
  
  # 按 NES 排序：负值在上，正值在下
  kegg_dat <- kegg_dat[order(kegg_dat$NES, decreasing = FALSE), ]
  kegg_dat$pathway <- factor(kegg_dat$Description, levels = kegg_dat$Description)
  kegg_dat$direction <- ifelse(kegg_dat$NES > 0, "Up", "Down")
  
  # 定义颜色
  up_color   <- "#FFCD00"   # 金黄
  down_color <- "#0B8BEE"   # 蓝色
  root_color <- "#E0E0E0"   # 棒棒根部颜色（浅灰）
  
  # ===== 关键修改：生成渐变棒棒数据 =====
  # 将每条棒棒拆分成 n 段，每段一个颜色，形成渐变
  
  n_segments <- 60  # 每根棒棒的分段数，越大越平滑
  
  # 为每条通路生成渐变色棒棒
  grad_segments <- do.call(rbind, lapply(seq_len(nrow(kegg_dat)), function(i) {
    row <- kegg_dat[i, ]
    nes <- row$NES
    
    # 确定端点颜色
    end_color <- ifelse(nes > 0, up_color, down_color)
    
    # 生成从根部到端点的颜色渐变
    col_ramp <- colorRampPalette(c(root_color, end_color))(n_segments)
    
    # 将 NES 范围分成 n 段
    x_seq <- seq(0, nes, length.out = n_segments + 1)
    
    data.frame(
      pathway = row$pathway,
      NES     = nes,
      direction = row$direction,
      x_start = head(x_seq, -1),
      x_end   = tail(x_seq, -1),
      color   = col_ramp,
      stringsAsFactors = FALSE
    )
  }))
  
  # 计算 y 轴颜色（与原始代码一致）
  loli_y_colors <- ifelse(kegg_dat$NES > 0, up_color, down_color)
  
  p3 <- ggplot() +
    # 渐变棒棒（用多个小段模拟）
    geom_segment(
      data = grad_segments,
      aes(x = x_start, xend = x_end, y = pathway, yend = pathway, color = color),
      linewidth = 1.2, lineend = "butt"
    ) +
    # 手动设置颜色（覆盖 aes 中的 color 映射）
    scale_color_identity() +
    # 圆圈（端点）
    geom_point(
      data = kegg_dat,
      aes(NES, pathway, fill = direction),
      shape = 21, size = 3.5, color = "white", stroke = 0.4
    ) +
    scale_fill_manual(values = c("Up" = up_color, "Down" = down_color)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      title = paste("KEGG GSEA -", comparison)
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_text(color = loli_y_colors, size = 12),
      axis.text.x = element_text(color = "black", size = 11),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      plot.margin = margin(4, 4, 4, 4)
    )
  print(p3)
}
if (!is.null(gse_go) && nrow(gse_go) > 0) {
  go_dat <- as.data.frame(gse_go)
  go_dat <- go_dat[!grepl("^obsolete", go_dat$Description, ignore.case = TRUE), ]   # ← 插入这一行
  go_dat <- go_dat[order(go_dat$p.adjust), ]
  go_dat <- head(go_dat, 15)
  go_dat$direction <- ifelse(go_dat$NES > 0, "Up", "Down")
  
  # ===== 缩放 NES 值以压缩条形长度（不改变文字布局）=====
  scale_factor <- 0.4  # 0.3~0.6 可调，越小条形越短
  go_dat$NES_scaled <- go_dat$NES * scale_factor
  
  cat(">>> GO GSEA 条形图...\n")
  p_go_bar <- ggplot(go_dat, aes(NES_scaled, reorder(Description, NES))) +
    geom_col(aes(fill = direction), width = 0.66, alpha = 0.85) +
    geom_text(
      aes(
        x     = ifelse(NES > 0, -0.2 * scale_factor, 0.2 * scale_factor),
        hjust = ifelse(NES > 0,  1,   0),
        label = Description
      ),
      size = 4.1,
      color = "black"
    ) +
    scale_fill_manual(values = go_cols) +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.02)),
      breaks = pretty(go_dat$NES_scaled, n = 5),
      labels = function(x) round(x / scale_factor, 2)  # 标签显示原始 NES 值
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      title = "GO GSEA"
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 11),
      axis.title.x = element_text(face = "bold", size = 12),
      # ===== 去除边框，只保留 x 轴 =====
      axis.line.x = element_line(color = "black", linewidth = 0.5),
      axis.line.y = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
  print(p_go_bar)
} else {
  cat(">>> GO: 无显著结果，跳过绘图\n")
}
# ------------------ 9. 论文统计输出 ------------------
cat("\n============================================================\n")
cat("========== 论文 Materials & Methods / Results 可用统计 ==========\n")
cat("============================================================\n\n")

cat("【Materials & Methods】\n")
cat("----------------------------------------\n")
cat("Gene Set Enrichment Analysis (GSEA) was performed using the clusterProfiler\n")
cat("package (v", as.character(packageVersion("clusterProfiler")), ") with the fgsea algorithm. All genes from DESeq2\n", sep = "")
cat("results (", comparison, ") were ranked by the Wald statistic (stat) without\n", sep = "")
cat("pre-filtering by significance thresholds. Gene IDs from DESeq2 results\n")
cat("(locus_tag format) were mapped to eggNOG annotations via the NCBI feature\n")
cat("table (GCA_011752425.2_MCOR1.1_feature_table.txt). KEGG pathway and GO term\n")
cat("gene sets were constructed from eggNOG-mapper annotations. Gene sets with size\n")
cat("between 10 and 500 were included. P-values were adjusted using the\n")
cat("Benjamini-Hochberg method, and pathways with padj < 0.05 were considered\n")
cat("significantly enriched.\n\n")

cat("【Results - GSEA (", comparison, ")】\n", sep = "")
cat("----------------------------------------\n")

if (!is.null(gse_kegg) && nrow(gse_kegg) > 0) {
  cat(sprintf("KEGG GSEA identified %d significantly enriched pathways (padj < 0.05).\n", nrow(gse_kegg)))
  cat("Top enriched pathways included:\n")
  k_df <- as.data.frame(gse_kegg)
  for (i in 1:min(5, nrow(k_df))) {
    desc <- ifelse("Description" %in% colnames(k_df), k_df$Description[i], k_df$ID[i])
    cat(sprintf("  - %s (NES = %.2f, padj = %.2e, gene set size = %d)\n",
                desc, k_df$NES[i], k_df$p.adjust[i], k_df$setSize[i]))
  }
} else {
  cat("No significantly enriched KEGG pathways were identified (padj < 0.05).\n")
}

cat("\n")
if (!is.null(gse_go) && nrow(gse_go) > 0) {
  cat(sprintf("GO GSEA identified %d significantly enriched terms (padj < 0.05).\n", nrow(gse_go)))
  go_df_out <- as.data.frame(gse_go)
  for (i in 1:min(5, nrow(go_df_out))) {
    desc <- ifelse("Description" %in% colnames(go_df_out), go_df_out$Description[i], go_df_out$ID[i])
    cat(sprintf("  - %s (NES = %.2f, padj = %.2e)\n", desc, go_df_out$NES[i], go_df_out$p.adjust[i]))
  }
} else {
  cat("No significantly enriched GO terms were identified (padj < 0.05).\n")
}

cat("\n============================================================\n")
cat("统计输出完毕。\n")
cat("【后续】将 comparison 和 res_csv 改为 T32_vs_Ctrl / T32_vs_T24 即可重复\n")
cat("============================================================\n")