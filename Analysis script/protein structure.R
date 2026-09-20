library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 1. 读取注释 ==========
header_line <- readLines("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", n = 1)
col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]
feat <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", 
                   stringsAsFactors = FALSE, skip = 1, header = FALSE, col.names = col_names)

eggnog <- read.delim("annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular", 
                     stringsAsFactors = FALSE, comment.char = "")
names(eggnog) <- gsub("^#", "", names(eggnog))

# 建立映射（注意：eggNOG 列名是 X.query，不是 query）
mrna_map <- feat %>% filter(feature == "mRNA", related_accession != "") %>% 
  select(locus_tag, protein_id = related_accession) %>% distinct()

sym_map <- eggnog %>% 
  select(query = X.query, Preferred_name) %>% 
  filter(Preferred_name != "" & Preferred_name != "-") %>% 
  distinct()

id_map <- mrna_map %>% 
  left_join(sym_map, by = c("protein_id" = "query")) %>% 
  filter(!is.na(Preferred_name)) %>% 
  distinct()

# 从 feature_table 的 mRNA 行获取官方蛋白长度（用于验证）
mrna_len <- feat %>% 
  filter(feature == "mRNA") %>% 
  select(locus_tag, official_aa = product_length) %>% 
  distinct()

# 目标基因
TARGET <- c("DNAJA4", "DNAJB5", "HSPBP1", "BAX", "HSPA5")
target_map <- id_map %>% 
  filter(Preferred_name %in% TARGET) %>% 
  group_by(Preferred_name) %>% slice(1) %>% ungroup()

# ========== 2. 读取基因组 ==========
read_fasta_sel <- function(f, chrs) {
  con <- if (grepl("\\.gz$", f)) gzfile(f, "r") else file(f, "r")
  seqs <- list(); cur <- NULL; buf <- character()
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (grepl("^>", line)) {
      if (!is.null(cur) && cur %in% chrs) seqs[[cur]] <- paste(buf, collapse = "")
      cur <- sub("^>(\\S+).*", "\\1", line); buf <- character()
    } else { buf <- c(buf, line) }
  }
  if (!is.null(cur) && cur %in% chrs) seqs[[cur]] <- paste(buf, collapse = "")
  close(con); return(seqs)
}

need_chr <- character()
for (s in target_map$Preferred_name) {
  loc <- target_map$locus_tag[target_map$Preferred_name == s]
  cds <- feat %>% filter(feature == "CDS", locus_tag == loc)
  if (nrow(cds) > 0) need_chr <- c(need_chr, cds$genomic_accession[1])
}
need_chr <- unique(need_chr)

cat("读取基因组... ")
genome <- read_fasta_sel("./genome/GCA_011752425.2_MCOR1.1_genomic.fna", need_chr)
cat(length(genome), "条序列\n")

# ========== 3. 三重验证 ==========
CODON_TABLE <- list(
  TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S", CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R", GGT="G", GGC="G", GGA="G", GGG="G"
)

translate_exact <- function(dna) {
  dna <- toupper(dna)
  atg <- regexpr("ATG", dna)
  if (atg < 0) return(list(prot = NULL, cds_len = nchar(dna), aa_len = 0))
  orf <- substr(dna, atg, nchar(dna))
  len <- nchar(orf) - (nchar(orf) %% 3)
  orf <- substr(orf, 1, len)
  codons <- sapply(seq(1, len, 3), function(i) substr(orf, i, i+2))
  aa <- sapply(codons, function(c) ifelse(is.null(CODON_TABLE[[c]]), "X", CODON_TABLE[[c]]))
  prot <- paste(aa, collapse = "")
  stop_pos <- regexpr("[*]", prot)
  if (stop_pos > 0) prot <- substr(prot, 1, stop_pos - 1)
  return(list(prot = prot, cds_len = nchar(dna), aa_len = nchar(prot)))
}

gtf_lines <- readLines("annotation/MCOR1.1.gtf")

cat("\n========== 三重验证报告 ==========\n")
cat(sprintf("%-10s %10s %10s %8s %10s %10s %s\n", 
            "Gene", "FT_CDS_bp", "GTF_bp", "Diff%", "Official", "My_aa", "Status"))

fasta_out <- character()

for (s in target_map$Preferred_name) {
  loc <- target_map$locus_tag[target_map$Preferred_name == s]
  
  # --- feature_table CDS ---
  cds_ft <- feat %>% filter(feature == "CDS", locus_tag == loc) %>% 
    select(chr = genomic_accession, start, end, strand) %>% arrange(start)
  if (nrow(cds_ft) == 0) next
  
  if (cds_ft$strand[1] == "-") cds_ft <- cds_ft %>% arrange(desc(start))
  seq_chr <- genome[[cds_ft$chr[1]]]
  pieces <- sapply(1:nrow(cds_ft), function(i) substr(seq_chr, cds_ft$start[i], cds_ft$end[i]))
  if (cds_ft$strand[1] == "-") {
    pieces <- sapply(pieces, function(x) chartr("ATCGatcg", "TAGCtagc", 
                                                paste(rev(strsplit(x, "")[[1]]), collapse = "")))
  }
  ft_dna <- toupper(paste(pieces, collapse = ""))
  ft_cds_len <- nchar(ft_dna)
  
  # --- GTF exon 长度 ---
  gid <- paste0("gene-", loc)
  exon_lines <- gtf_lines[grepl(paste0('transcript_id "', gid, '"'), gtf_lines) & grepl("\texon\t", gtf_lines)]
  gtf_len <- sum(sapply(exon_lines, function(line) {
    p <- strsplit(line, "\t")[[1]]
    as.numeric(p[5]) - as.numeric(p[4]) + 1
  }))
  
  # --- 翻译 ---
  res <- translate_exact(ft_dna)
  my_aa <- res$aa_len
  
  # --- Official 长度 from feature_table mRNA ---
  official_aa <- mrna_len$official_aa[mrna_len$locus_tag == loc]
  if (length(official_aa) == 0 || is.na(official_aa)) official_aa <- 0
  
  # --- 差异 ---
  diff_pct <- ifelse(gtf_len > 0, round(abs(ft_cds_len - gtf_len) / gtf_len * 100, 1), NA)
  
  status <- "✅ PASS"
  if (!is.na(diff_pct) && diff_pct > 5) status <- "⚠️ CDS vs exon 差异大"
  if (official_aa > 0 && abs(my_aa - official_aa) / official_aa > 0.15) status <- "⚠️ 蛋白长度 vs official 差异大"
  if (my_aa == 0) status <- "❌ 翻译失败"
  
  cat(sprintf("%-10s %10d %10d %7s%% %10d %10d %s\n", 
              s, ft_cds_len, gtf_len, ifelse(is.na(diff_pct), "-", diff_pct), 
              official_aa, my_aa, status))
  
  # --- 保存 FASTA ---
  if (my_aa > 0) {
    fasta_out <- c(fasta_out, 
                   paste0(">", s, " | locus=", loc, " | CDS=", ft_cds_len, "bp | aa=", my_aa),
                   res$prot)
  }
}

# 保存
dir.create("results/RNAseq_analysis/Results/structure", showWarnings = FALSE, recursive = TRUE)
writeLines(fasta_out, "results/RNAseq_analysis/Results/structure/core_proteins_exact.fasta")
cat("\n✅ 已保存: core_proteins_exact.fasta\n")





setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 标准密码子表
CODON_TABLE <- list(
  TTT="F", TTC="F", TTA="L", TTG="L", CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M", GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S", CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T", GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*", CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K", GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W", CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R", GGT="G", GGC="G", GGA="G", GGG="G"
)

# 找最长 ORF（从所有 ATG 开始，取最长）
find_longest_orf <- function(dna) {
  dna <- toupper(dna)
  atg_positions <- gregexpr("ATG", dna)[[1]]
  if (atg_positions[1] < 0) return(NULL)
  
  best_orf <- ""
  best_len <- 0
  best_start <- 0
  
  for (pos in atg_positions) {
    orf <- substr(dna, pos, nchar(dna))
    len <- nchar(orf) - (nchar(orf) %% 3)
    orf <- substr(orf, 1, len)
    codons <- sapply(seq(1, len, 3), function(i) substr(orf, i, i+2))
    aa <- sapply(codons, function(c) ifelse(is.null(CODON_TABLE[[c]]), "X", CODON_TABLE[[c]]))
    prot <- paste(aa, collapse = "")
    stop_pos <- regexpr("[*]", prot)
    if (stop_pos > 0) prot <- substr(prot, 1, stop_pos - 1)
    
    if (nchar(prot) > best_len) {
      best_len <- nchar(prot)
      best_orf <- prot
      best_start <- pos
    }
  }
  return(list(prot = best_orf, start = best_start, len = best_len))
}

# 读取 cDNA FASTA
lines <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA.fasta")

records <- list()
cur_name <- NULL
cur_seq <- character()

for (line in lines) {
  if (grepl("^>", line)) {
    if (!is.null(cur_name)) records[[cur_name]] <- paste(cur_seq, collapse = "")
    cur_name <- sub("^>(\\S+).*", "\\1", line)
    cur_seq <- character()
  } else { cur_seq <- c(cur_seq, line) }
}
if (!is.null(cur_name)) records[[cur_name]] <- paste(cur_seq, collapse = "")

# 翻译 + 验证
cat(sprintf("%-10s %10s %10s %s\n", "Gene", "cDNA_bp", "Protein_aa", "Status"))
fasta_out <- character()

for (s in names(records)) {
  dna <- toupper(records[[s]])
  res <- find_longest_orf(dna)
  
  if (is.null(res)) {
    cat(sprintf("%-10s %10d %10s %s\n", s, nchar(dna), "-", "❌ 无ATG"))
    next
  }
  
  # 与近缘物种预期长度核对
  expected <- switch(s,
                     "DNAJA4" = 400, "DNAJB5" = 300, "HSPBP1" = 300,
                     "BAG3" = 350, "BCL2" = 200, "CASP9" = 400,
                     "NFKBIA" = 300, "HSPA5" = 650, "BAX" = 190, 0
  )
  
  status <- ifelse(expected > 0 && abs(res$len - expected) / expected < 0.2, "✅ PASS", "⚠️ 核对")
  cat(sprintf("%-10s %10d %10d %s\n", s, nchar(dna), res$len, status))
  
  fasta_out <- c(fasta_out, 
                 paste0(">", s, " | cDNA=", nchar(dna), "bp | ORF_start=", res$start, " | aa=", res$len),
                 res$prot)
}

# 保存
writeLines(fasta_out, "results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein_v2.fasta")
cat("\n✅ 已保存: target_genes_protein_v2.fasta\n")





library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 1. 找到 HSPA5 的 locus_tag ==========
header_line <- readLines("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", n = 1)
col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]
feat <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", 
                   stringsAsFactors = FALSE, skip = 1, header = FALSE, col.names = col_names)

eggnog <- read.delim("annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular", 
                     stringsAsFactors = FALSE, comment.char = "")
names(eggnog) <- gsub("^#", "", names(eggnog))

mrna_map <- feat %>% filter(feature == "mRNA", related_accession != "") %>% 
  select(locus_tag, protein_id = related_accession) %>% distinct()

sym_map <- eggnog %>% 
  select(query = X.query, Preferred_name) %>% 
  filter(Preferred_name != "" & Preferred_name != "-") %>% 
  distinct()

id_map <- mrna_map %>% left_join(sym_map, by = c("protein_id" = "query")) %>% 
  filter(!is.na(Preferred_name)) %>% distinct()

# HSPA5 可能叫 HSPA5, GRP78, 或 BIP
hspa5_hits <- id_map %>% 
  filter(Preferred_name == "HSPA5" | grepl("GRP78|BIP", Preferred_name, ignore.case = TRUE))
cat("HSPA5 匹配结果:\n"); print(hspa5_hits)

if (nrow(hspa5_hits) == 0) {
  cat("❌ eggNOG 中未找到 HSPA5/GRP78/BIP\n")
} else {
  loc <- hspa5_hits$locus_tag[1]
  gid <- paste0("gene-", loc)
  
  # ========== 2. 从 GTF 提取外显子 ==========
  gtf_lines <- readLines("annotation/MCOR1.1.gtf")
  exon_lines <- gtf_lines[grepl(paste0('transcript_id "', gid, '"'), gtf_lines) & grepl("\texon\t", gtf_lines)]
  
  cat("HSPA5 外显子数:", length(exon_lines), "\n")
  
  if (length(exon_lines) > 0) {
    exons <- do.call(rbind, lapply(exon_lines, function(line) {
      p <- strsplit(line, "\t")[[1]]
      data.frame(chr = p[1], start = as.numeric(p[4]), end = as.numeric(p[5]), 
                 strand = p[7], stringsAsFactors = FALSE)
    }))
    exons <- exons[order(exons$start), ]
    if (exons$strand[1] == "-") exons <- exons[order(exons$start, decreasing = TRUE), ]
    
    chr <- exons$chr[1]
    
    # ========== 3. 读取基因组（只读该染色体）==========
    read_fasta_sel <- function(f, chrs) {
      con <- if (grepl("\\.gz$", f)) gzfile(f, "r") else file(f, "r")
      seqs <- list(); cur <- NULL; buf <- character()
      while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
        if (grepl("^>", line)) {
          if (!is.null(cur) && cur %in% chrs) seqs[[cur]] <- paste(buf, collapse = "")
          cur <- sub("^>(\\S+).*", "\\1", line); buf <- character()
        } else { buf <- c(buf, line) }
      }
      if (!is.null(cur) && cur %in% chrs) seqs[[cur]] <- paste(buf, collapse = "")
      close(con); return(seqs)
    }
    
    genome <- read_fasta_sel("./genome/GCA_011752425.2_MCOR1.1_genomic.fna", chr)
    seq_chr <- genome[[chr]]
    
    exon_seqs <- sapply(1:nrow(exons), function(i) substr(seq_chr, exons$start[i], exons$end[i]))
    if (exons$strand[1] == "-") {
      exon_seqs <- sapply(exon_seqs, function(s) chartr("ATCGatcg", "TAGCtagc", 
        paste(rev(strsplit(s, "")[[1]]), collapse = "")))
    }
    
    hspa5_cdna <- toupper(paste(exon_seqs, collapse = ""))
    cat("HSPA5 cDNA 长度:", nchar(hspa5_cdna), "bp\n")
    
    # ========== 4. 最长 ORF 翻译 ==========
    find_longest_orf <- function(dna) {
      dna <- toupper(dna)
      atg_positions <- gregexpr("ATG", dna)[[1]]
      if (atg_positions[1] < 0) return(NULL)
      best_orf <- ""; best_len <- 0; best_start <- 0
      for (pos in atg_positions) {
        orf <- substr(dna, pos, nchar(dna))
        len <- nchar(orf) - (nchar(orf) %% 3)
        orf <- substr(orf, 1, len)
        codons <- sapply(seq(1, len, 3), function(i) substr(orf, i, i+2))
        aa <- sapply(codons, function(c) {
          ct <- list(TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
                     ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
                     TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
                     ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
                     TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
                     AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
                     TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
                     AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G")
          ifelse(is.null(ct[[c]]), "X", ct[[c]])
        })
        prot <- paste(aa, collapse = "")
        stop_pos <- regexpr("[*]", prot)
        if (stop_pos > 0) prot <- substr(prot, 1, stop_pos - 1)
        if (nchar(prot) > best_len) { best_len <- nchar(prot); best_orf <- prot; best_start <- pos }
      }
      return(list(prot = best_orf, start = best_start, len = best_len))
    }
    
    res <- find_longest_orf(hspa5_cdna)
    cat("HSPA5: 蛋白长度", res$len, "aa\n")
    
    # ========== 5. 追加到 v2 ==========
    hspa5_fasta <- c(
      paste0(">HSPA5 | cDNA=", nchar(hspa5_cdna), "bp | ORF_start=", res$start, " | aa=", res$len),
      res$prot
    )
    existing <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein_v2.fasta")
    writeLines(c(existing, hspa5_fasta), "results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein_v2.fasta")
    cat("✅ HSPA5 已追加到 target_genes_protein_v2.fasta\n")
  }
}



setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 1. 删除错误的 342 aa HSPA5，保留正确的 658 aa ==========
lines <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein_v2.fasta")
clean <- character()
skip <- FALSE
hspa5_count <- 0

for (line in lines) {
  if (grepl("^>HSPA5", line)) {
    hspa5_count <- hspa5_count + 1
    if (hspa5_count == 1) {
      skip <- TRUE   # 跳过第一个（342 aa 错误版本）
      next
    } else {
      skip <- FALSE  # 保留第二个（658 aa 正确版本）
    }
  }
  if (skip && grepl("^>", line)) {
    skip <- FALSE  # 如果跳过时遇到其他基因 header，停止跳过
  }
  if (!skip) {
    clean <- c(clean, line)
  }
}

writeLines(clean, "results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein_v2.fasta")
cat("✅ 已删除 342 aa 错误 HSPA5，保留 658 aa 正确版本\n")

# ========== 2. 提取核心 5 蛋白到专用文件 ==========
core_genes <- c("DNAJA4", "DNAJB5", "HSPBP1", "HSPA5", "BAX")
core_lines <- character()

for (i in seq_along(lines)) {
  if (grepl("^>", lines[i])) {
    gene_name <- sub("^>(\\S+).*", "\\1", lines[i])
    if (gene_name %in% core_genes) {
      # 收集该基因的所有序列行
      seq_lines <- character()
      j <- i + 1
      while (j <= length(lines) && !grepl("^>", lines[j])) {
        seq_lines <- c(seq_lines, lines[j])
        j <- j + 1
      }
      # 只保留 658 aa 的 HSPA5（跳过 342 aa 的）
      if (gene_name == "HSPA5" && nchar(paste(seq_lines, collapse = "")) < 400) next
      core_lines <- c(core_lines, lines[i], seq_lines)
    }
  }
}

dir.create("results/RNAseq_analysis/Results/structure", showWarnings = FALSE, recursive = TRUE)
writeLines(core_lines, "results/RNAseq_analysis/Results/structure/core_5_proteins.fasta")
cat("✅ 核心 5 蛋白已保存到 core_5_proteins.fasta\n")

# ========== 3. 最终核对 ==========
cat("\n========== 最终核对 ==========\n")
final_lines <- readLines("results/RNAseq_analysis/Results/structure/core_5_proteins.fasta")
for (i in seq_along(final_lines)) {
  if (grepl("^>", final_lines[i])) {
    gene <- sub("^>(\\S+).*", "\\1", final_lines[i])
    seq <- paste(final_lines[(i+1):length(final_lines)][!grepl("^>", final_lines[(i+1):length(final_lines)])][1:100], collapse = "")
    # 找到该基因序列的准确长度
    seq_len <- 0
    j <- i + 1
    while (j <= length(final_lines) && !grepl("^>", final_lines[j])) {
      seq_len <- seq_len + nchar(final_lines[j])
      j <- j + 1
    }
    cat(sprintf("%-10s %6d aa\n", gene, seq_len))
  }
}




library(jsonlite)
library(ggplot2)
library(reshape2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure")

# 列出所有子文件夹
all_items <- list.files()
folders <- all_items[file.info(all_items)$isdir]
cat("发现", length(folders), "个文件夹:\n")
print(folders)

cat("\n========== 开始批量处理 ==========\n\n")

for (folder in folders) {
  
  # 提取基因名（去掉 fold_ 前缀和括号后缀）
  gene <- sub("^fold_", "", folder)
  gene <- sub("\\(.*", "", gene)
  gene <- toupper(gene)
  
  cat("处理:", gene, "| 文件夹:", folder, "\n")
  
  # 列出文件夹内所有文件
  files <- list.files(folder)
  cat("  文件列表:", paste(files, collapse = ", "), "\n")
  
  # 找关键文件
  sum_json_file <- files[grepl("summary_confidences_0\\.json$", files)]
  full_json_file <- files[grepl("full_data_0\\.json$", files)]
  cif_file <- files[grepl("model_0\\.cif$", files)]
  
  # --- 1. 生成 pLDDT 图 ---
  if (length(sum_json_file) == 1) {
    sum_path <- file.path(folder, sum_json_file)
    sum_json <- fromJSON(sum_path)
    plddt <- sum_json$plddt
    
    df_plddt <- data.frame(Residue = 1:length(plddt), pLDDT = plddt)
    
    p1 <- ggplot(df_plddt, aes(x = Residue, y = pLDDT)) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
      geom_hline(yintercept = 70, linetype = "dashed", color = "orange") +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      annotate("text", x = length(plddt)*0.95, y = 92, label = "Very high", color = "darkgreen", hjust = 1, size = 3) +
      annotate("text", x = length(plddt)*0.95, y = 72, label = "Confident", color = "orange", hjust = 1, size = 3) +
      ylim(0, 100) +
      theme_minimal() +
      labs(title = paste0(gene, " - Predicted LDDT"), x = "Residue", y = "pLDDT")
    
    out_plddt <- paste0(gene, "_pLDDT.pdf")
    ggsave(out_plddt, p1, width = 8, height = 4)  # 去掉 dpi
    cat("  ✅ 生成:", out_plddt, "| 残基数:", length(plddt), "\n")
  } else {
    cat("  ⚠️ 未找到 summary_confidences_0.json\n")
  }
  
  # --- 2. 生成 PAE 图 ---
  if (length(full_json_file) == 1) {
    full_path <- file.path(folder, full_json_file)
    full_json <- fromJSON(full_path)
    pae <- full_json$pae
    
    pae_df <- melt(pae)
    colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE")
    
    p2 <- ggplot(pae_df, aes(x = Residue_j, y = Residue_i, fill = PAE)) +
      geom_tile() +
      scale_fill_gradient(low = "darkgreen", high = "white", limits = c(0, 30)) +
      coord_fixed() +
      theme_minimal() +
      labs(title = paste0(gene, " - Predicted Aligned Error"), 
           x = "Scored Residue", y = "Aligned Residue")
    
    out_pae <- paste0(gene, "_PAE.pdf")
    ggsave(out_pae, p2, width = 6, height = 6)   # 去掉 dpi
    cat("  ✅ 生成:", out_pae, "| 矩阵大小:", nrow(pae), "x", ncol(pae), "\n")
  } else {
    cat("  ⚠️ 未找到 full_data_0.json\n")
  }
  
  # --- 3. 整理结构文件 ---
  if (length(cif_file) == 1) {
    cif_path <- file.path(folder, cif_file)
    out_cif <- paste0(gene, "_AF.cif")
    file.copy(cif_path, out_cif, overwrite = TRUE)
    cat("  ✅ 复制:", out_cif, "\n")
  } else {
    cat("  ⚠️ 未找到 model_0.cif\n")
  }
  
  cat("\n")
}

cat("========== 全部完成 ==========\n")
cat("structure 文件夹现在包含:\n")
final_files <- list.files(pattern = "(_pLDDT\\.png$|_PAE\\.png$|_AF\\.cif$)")
print(final_files)



library(jsonlite)

# ========== 配置路径 ==========
base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure"

folders <- list(
  DNAJA4_HSPA5 = file.path(base_dir, "fold_dnaja4_hspa5"),
  HSPA5_BAX    = file.path(base_dir, "fold_hspa5_bax")
)

# ========== 1. 打印每个文件夹的所有文件 ==========
cat("=" , rep("=", 59), "\n", sep = "")
cat("【1】文件夹内容扫描\n")
cat("=" , rep("=", 59), "\n", sep = "")

for (job_name in names(folders)) {
  folder_path <- folders[[job_name]]
  cat("\n📁", job_name, ":", folder_path, "\n")
  
  if (!dir.exists(folder_path)) {
    cat("   ⚠️  路径不存在！\n")
    next
  }
  
  # 递归列出所有文件
  all_files <- list.files(folder_path, recursive = TRUE, full.names = FALSE)
  all_files <- sort(all_files)
  
  for (f in all_files) {
    marker <- ""
    if (grepl("\\.pdb$", f, ignore.case = TRUE)) {
      marker <- "  ← PDB模型"
    } else if (grepl("pae|PAE", f)) {
      marker <- "  ← PAE文件"
    } else if (grepl("plddt|pLDDT", f)) {
      marker <- "  ← pLDDT图"
    } else if (grepl("confidence", f, ignore.case = TRUE) && grepl("\\.json$", f, ignore.case = TRUE)) {
      marker <- "  ← 置信度JSON"
    } else if (grepl("\\.png$", f, ignore.case = TRUE)) {
      marker <- "  ← 图片"
    }
    cat("   ├──", f, marker, "\n")
  }
}

# ========== 2. 提取关键指标（ipTM, pTM, pLDDT） ==========
cat("\n", "=" , rep("=", 59), "\n", sep = "")
cat("【2】关键指标提取\n")
cat("=" , rep("=", 59), "\n", sep = "")

for (job_name in names(folders)) {
  folder_path <- folders[[job_name]]
  cat("\n📊", job_name, ":\n")
  
  if (!dir.exists(folder_path)) {
    cat("   ⚠️  路径不存在\n")
    next
  }
  
  # 找所有 JSON 文件
  json_files <- list.files(folder_path, pattern = "\\.json$", recursive = TRUE, full.names = TRUE)
  conf_files <- json_files[grepl("confidence|summary", basename(json_files), ignore.case = TRUE)]
  
  found_metrics <- FALSE
  
  for (cf in conf_files) {
    tryCatch({
      data <- fromJSON(cf)
      cat("   读取:", basename(cf), "\n")
      
      # 尝试不同可能的键名
      iptm <- data$iptm %||% data$ipTM %||% data$interface_pTM
      ptm  <- data$ptm %||% data$pTM %||% data$predictedTMscore
      plddt <- data$mean_plddt %||% data$mean_pLDDT %||% data$plddt
      
      # 嵌套结构
      if (is.null(iptm) && !is.null(data$confidence)) {
        iptm <- data$confidence$iptm %||% data$confidence$ipTM
      }
      if (is.null(ptm) && !is.null(data$confidence)) {
        ptm <- data$confidence$ptm %||% data$confidence$pTM
      }
      if (is.null(plddt) && !is.null(data$confidence)) {
        plddt <- data$confidence$mean_plddt %||% data$confidence$plddt
      }
      
      # 如果还是找不到，打印顶层键
      if (is.null(iptm) && is.null(ptm)) {
        cat("   顶层键:", paste(names(data)[1:min(10, length(names(data)))], collapse = ", "), "\n")
      }
      
      if (!is.null(iptm)) {
        cat(sprintf("   ipTM (界面置信度): %.4f\n", iptm))
      }
      if (!is.null(ptm)) {
        cat(sprintf("   pTM  (整体置信度): %.4f\n", ptm))
      }
      if (!is.null(plddt)) {
        cat(sprintf("   pLDDT(平均): %.2f\n", plddt))
      }
      
      if (!is.null(iptm) || !is.null(ptm)) {
        found_metrics <- TRUE
      }
      
    }, error = function(e) {
      cat("   ⚠️  读取", basename(cf), "出错:", conditionMessage(e), "\n")
    })
  }
  
  if (!found_metrics) {
    cat("   ⚠️  未找到 ipTM/pTM 数据，请检查 JSON 文件内容\n")
  }
}

# ========== 3. 定位关键文件路径（供后续使用） ==========
cat("\n", "=" , rep("=", 59), "\n", sep = "")
cat("【3】关键文件路径汇总\n")
cat("=" , rep("=", 59), "\n", sep = "")

for (job_name in names(folders)) {
  folder_path <- folders[[job_name]]
  cat("\n📁", job_name, ":\n")
  
  # PDB 文件
  pdb_files <- list.files(folder_path, pattern = "\\.pdb$", recursive = TRUE, full.names = TRUE)
  cat("   PDB 模型:", length(pdb_files), "个\n")
  if (length(pdb_files) > 0) {
    for (p in head(pdb_files, 3)) {
      cat("      -", basename(p), "\n")
    }
  }
  
  # PAE 文件
  pae_files <- list.files(folder_path, pattern = "pae|PAE", recursive = TRUE, full.names = TRUE)
  cat("   PAE 文件:", length(pae_files), "个\n")
  if (length(pae_files) > 0) {
    for (p in head(pae_files, 3)) {
      cat("      -", basename(p), "\n")
    }
  }
  
  # pLDDT 图
  plddt_files <- list.files(folder_path, pattern = "plddt|pLDDT", recursive = TRUE, full.names = TRUE)
  cat("   pLDDT 图:", length(plddt_files), "个\n")
  if (length(plddt_files) > 0) {
    for (p in head(plddt_files, 3)) {
      cat("      -", basename(p), "\n")
    }
  }
}

cat("\n", "=" , rep("=", 59), "\n", sep = "")
cat("✅ 扫描完成。如果路径不存在，请确认文件夹名是否正确。\n")
cat("=" , rep("=", 59), "\n", sep = "")





library(jsonlite)
library(ggplot2)
library(reshape2)

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure"

# ========== 定义两个任务（根据你的实际序列长度修改）==========
jobs <- list(
  DNAJA4_HSPA5 = list(
    folder      = file.path(base_dir, "fold_dnaja4_hspa5"),
    name        = "DNAJA4_HSPA5",
    chain_a_len = 479,   # DNAJA4 长度，请根据你的实际序列修改
    chain_b_len = 658    # HSPA5 长度
  ),
  HSPA5_BAX = list(
    folder      = file.path(base_dir, "fold_hspa5_bax"),
    name        = "HSPA5_BAX",
    chain_a_len = 658,   # HSPA5 长度
    chain_b_len = 209    # BAX 长度
  )
)

# ========== 逐个处理 ==========
for (job in jobs) {
  
  cat("\n========================================\n")
  cat("处理:", job$name, "\n")
  cat("========================================\n")
  
  # --- 1. 读取 summary_confidences_0.json（提取 ipTM/pTM）---
  summary_files <- list.files(job$folder, pattern = "summary_confidences_0\\.json$", full.names = TRUE)
  
  if (length(summary_files) == 0) {
    cat("⚠️ 未找到 summary_confidences_0.json\n")
    next
  }
  
  summary_data <- fromJSON(summary_files[1])
  iptm_val <- summary_data$iptm
  ptm_val  <- summary_data$ptm
  cat("ipTM:", iptm_val, "| pTM:", ptm_val, "\n")
  
  # --- 2. 读取 full_data_0.json（提取 pLDDT 和 PAE）---
  json_files <- list.files(job$folder, pattern = "full_data_0\\.json$", full.names = TRUE)
  
  if (length(json_files) == 0) {
    cat("⚠️ 未找到 full_data_0.json\n")
    next
  }
  
  cat("读取:", basename(json_files[1]), "\n")
  data <- fromJSON(json_files[1])
  cat("JSON 顶层键:", paste(names(data), collapse = ", "), "\n")
  
  # --- 提取 pLDDT ---
  plddt <- NULL
  if (!is.null(data$atom_plddts)) {
    plddt <- data$atom_plddts
  } else if (!is.null(data$plddt)) {
    plddt <- data$plddt
  }
  
  if (!is.null(plddt)) {
    cat("pLDDT 长度:", length(plddt), "\n")
    cat("pLDDT 均值:", round(mean(plddt), 2), "\n")
  } else {
    cat("⚠️ 未找到 pLDDT\n")
  }
  
  # --- 提取 PAE ---
  pae <- NULL
  if (!is.null(data$pae)) {
    pae <- as.matrix(data$pae)
    cat("PAE 矩阵维度:", nrow(pae), "x", ncol(pae), "\n")
  } else {
    cat("⚠️ 未找到 PAE\n")
  }
  
  # --- 3. 绘制 pLDDT 曲线 ---
  if (!is.null(plddt)) {
    total_len <- length(plddt)
    chain_boundary <- job$chain_a_len
    
    plddt_df <- data.frame(
      Residue = 1:total_len,
      pLDDT   = plddt,
      Chain   = ifelse(1:total_len <= chain_boundary, "Chain A", "Chain B")
    )
    
    subtitle_text <- paste0("ipTM = ", round(iptm_val, 2), " | pTM = ", round(ptm_val, 2),
                            " | Mean pLDDT = ", round(mean(plddt), 1))
    
    p_plddt <- ggplot(plddt_df, aes(x = Residue, y = pLDDT, color = Chain)) +
      geom_line(linewidth = 0.8) +
      geom_hline(yintercept = 70, linetype = "dashed", color = "orange", alpha = 0.7) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "blue", alpha = 0.7) +
      annotate("text", x = total_len * 0.02, y = 72, label = "pLDDT = 70", color = "orange", hjust = 0, size = 3) +
      annotate("text", x = total_len * 0.02, y = 92, label = "pLDDT = 90", color = "blue", hjust = 0, size = 3) +
      geom_vline(xintercept = chain_boundary + 0.5, linetype = "dotted", color = "red") +
      annotate("text", x = chain_boundary, y = 102, label = "Chain A | Chain B", color = "red", hjust = 0.5, size = 3) +
      scale_color_manual(values = c("Chain A" = "#2E86AB", "Chain B" = "#A23B72")) +
      labs(title = paste(job$name, "- pLDDT Profile (Model 0)"),
           subtitle = subtitle_text,
           x = "Residue Index", y = "pLDDT") +
      theme_bw() +
      theme(legend.position = "bottom") +
      ylim(0, 105)
    
    out_plddt <- file.path(job$folder, paste0(job$name, "_pLDDT.pdf"))
    ggsave(out_plddt, p_plddt, width = 10, height = 5)
    cat("✅ pLDDT 图已保存:", basename(out_plddt), "\n")
  }
  
  # --- 4. 绘制 PAE 热图 ---
  if (!is.null(pae)) {
    pae_df <- melt(pae)
    colnames(pae_df) <- c("Residue_i", "Residue_j", "PAE")
    
    boundary <- job$chain_a_len
    
    p_pae <- ggplot(pae_df, aes(x = Residue_j, y = Residue_i, fill = PAE)) +
      geom_raster() +
      scale_fill_gradientn(colors = c("darkblue", "blue", "cyan", "yellow", "orange", "red"),
                           limits = c(0, 30), name = "PAE (Å)") +
      geom_vline(xintercept = boundary + 0.5, color = "white", linewidth = 0.5) +
      geom_hline(yintercept = boundary + 0.5, color = "white", linewidth = 0.5) +
      annotate("text", x = boundary/2, y = nrow(pae) * 0.98, label = "Chain A", color = "white", fontface = "bold", size = 4) +
      annotate("text", x = boundary + (ncol(pae)-boundary)/2, y = boundary * 0.02, label = "Chain B", color = "white", fontface = "bold", size = 4) +
      labs(title = paste(job$name, "- Predicted Aligned Error (Model 0)"),
           subtitle = "White crosshair = chain boundary; Dark = low error, Light = high error",
           x = "Scored Residue", y = "Aligned Residue") +
      coord_fixed() +
      theme_minimal() +
      theme(legend.position = "right")
    
    out_pae <- file.path(job$folder, paste0(job$name, "_PAE.pdf"))
    ggsave(out_pae, p_pae, width = 8, height = 7)
    cat("✅ PAE 图已保存:", basename(out_pae), "\n")
  }
  
  cat("📦 完成，文件保存在:", job$folder, "\n")
}

cat("\n========================================\n")
cat("全部完成！\n")
cat("========================================\n")


library(jsonlite)
library(ggplot2)

setwd("/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure")

# ========== 纯 R CIF 解析（提取 B-factor = pLDDT）==========
read_cif_plddt <- function(cif_path) {
  lines <- readLines(cif_path, warn = FALSE)
  
  # 找到 _atom_site 列名
  atom_site_lines <- grep("_atom_site\\.", lines, value = TRUE)
  col_names <- sub("_atom_site\\.", "", atom_site_lines)
  
  # 找关键列
  b_idx <- which(col_names == "B_iso_or_equiv")[1]
  seq_idx <- which(col_names == "label_seq_id")[1]
  if (is.na(seq_idx)) seq_idx <- which(col_names == "auth_seq_id")[1]
  
  # 找数据行（_atom_site 之后，非空、非#、非_ 开头的行）
  header_end <- max(grep("_atom_site\\.", lines))
  data_start <- header_end + 1
  while (data_start <= length(lines) && trimws(lines[data_start]) == "") {
    data_start <- data_start + 1
  }
  
  data_end <- data_start
  for (i in data_start:length(lines)) {
    line <- trimws(lines[i])
    if (line == "" || substr(line, 1, 1) == "#" || substr(line, 1, 1) == "_") break
    data_end <- i
  }
  
  # 解析
  res_b <- list()
  for (line in lines[data_start:data_end]) {
    line <- trimws(line)
    if (line == "") next
    # 处理引号包裹的空格
    line <- gsub('"([^"]*)"', '\\1', line)
    parts <- strsplit(line, "\\s+")[[1]]
    parts <- parts[parts != ""]
    
    if (length(parts) >= max(b_idx, seq_idx, na.rm = TRUE)) {
      res_id <- suppressWarnings(as.numeric(parts[seq_idx]))
      b_val <- suppressWarnings(as.numeric(parts[b_idx]))
      if (!is.na(res_id) && !is.na(b_val)) {
        key <- as.character(res_id)
        res_b[[key]] <- c(res_b[[key]], b_val)
      }
    }
  }
  
  sorted_ids <- sort(as.numeric(names(res_b)))
  sapply(sorted_ids, function(r) mean(res_b[[as.character(r)]]))
}

# ========== 5 个单体 ==========
genes <- c("dnaja4", "dnajb5", "hspbp1", "hspa5", "bax")

cat("========== 从 CIF 重新生成 pLDDT PDF ==========\n\n")

for (g in genes) {
  
  gene_upper <- toupper(g)
  cif_name <- paste0("fold_", g, "_AF.cif")
  cif_path <- file.path(paste0("fold_", g), cif_name)
  
  cat("处理:", gene_upper, "| CIF:", cif_path, "\n")
  
  if (!file.exists(cif_path)) {
    cat("  ⚠️ 文件不存在，跳过\n\n")
    next
  }
  
  # 读取 pLDDT
  plddt <- read_cif_plddt(cif_path)
  cat("  ℹ️ 残基数:", length(plddt), 
      "| pLDDT 范围:", paste(round(range(plddt), 1), collapse = " ~ "), "\n")
  
  # 绘图
  df <- data.frame(Residue = seq_along(plddt), pLDDT = plddt)
  
  p <- ggplot(df, aes(x = Residue, y = pLDDT)) +
    geom_line(color = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "darkgreen") +
    geom_hline(yintercept = 70, linetype = "dashed", color = "orange") +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    annotate("text", x = length(plddt)*0.95, y = 92, label = "Very high", color = "darkgreen", hjust = 1, size = 3) +
    annotate("text", x = length(plddt)*0.95, y = 72, label = "Confident", color = "orange", hjust = 1, size = 3) +
    ylim(0, 100) +
    theme_minimal() +
    labs(title = paste0(gene_upper, " - Predicted LDDT"), 
         subtitle = "Extracted from CIF B-factor (residue-level)",
         x = "Residue", y = "pLDDT")
  
  # 保存到根目录，覆盖旧文件
  out_file <- paste0(gene_upper, "_pLDDT.pdf")
  ggsave(out_file, p, width = 8, height = 4)
  cat("  ✅ 生成:", out_file, "\n\n")
}

cat("========== 完成 ==========\n")





library(jsonlite)

# 定义两个复合物的路径
paths <- c(
  dnaja4_hspa5 = "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/fold_dnaja4_hspa5/summary_confidences.json",
  hspa5_bax = "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/fold_hspa5_bax/summary_confidences.json"
)

# 提取关键数值
results <- lapply(paths, function(p) {
  if (file.exists(p)) {
    j <- fromJSON(p)
    data.frame(
      complex = basename(dirname(p)),
      ipTM = ifelse(is.null(j$iptm), NA, j$iptm),
      pTM = ifelse(is.null(j$ptm), NA, j$ptm),
      pDockQ = ifelse(is.null(j$pdockq), NA, j$pdockq),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      complex = basename(dirname(p)),
      ipTM = NA, pTM = NA, pDockQ = NA,
      stringsAsFactors = FALSE
    )
  }
})

do.call(rbind, results)



library(jsonlite)

# 正确的文件路径（_0 是排名第一的模型）
path1 <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/fold_dnaja4_hspa5/fold_dnaja4_hspa5_summary_confidences_0.json"
path2 <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/fold_hspa5_bax/fold_hspa5_bax_summary_confidences_0.json"

# 提取函数
extract_scores <- function(path, name) {
  j <- fromJSON(path)
  data.frame(
    complex = name,
    ipTM = j$ipTM,
    pTM = j$pTM,
    pDockQ = ifelse(is.null(j$pDockQ), NA, j$pDockQ),
    stringsAsFactors = FALSE
  )
}

r1 <- extract_scores(path1, "DNAJA4-HSPA5")
r2 <- extract_scores(path2, "HSPA5-BAX")

results <- rbind(r1, r2)
print(results)




library(jsonlite)

# ========== 提取真实数据 ==========
path1 <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/fold_dnaja4_hspa5/fold_dnaja4_hspa5_summary_confidences_0.json"
path2 <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/fold_hspa5_bax/fold_hspa5_bax_summary_confidences_0.json"

j1 <- fromJSON(path1)
j2 <- fromJSON(path2)

complexes <- c("DNAJA4-HSPA5", "HSPA5-BAX")
ipTM <- c(j1$iptm, j2$iptm)

# ========== 画柱状图 ==========
out_path <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/complex_quality_scores.png"

png(out_path, width = 8, height = 6, units = "in", res = 600)

par(mar = c(5, 5, 4, 2) + 0.1, family = "Arial")
x <- barplot(ipTM, names.arg = complexes, col = "#3070B3", 
             ylim = c(0, 1.0), ylab = "ipTM Score", 
             main = "AlphaFold-Multimer Complex Prediction Quality",
             cex.main = 1.3, cex.lab = 1.1, cex.names = 1.1,
             border = "black", space = 0.6, lwd = 0.8)

# 阈值线
abline(h = 0.75, col = "#C74375", lty = 2, lwd = 2)
abline(h = 0.50, col = "grey60", lty = 2, lwd = 1.5)

# 标注数值
text(x, ipTM + 0.03, labels = sprintf("%.2f", ipTM), 
     cex = 1.3, font = 2, col = "black")

# 图例
legend("topright", 
       legend = c("High confidence (0.75)", "Acceptable (0.50)", "Actual scores"), 
       col = c("#C74375", "grey60", "#3070B3"), 
       lty = c(2, 2, NA), lwd = c(2, 1.5, NA),
       pch = c(NA, NA, 15), pt.cex = c(NA, NA, 2),
       bty = "n", cex = 1.0)

dev.off()
cat("柱状图已保存至：", out_path, "\n")




library(dplyr)
library(tidyr)

# ========== 路径 ==========
base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results"
out_csv  <- file.path(base_dir, "structure", "core_5_genes_RNAseq_summary.csv")

# 映射表（从之前成功的运行结果中确认）
id_map <- c(
  "gene-MCOR_13638" = "DNAJA4",
  "gene-MCOR_7799"  = "DNAJB5",
  "gene-MCOR_26887" = "HSPBP1",
  "gene-MCOR_56583" = "HSPA5",
  "gene-MCOR_18"    = "BAX"
)

files <- list(
  "T24_vs_Ctrl" = file.path(base_dir, "DESeq2_T24_vs_Ctrl.csv"),
  "T32_vs_Ctrl" = file.path(base_dir, "DESeq2_T32_vs_Ctrl.csv"),
  "T32_vs_T24"  = file.path(base_dir, "DESeq2_T32_vs_T24.csv")
)

# ========== 提取函数 ==========
extract <- function(path, contrast) {
  res <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(res)[1] <- "gene_id"
  res$gene_symbol <- id_map[res$gene_id]
  df <- res[res$gene_id %in% names(id_map), 
            c("gene_id", "gene_symbol", "baseMean", "log2FoldChange", "padj")]
  df$contrast <- contrast
  return(df)
}

# ========== 合并 ==========
all_data <- bind_rows(lapply(names(files), function(nm) extract(files[[nm]], nm)))

# 转成宽格式（每个基因一行，各对比组分开列）
wide <- all_data %>%
  select(gene_id, gene_symbol, contrast, log2FoldChange, padj) %>%
  pivot_wider(
    id_cols = c(gene_id, gene_symbol),
    names_from = contrast,
    values_from = c(log2FoldChange, padj)
  )

# 按目标基因顺序排列
wide$gene_symbol <- factor(wide$gene_symbol, levels = c("DNAJA4", "DNAJB5", "HSPBP1", "HSPA5", "BAX"))
wide <- wide[order(wide$gene_symbol), ]

# 保存
write.csv(wide, out_csv, row.names = FALSE)
cat("CSV 已保存至：", out_csv, "\n")
print(wide)




# ========== 验证脚本 ==========
# 读取你的结构基因序列文件（从文件列表里看到的）
cds_csv <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA_with_sequence.csv"
prot_fasta <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein.fasta"

cat("=== 方法1：检查 CDS CSV ===\n")
if (file.exists(cds_csv)) {
  cds <- read.csv(cds_csv, stringsAsFactors = FALSE)
  print(head(cds, 3))
  cat("\n列名：", colnames(cds), "\n")
  # 看看有没有 gene_id 或基因名列
  for (col in colnames(cds)) {
    hits <- grep("DNAJA4|DNAJB5|HSPBP1|HSPA5|BAX", cds[[col]], ignore.case = TRUE, value = TRUE)
    if (length(hits) > 0) cat("列 [", col, "] 找到匹配：", head(hits, 3), "\n")
  }
}

cat("\n=== 方法2：检查蛋白 FASTA 的 header ===\n")
if (file.exists(prot_fasta)) {
  lines <- readLines(prot_fasta, n = 50)  # 读前50行
  headers <- lines[grep("^>", lines)]
  print(headers)
}





library(dplyr)

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results"
anno_path <- file.path(base_dir, "gene_annotation_table.csv")

anno <- read.csv(anno_path, stringsAsFactors = FALSE, check.names = FALSE)
targets <- c("DNAJA4", "DNAJB5", "HSPBP1", "HSPA5", "BAX", "BAG3", "BCL2", "CASP9", "NFKBIA")

# 找到所有匹配
hits <- anno[anno$Preferred_name %in% targets, c("gene_id", "Preferred_name")]
hits <- hits[!duplicated(hits), ]
print(hits[order(hits$Preferred_name), ])













library(dplyr)
library(ggplot2)

grades_file <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/DNAJA4_ConSurf/DNAJA4_AF_consurf_grades.txt"

# ========== 1. 智能读取（不硬编码列数）==========
lines <- readLines(grades_file)

# 找表头行（包含 POS 的行）
header_idx <- grep("^\\s*POS\\s+SEQ", lines)[1]
cat("表头在第", header_idx, "行\n")

# 读取数据，让 R 自动推断列数和列名
grades <- read.table(text = lines[header_idx:length(lines)], 
                     header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE, fill = TRUE)

cat("实际列名:", colnames(grades), "\n")
cat("列数:", ncol(grades), "\n")
cat("行数:", nrow(grades), "\n")

# ========== 2. 提取关键列（兼容不同版本）==========
# 找 POS, SEQ, SCORE, COLOR 列（不管它们在哪个位置）
pos_col   <- grep("^(POS|pos)$", colnames(grades), value = TRUE)[1]
seq_col   <- grep("^(SEQ|seq|SEQUENCE)$", colnames(grades), value = TRUE)[1]
score_col <- grep("^(SCORE|score|ConSurf.Score)$", colnames(grades), value = TRUE)[1]
color_col <- grep("^(COLOR|color|ConSurf.Grade)$", colnames(grades), value = TRUE)[1]

cat("POS列:", pos_col, "\n")
cat("SEQ列:", seq_col, "\n")
cat("SCORE列:", score_col, "\n")
cat("COLOR列:", color_col, "\n")

# 重命名简化
grades <- grades %>%
  select(POS = all_of(pos_col), 
         SEQ = all_of(seq_col), 
         SCORE = all_of(score_col), 
         COLOR = all_of(color_col))

# 强制转数值
grades$POS   <- suppressWarnings(as.numeric(grades$POS))
grades$SCORE <- suppressWarnings(as.numeric(grades$SCORE))
grades$COLOR <- suppressWarnings(as.numeric(grades$COLOR))

# 去掉 NA
grades <- grades %>% filter(!is.na(POS), !is.na(SCORE))

cat("\n总残基数:", nrow(grades), "\n")
cat("保守性分数范围:", round(range(grades$SCORE, na.rm = TRUE), 2), "\n")

# ========== 3. 找 HPD motif ==========
hpd_candidates <- which(grades$SEQ == "H" & 
                          dplyr::lead(grades$SEQ, 1) == "P" & 
                          dplyr::lead(grades$SEQ, 2) == "D")

if (length(hpd_candidates) > 0) {
  hpd_idx <- hpd_candidates[1]
  cat("\n✅ 找到 HPD motif 位置:", hpd_idx, "-", hpd_idx+2, "\n")
  cat("对应残基:", paste(grades$SEQ[hpd_idx:(hpd_idx+2)], collapse = ""), "\n")
  print(grades[hpd_idx:(hpd_idx+2), ])
} else {
  cat("\n⚠️ 未自动找到 HPD\n")
  hpd_idx <- integer(0)
}

# ========== 4. 区域定义 + 汇总 ==========
grades <- grades %>%
  mutate(
    Region = case_when(
      POS <= 70  ~ "J-domain",
      POS <= 100 ~ "Linker",
      TRUE       ~ "C-terminal"
    ),
    Region = factor(Region, levels = c("J-domain", "Linker", "C-terminal"))
  )

summary_stats <- grades %>%
  group_by(Region) %>%
  summarise(
    n = n(),
    mean_score = round(mean(SCORE, na.rm = TRUE), 2),
    median_color = round(median(COLOR, na.rm = TRUE), 1),
    .groups = "drop"
  )
cat("\n========== 区域保守性汇总 ==========\n")
print(summary_stats)

# ========== 5. 保存关键残基表 ==========
key_table <- grades %>% filter(POS <= 100) %>% select(POS, SEQ, SCORE, COLOR, Region)
write.csv(key_table, 
          "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/DNAJA4_consurf_key_residues.csv", 
          row.names = FALSE)
cat("\n✅ 关键残基表已保存\n")

# ========== 6. 画图 ==========
p <- ggplot(grades, aes(x = POS, y = SCORE)) +
  annotate("rect", xmin = -Inf, xmax = 70.5, ymin = -Inf, ymax = Inf, fill = "#8E44AD", alpha = 0.08) +
  annotate("rect", xmin = 70.5, xmax = 100.5, ymin = -Inf, ymax = Inf, fill = "#E74C3C", alpha = 0.08) +
  annotate("rect", xmin = 100.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#3498DB", alpha = 0.08) +
  geom_line(color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = COLOR), size = 1) +
  scale_color_gradientn(colors = c("cyan", "white", "magenta", "purple"), limits = c(1, 9), name = "Grade") +
  {if (length(hpd_idx) > 0) 
    geom_point(data = grades[hpd_idx:(hpd_idx+2), ], aes(x = POS, y = SCORE), 
               color = "red", size = 3, shape = 21, stroke = 1.5)
  } +
  annotate("text", x = 35,  y = max(grades$SCORE, na.rm = TRUE)*1.08, label = "J-domain", color = "#8E44AD", fontface = "bold", size = 3.5) +
  annotate("text", x = 85,  y = max(grades$SCORE, na.rm = TRUE)*1.08, label = "Linker", color = "#E74C3C", fontface = "bold", size = 3.5) +
  annotate("text", x = 290, y = max(grades$SCORE, na.rm = TRUE)*1.08, label = "C-terminal", color = "#3498DB", fontface = "bold", size = 3.5) +
  labs(title = "DNAJA4 - Evolutionary Conservation Profile (ConSurf)",
       subtitle = ifelse(length(hpd_idx) > 0, "HPD motif marked in red", ""),
       x = "Residue Position", y = "Conservation Score") +
  theme_bw() +
  ylim(NA, max(grades$SCORE, na.rm = TRUE)*1.15)

ggsave("/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/DNAJA4_consurf_scores.pdf", 
       p, width = 10, height = 5)
cat("✅ 保守性曲线图已保存\n")




library(dplyr)
library(ggplot2)

gene_name <- "HSPA5"
base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure"
grades_file <- file.path(base_dir, paste0(gene_name, "_ConSurf"), 
                         paste0(gene_name, "_AF_consurf_grades.txt"))

# 读取（通用解析，兼容任何格式）
lines <- readLines(grades_file)
header_idx <- grep("^\\s*POS\\s+SEQ", lines)[1]
data_lines <- lines[(header_idx + 1):length(lines)]
data_lines <- data_lines[nchar(trimws(data_lines)) > 0]

parsed <- lapply(data_lines, function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(parts) >= 5) {
    data.frame(POS = as.numeric(parts[1]), SEQ = parts[2],
               SCORE = as.numeric(parts[4]), COLOR = as.numeric(parts[5]),
               stringsAsFactors = FALSE)
  } else NULL
})

grades <- bind_rows(parsed) %>% filter(!is.na(POS), !is.na(SCORE))

# HSPA5 结构域：NBD ~1-386, Hinge ~387-420, SBD ~421-658
grades <- grades %>% mutate(
  Region = case_when(
    POS <= 386 ~ "NBD",
    POS <= 420 ~ "Hinge",
    TRUE       ~ "SBD"
  ),
  Region = factor(Region, levels = c("NBD", "Hinge", "SBD"))
)

summary_stats <- grades %>% group_by(Region) %>% 
  summarise(n = n(), mean_score = round(mean(SCORE), 2), .groups = "drop")
print(summary_stats)

# 画图
p <- ggplot(grades, aes(x = POS, y = SCORE)) +
  annotate("rect", xmin = -Inf, xmax = 386.5, ymin = -Inf, ymax = Inf, fill = "#2E86AB", alpha = 0.08) +
  annotate("rect", xmin = 386.5, xmax = 420.5, ymin = -Inf, ymax = Inf, fill = "#E74C3C", alpha = 0.08) +
  annotate("rect", xmin = 420.5, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#27AE60", alpha = 0.08) +
  geom_line(color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = COLOR), size = 1) +
  scale_color_gradientn(colors = c("cyan", "white", "magenta", "purple"), limits = c(1, 9), name = "Grade") +
  annotate("text", x = 200, y = max(grades$SCORE)*1.08, label = "NBD", color = "#2E86AB", fontface = "bold", size = 3.5) +
  annotate("text", x = 403, y = max(grades$SCORE)*1.08, label = "Hinge", color = "#E74C3C", fontface = "bold", size = 3.5) +
  annotate("text", x = 550, y = max(grades$SCORE)*1.08, label = "SBD", color = "#27AE60", fontface = "bold", size = 3.5) +
  labs(title = "HSPA5 - Evolutionary Conservation (ConSurf)",
       x = "Residue Position", y = "Conservation Score") +
  theme_bw()

ggsave(file.path(base_dir, "HSPA5_consurf_scores.pdf"), p, width = 10, height = 5)
cat







library(dplyr)

# ========== 改这里 1：文件路径 ==========
grades_file <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/HSPA5_ConSurf/HSPA5_AF_consurf_grades.txt"

# 读取（通用解析，与 DNAJA4 完全相同）
lines <- readLines(grades_file)
header_idx <- grep("^\\s*POS\\s+SEQ", lines)[1]
data_lines <- lines[(header_idx + 1):length(lines)]
data_lines <- data_lines[nchar(trimws(data_lines)) > 0]

parsed <- lapply(data_lines, function(line) {
  parts <- strsplit(trimws(line), "\\s+")[[1]]
  if (length(parts) >= 5) {
    data.frame(
      POS = as.numeric(parts[1]), SEQ = parts[2],
      SCORE = as.numeric(parts[4]), COLOR = as.numeric(parts[5]),
      stringsAsFactors = FALSE
    )
  } else NULL
})

grades <- bind_rows(parsed) %>% 
  filter(!is.na(POS), !is.na(SCORE))

# ========== 改这里 2：结构域定义 ==========
grades <- grades %>% mutate(
  Region = case_when(
    POS <= 386 ~ "NBD",
    POS <= 420 ~ "Hinge",
    TRUE       ~ "SBD"
  )
)

# ========== 改这里 3：保存路径 ==========
out_csv <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/HSPA5_consurf_key_residues.csv"
write.csv(grades, out_csv, row.names = FALSE)
cat("✅ HSPA5 CSV 已保存:", out_csv, "\n总行数:", nrow(grades), "\n")

# 区域汇总（论文要用的数字）
summary_stats <- grades %>% 
  group_by(Region) %>% 
  summarise(n = n(), mean_score = round(mean(SCORE), 2), .groups = "drop")
print(summary_stats)





cabs_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure/DNAJA4_CABS"

cat("=== DNAJA4_CABS 根目录内容 ===\n")
print(list.files(cabs_dir, full.names = FALSE))

cat("\n=== 递归列出所有子目录和文件 ===\n")
all_files <- list.files(cabs_dir, recursive = TRUE, full.names = FALSE)
for (f in all_files) {
  cat(f, "\n")
}

cat("\n=== 带完整路径 ===\n")
print(list.files(cabs_dir, recursive = TRUE, full.names = TRUE))






library(dplyr)
library(ggplot2)
library(readr)

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure"
cabs_dir <- file.path(base_dir, "DNAJA4_CABS")

# ---------- 1. 解析 RMSF ----------
lines <- readLines(file.path(cabs_dir, "plots", "RMSF.csv"))
data_lines <- lines[grep("^A[0-9]+", lines)]
parsed <- strsplit(data_lines, "\t")

rmsf <- data.frame(
  POS = as.numeric(gsub("A", "", sapply(parsed, `[`, 1))),
  RMSF = as.numeric(sapply(parsed, `[`, 2)),
  stringsAsFactors = FALSE
)

# ---------- 2. 结构域标注 ----------
rmsf <- rmsf %>%
  mutate(
    Region = case_when(
      POS <= 70  ~ "J-domain",
      POS <= 100 ~ "Linker",
      TRUE       ~ "C-terminal"
    ),
    Region = factor(Region, levels = c("J-domain", "Linker", "C-terminal")),
    High_Flex = RMSF > (mean(RMSF) + sd(RMSF))
  )

# ---------- 3. 高柔性位点 ----------
high_flex <- rmsf %>% filter(High_Flex == TRUE) %>% arrange(desc(RMSF))

# ---------- 4. 画图（与 HSPA5 统一风格）----------
p <- ggplot(rmsf, aes(x = POS, y = RMSF)) +
  # 背景色块：三个结构域
  annotate("rect", xmin = -Inf, xmax = 70.5,  ymin = -Inf, ymax = Inf, fill = "#8E44AD", alpha = 0.08) +
  annotate("rect", xmin = 70.5,  xmax = 100.5, ymin = -Inf, ymax = Inf, fill = "#E74C3C", alpha = 0.08) +
  annotate("rect", xmin = 100.5, xmax = Inf,   ymin = -Inf, ymax = Inf, fill = "#2A9D8F", alpha = 0.08) +
  # 数据线和点
  geom_line(color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = Region), size = 0.8, alpha = 0.8) +
  # 手动指定颜色，与背景块一致
  scale_color_manual(values = c("J-domain" = "#8E44AD", "Linker" = "#E74C3C", "C-terminal" = "#3498DB")) +
  # 高柔性位点（品红空心圆）
  geom_point(data = high_flex, aes(x = POS, y = RMSF), color = "magenta", size = 1.5, shape = 21, fill = NA, stroke = 0.8) +
  # 结构域标签
  annotate("text", x = 35,  y = max(rmsf$RMSF)*1.08, label = "J-domain",   color = "#8E44AD", fontface = "bold", size = 3.5) +
  annotate("text", x = 85,  y = max(rmsf$RMSF)*1.08, label = "Linker",     color = "#E74C3C", fontface = "bold", size = 3.5) +
  annotate("text", x = 290, y = max(rmsf$RMSF)*1.08, label = "C-terminal", color = "#3498DB", fontface = "bold", size = 3.5) +
  # 高柔性位点图例说明（右上角，与 HSPA5 统一位置）
  annotate("point", x = 450, y = max(rmsf$RMSF)*0.92, color = "magenta", size = 2, shape = 21, stroke = 0.8) +
  annotate("text", x = 465, y = max(rmsf$RMSF)*0.92, label = "High flexibility\n(> mean + 1SD)", 
           hjust = 0, size = 2.5, color = "magenta") +
  # 坐标轴和标题
  labs(title = "DNAJA4 - Residue Flexibility (CABS-flex RMSF)",
       x = "Residue Position", y = "RMSF (Å)") +
  theme_bw() +
  # 图例设置：与 HSPA5 统一
  theme(legend.position = c(0.88, 0.25), 
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.4, "cm")) +
  guides(color = guide_legend(title = "Domain"))

ggsave(file.path(base_dir, "DNAJA4_CABS_RMSF_v2.pdf"), p, width = 10, height = 5, dpi = 300)
cat("\n✅ 带图注的图已保存: DNAJA4_CABS_RMSF_v2.pdf\n")

library(dplyr)
library(ggplot2)
library(readr)

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure"
cabs_dir <- file.path(base_dir, "DNAJA4_CABS")

# ---------- 1. 正确解析 RMSF.csv ----------
lines <- readLines(file.path(cabs_dir, "plots", "RMSF.csv"))

# 跳过版本号头行，提取数据行（格式：A2\t0.835）
data_lines <- lines[grep("^A[0-9]+", lines)]
parsed <- strsplit(data_lines, "\t")

rmsf <- data.frame(
  POS = as.numeric(gsub("A", "", sapply(parsed, `[`, 1))),
  RMSF = as.numeric(sapply(parsed, `[`, 2)),
  stringsAsFactors = FALSE
)

cat("解析成功，共", nrow(rmsf), "个残基\n")
head(rmsf)

# ---------- 2. 按结构域标注 ----------
rmsf <- rmsf %>%
  mutate(
    Region = case_when(
      POS <= 70  ~ "J-domain",
      POS <= 100 ~ "Linker",
      TRUE       ~ "C-terminal"
    ),
    Region = factor(Region, levels = c("J-domain", "Linker", "C-terminal")),
    High_Flex = RMSF > (mean(RMSF) + sd(RMSF))
  )

# ---------- 3. 区域统计 ----------
rmsf_summary <- rmsf %>%
  group_by(Region) %>%
  summarise(
    n = n(),
    mean_RMSF = round(mean(RMSF), 3),
    median_RMSF = round(median(RMSF), 3),
    max_RMSF = round(max(RMSF), 3),
    sd_RMSF = round(sd(RMSF), 3),
    .groups = "drop"
  )
cat("\n=== DNAJA4 CABS-flex 区域 RMSF (Å) ===\n")
print(rmsf_summary)

# ---------- 4. 高柔性位点 ----------
high_flex <- rmsf %>% filter(High_Flex == TRUE) %>% arrange(desc(RMSF))
cat("\n=== 高柔性位点 Top 10 ===\n")
print(head(high_flex, 10))

# ---------- 5. 画图 ----------
p <- ggplot(rmsf, aes(x = POS, y = RMSF)) +
  annotate("rect", xmin = -Inf, xmax = 70.5,  ymin = -Inf, ymax = Inf, fill = "#8E44AD", alpha = 0.08) +
  annotate("rect", xmin = 70.5,  xmax = 100.5, ymin = -Inf, ymax = Inf, fill = "#E74C3C", alpha = 0.08) +
  annotate("rect", xmin = 100.5, xmax = Inf,   ymin = -Inf, ymax = Inf, fill = "#27AE60", alpha = 0.08) +
  geom_line(color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = Region), size = 0.8, alpha = 0.8) +
  geom_point(data = high_flex, aes(x = POS, y = RMSF), color = "magenta", size = 1.5, shape = 21, fill = NA, stroke = 0.8) +
  annotate("text", x = 35,  y = max(rmsf$RMSF)*1.08, label = "J-domain",   color = "#8E44AD", fontface = "bold", size = 3.5) +
  annotate("text", x = 85,  y = max(rmsf$RMSF)*1.08, label = "Linker",     color = "#E74C3C", fontface = "bold", size = 3.5) +
  annotate("text", x = 290, y = max(rmsf$RMSF)*1.08, label = "C-terminal", color = "#27AE60", fontface = "bold", size = 3.5) +
  labs(title = "DNAJA4 - Residue Flexibility (CABS-flex RMSF)",
       x = "Residue Position", y = "RMSF (Å)") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(base_dir, "DNAJA4_CABS_RMSF.pdf"), p, width = 10, height = 5, dpi = 300)
cat("\n✅ 图已保存: DNAJA4_CABS_RMSF.pdf\n")

# ---------- 6. 保存数据 ----------
write_csv(rmsf, file.path(cabs_dir, "DNAJA4_CABS_RMSF_processed.csv"))
write_csv(rmsf_summary, file.path(cabs_dir, "DNAJA4_CABS_RMSF_summary.csv"))
write_csv(high_flex, file.path(cabs_dir, "DNAJA4_CABS_high_flexibility_sites.csv"))
cat("✅ 数据已保存\n")



library(dplyr)
library(ggplot2)
library(readr)

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/structure"
hsp_dir <- file.path(base_dir, "HSPA5_CABS")

# ---------- 1. 解析 RMSF ----------
lines <- readLines(file.path(hsp_dir, "plots", "RMSF.csv"))
data_lines <- lines[grep("^A[0-9]+", lines)]
parsed <- strsplit(data_lines, "\t")

rmsf <- data.frame(
  POS = as.numeric(gsub("A", "", sapply(parsed, `[`, 1))),
  RMSF = as.numeric(sapply(parsed, `[`, 2)),
  stringsAsFactors = FALSE
)

# ---------- 2. 结构域标注 ----------
rmsf <- rmsf %>%
  mutate(
    Region = case_when(
      POS <= 386 ~ "NBD",
      POS <= 420 ~ "Hinge",
      TRUE       ~ "SBD"
    ),
    Region = factor(Region, levels = c("NBD", "Hinge", "SBD")),
    High_Flex = RMSF > (mean(RMSF) + sd(RMSF))
  )

# ---------- 3. 区域统计 ----------
rmsf_summary <- rmsf %>%
  group_by(Region) %>%
  summarise(
    n = n(),
    mean_RMSF = round(mean(RMSF), 3),
    median_RMSF = round(median(RMSF), 3),
    max_RMSF = round(max(RMSF), 3),
    sd_RMSF = round(sd(RMSF), 3),
    .groups = "drop"
  )
cat("=== HSPA5 CABS-flex 区域 RMSF (Å) ===\n")
print(rmsf_summary)

# ---------- 4. 高柔性位点 ----------
high_flex <- rmsf %>% filter(High_Flex == TRUE) %>% arrange(desc(RMSF))
cat("\n=== 高柔性位点 Top 10 ===\n")
print(head(high_flex, 10))

# ---------- 5. 画图 ----------
p <- ggplot(rmsf, aes(x = POS, y = RMSF)) +
  annotate("rect", xmin = -Inf, xmax = 386.5,  ymin = -Inf, ymax = Inf, fill = "#2E86AB", alpha = 0.08) +
  annotate("rect", xmin = 386.5,  xmax = 420.5, ymin = -Inf, ymax = Inf, fill = "#E74C3C", alpha = 0.08) +
  annotate("rect", xmin = 420.5, xmax = Inf,   ymin = -Inf, ymax = Inf, fill = "#27AE60", alpha = 0.08) +
  geom_line(color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = Region), size = 0.8, alpha = 0.8) +
  geom_point(data = high_flex, aes(x = POS, y = RMSF), color = "magenta", size = 1.5, shape = 21, fill = NA, stroke = 0.8) +
  # 加在 geom_point 后面，theme 前面：
  scale_color_manual(values = c("NBD" = "#2E86AB", "Hinge" = "#E74C3C", "SBD" = "#27AE60"))
  annotate("text", x = 200, y = max(rmsf$RMSF)*1.08, label = "NBD",   color = "#2E86AB", fontface = "bold", size = 3.5) +
  annotate("text", x = 403, y = max(rmsf$RMSF)*1.08, label = "Hinge", color = "#E74C3C", fontface = "bold", size = 3.5) +
  annotate("text", x = 550, y = max(rmsf$RMSF)*1.08, label = "SBD",   color = "#27AE60", fontface = "bold", size = 3.5) +
    # 加在 annotate("text", x = 550...) 后面：
    annotate("point", x = 620, y = max(rmsf$RMSF)*0.95, color = "magenta", size = 2, shape = 21, stroke = 0.8) +
    annotate("text", x = 635, y = max(rmsf$RMSF)*0.95, label = "High flexibility\n(> mean + 1SD)", 
             hjust = 0, size = 2.5, color = "magenta")
  labs(title = "HSPA5 - Residue Flexibility (CABS-flex RMSF)",
       x = "Residue Position", y = "RMSF (Å)") +
  theme_bw() +
  theme(legend.position = c(0.85, 0.85), 
        legend.background = element_rect(fill = "white", color = "grey80"),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 8)) +
  guides(color = guide_legend(title = "Domain"))

ggsave(file.path(base_dir, "HSPA5_CABS_RMSF.pdf"), p, width = 10, height = 5, dpi = 300)
cat("\n✅ 图已保存: HSPA5_CABS_RMSF.pdf\n")

# ---------- 6. 保存 ----------
write_csv(rmsf, file.path(hsp_dir, "HSPA5_CABS_RMSF_processed.csv"))
write_csv(rmsf_summary, file.path(hsp_dir, "HSPA5_CABS_RMSF_summary.csv"))
write_csv(high_flex, file.path(hsp_dir, "HSPA5_CABS_high_flexibility_sites.csv"))
cat("✅ 数据已保存\n")