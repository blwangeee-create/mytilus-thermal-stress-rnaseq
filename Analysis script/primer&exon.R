library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

GENOME_FASTA <- "./genome/GCA_011752425.2_MCOR1.1_genomic.fna"
GTF_FILE     <- "annotation/MCOR1.1.gtf"
TARGET_SYMBOLS <- c("DNAJA4", "DNAJB5", "HSPBP1", "CASP9", "BCL2", "NFKBIA", "BAG3")

# ========== 核心优化：只读需要的染色体 ==========
read_fasta_selective <- function(fasta_file, target_chrs) {
  if (grepl("\\.gz$", fasta_file)) {
    con <- gzfile(fasta_file, "r")
  } else {
    con <- file(fasta_file, "r")
  }
  
  seqs <- list()
  current_name <- NULL
  current_seq <- character()
  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (grepl("^>", line)) {
      # 保存上一个序列（如果它是目标染色体）
      if (!is.null(current_name) && current_name %in% target_chrs) {
        seqs[[current_name]] <- paste(current_seq, collapse = "")
      }
      # 解析新header，只保留染色体名
      current_name <- sub("^>(\\S+).*", "\\1", line)
      current_seq <- character()
    } else {
      current_seq <- c(current_seq, line)
    }
  }
  # 最后一个序列
  if (!is.null(current_name) && current_name %in% target_chrs) {
    seqs[[current_name]] <- paste(current_seq, collapse = "")
  }
  close(con)
  return(seqs)
}

calc_tm <- function(seq) {
  seq <- toupper(seq); nA <- nchar(gsub("[^A]", "", seq)); nT <- nchar(gsub("[^T]", "", seq))
  nG <- nchar(gsub("[^G]", "", seq)); nC <- nchar(gsub("[^C]", "", seq)); len <- nchar(seq)
  if (len < 14) (nA + nT) * 2 + (nG + nC) * 4 else 64.9 + 41 * (nG + nC - 16.4) / len
}
calc_gc <- function(seq) { seq <- toupper(seq); round(nchar(gsub("[^GC]", "", seq)) / nchar(seq) * 100, 1) }

check_primer <- function(seq) {
  seq <- toupper(seq); issues <- c()
  if (grepl("AAAAA|TTTTT|CCCCC|GGGGG", seq)) issues <- c(issues, "polyN")
  if (calc_gc(seq) < 35 || calc_gc(seq) > 65) issues <- c(issues, "GC%")
  if (nchar(seq) < 18 || nchar(seq) > 25) issues <- c(issues, "length")
  paste(issues, collapse = ";")
}

# ========== 重建 id_map ==========
if (!exists("id_map")) {
  header_line <- readLines("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", n = 1)
  col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]
  feat <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", 
                     stringsAsFactors = FALSE, skip = 1, header = FALSE, col.names = col_names)
  mrna_map <- feat %>% dplyr::filter(feature == "mRNA", related_accession != "") %>% 
    dplyr::select(locus_tag, protein_id = related_accession) %>% dplyr::distinct()
  eggnog <- read.delim("annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular", 
                       stringsAsFactors = FALSE, comment.char = "")
  names(eggnog) <- gsub("^#", "", names(eggnog)); names(eggnog) <- gsub("^X\\.", "", names(eggnog))
  symbol_map <- eggnog %>% dplyr::select(query, Preferred_name) %>% 
    dplyr::filter(Preferred_name != "" & Preferred_name != "-") %>% dplyr::distinct()
  id_map <- mrna_map %>% dplyr::left_join(symbol_map, by = c("protein_id" = "query")) %>% 
    dplyr::filter(!is.na(Preferred_name)) %>% 
    dplyr::mutate(gene_id = paste0("gene-", locus_tag)) %>% 
    dplyr::select(gene_id, symbol = Preferred_name) %>% dplyr::distinct()
}

target_map <- id_map %>% dplyr::filter(symbol %in% TARGET_SYMBOLS) %>% 
  dplyr::group_by(symbol) %>% dplyr::slice(1) %>% dplyr::ungroup()
cat("目标基因:\n"); print(target_map)

# ========== GTF 外显子提取 + 确定所需染色体 ==========
gtf_lines <- readLines(GTF_FILE)
exon_list <- list()
needed_chrs <- character()

for (i in seq_len(nrow(target_map))) {
  gid <- target_map$gene_id[i]; sym <- target_map$symbol[i]
  pattern <- paste0('transcript_id "', gid, '"')
  exon_lines <- gtf_lines[grepl(pattern, gtf_lines) & grepl("\texon\t", gtf_lines)]
  cat(sym, "(", gid, "):", length(exon_lines), "个外显子\n")
  if (length(exon_lines) == 0) next
  
  exons <- do.call(rbind, lapply(exon_lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    data.frame(chr = parts[1], start = as.numeric(parts[4]), end = as.numeric(parts[5]), 
               strand = parts[7], stringsAsFactors = FALSE)
  }))
  exons <- exons[order(exons$start), ]
  if (exons$strand[1] == "-") exons <- exons[order(exons$start, decreasing = TRUE), ]
  exon_list[[gid]] <- exons
  needed_chrs <- c(needed_chrs, exons$chr[1])
}

needed_chrs <- unique(needed_chrs)
cat("\n7个基因分布在", length(needed_chrs), "条染色体上:\n")
print(needed_chrs)

# ========== 只读需要的染色体（秒级完成） ==========
cat("\n读取基因组... ")
genome_seqs <- read_fasta_selective(GENOME_FASTA, needed_chrs)
cat(length(genome_seqs), "条序列（仅目标染色体）\n")

# ========== 放宽版引物设计 ==========
design_primers_cDNA_relaxed <- function(gene_id, symbol, exons, genome_seqs, n_candidates = 3) {
  if (nrow(exons) < 2) return(NULL)
  
  chr <- exons$chr[1]
  if (!(chr %in% names(genome_seqs))) return(NULL)
  seq_chr <- genome_seqs[[chr]]
  
  exon_seqs <- sapply(seq_len(nrow(exons)), function(i) substr(seq_chr, exons$start[i], exons$end[i]))
  if (exons$strand[1] == "-") {
    exon_seqs <- sapply(exon_seqs, function(s) chartr("ATCGatcg", "TAGCtagc", paste(rev(strsplit(s, "")[[1]]), collapse = "")))
  }
  
  cDNA <- paste(exon_seqs, collapse = "")
  exon_lengths <- nchar(exon_seqs)
  cum_lengths <- cumsum(exon_lengths)
  junctions <- cum_lengths[-length(cum_lengths)]
  
  candidates <- data.frame()
  
  for (j_idx in seq_along(junctions)) {
    junc <- junctions[j_idx]
    upstream_seq <- substr(cDNA, max(1, junc - 499), junc)
    downstream_seq <- substr(cDNA, junc + 1, min(nchar(cDNA), junc + 500))
    
    for (f_len in 22:18) {
      if (nchar(upstream_seq) < f_len) next
      for (f_pos in (nchar(upstream_seq) - f_len + 1):1) {
        f_seq <- substr(upstream_seq, f_pos, f_pos + f_len - 1)
        if (check_primer(f_seq) != "") next
        f_tm <- calc_tm(f_seq); f_gc <- calc_gc(f_seq)
        if (f_tm < 55 || f_tm > 65) next
        
        for (r_len in 22:18) {
          if (nchar(downstream_seq) < r_len) next
          for (r_pos in 1:(nchar(downstream_seq) - r_len + 1)) {
            r_seq <- substr(downstream_seq, r_pos, r_pos + r_len - 1)
            r_seq_rc <- chartr("ATCGatcg", "TAGCtagc", paste(rev(strsplit(r_seq, "")[[1]]), collapse = ""))
            if (check_primer(r_seq_rc) != "") next
            r_tm <- calc_tm(r_seq_rc); r_gc <- calc_gc(r_seq_rc)
            if (r_tm < 55 || r_tm > 65) next
            
            f_start <- max(1, junc - 499) + f_pos - 1
            r_end <- (junc + 1) + r_pos + r_len - 2
            product_len <- r_end - f_start + 1
            if (product_len < 60 || product_len > 250) next
            
            f_3 <- substr(f_seq, nchar(f_seq)-2, nchar(f_seq))
            r_3 <- substr(r_seq_rc, nchar(r_seq_rc)-2, nchar(r_seq_rc))
            if (grepl(chartr("ATCG", "TAGC", f_3), r_3)) next
            
            f_3at <- substr(f_seq, nchar(f_seq), nchar(f_seq)) %in% c("A", "T")
            r_3at <- substr(r_seq_rc, nchar(r_seq_rc), nchar(r_seq_rc)) %in% c("A", "T")
            score <- as.integer(f_3at) + as.integer(r_3at)
            
            candidates <- rbind(candidates, data.frame(
              symbol = symbol, gene_id = gene_id, 
              junction = paste0("E", j_idx, "-E", j_idx+1),
              F_primer = f_seq, F_Tm = round(f_tm, 1), F_GC = f_gc,
              R_primer = r_seq_rc, R_Tm = round(r_tm, 1), R_GC = r_gc,
              product_bp = product_len, score = score,
              stringsAsFactors = FALSE
            ))
            if (nrow(candidates) >= n_candidates * length(junctions)) break
          }
          if (nrow(candidates) >= n_candidates * length(junctions)) break
        }
        if (nrow(candidates) >= n_candidates * length(junctions)) break
      }
      if (nrow(candidates) >= n_candidates * length(junctions)) break
    }
  }
  
  if (nrow(candidates) == 0) return(NULL)
  candidates <- candidates[order(candidates$score), ]
  return(head(candidates, n_candidates))
}

# ========== 批量设计 ==========
all_primers <- data.frame()
for (i in seq_len(nrow(target_map))) {
  sym <- target_map$symbol[i]; gid <- target_map$gene_id[i]
  if (!(gid %in% names(exon_list))) { cat(sym, ": 无外显子\n"); next }
  cat("\n设计", sym, "... ")
  res <- design_primers_cDNA_relaxed(gid, sym, exon_list[[gid]], genome_seqs, n_candidates = 3)
  if (!is.null(res) && nrow(res) > 0) {
    all_primers <- rbind(all_primers, res)
    cat(nrow(res), "对（score越低越好）\n")
    print(res[, c("junction", "F_primer", "F_Tm", "R_primer", "R_Tm", "product_bp", "score")])
  } else { cat("无合适引物\n") }
}

# ========== 输出 ==========
cat("\n\n========== 候选引物汇总 ==========\n")
if (nrow(all_primers) > 0) {
  print(all_primers[, c("symbol", "junction", "F_primer", "F_Tm", "R_primer", "R_Tm", "product_bp", "score")])
  write.csv(all_primers, "results/RNAseq_analysis/Results/qPCR_primers_relaxed.csv", row.names = FALSE)
  cat("\n已保存\n")
  
  cat("\n========== 每基因最优推荐（score=0 优先）==========\n")
  best <- all_primers %>% dplyr::group_by(symbol) %>% dplyr::slice_min(score, n = 1, with_ties = FALSE)
  print(best[, c("symbol", "junction", "F_primer", "F_Tm", "R_primer", "R_Tm", "product_bp", "score")])
} else {
  cat("仍未找到任何引物。\n")
}





library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

GENOME_FASTA <- "./genome/GCA_011752425.2_MCOR1.1_genomic.fna"
GTF_FILE <- "annotation/MCOR1.1.gtf"

# ========== 1. 重建 id_map ==========
header_line <- readLines("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", n = 1)
col_names <- strsplit(sub("^# ", "", header_line), "\t")[[1]]
feat <- read.delim("annotation/GCA_011752425.2_MCOR1.1_feature_table.txt", 
                   stringsAsFactors = FALSE, skip = 1, header = FALSE, col.names = col_names)
mrna_map <- feat %>% dplyr::filter(feature == "mRNA", related_accession != "") %>% 
  dplyr::select(locus_tag, protein_id = related_accession) %>% dplyr::distinct()

eggnog <- read.delim("annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular", 
                     stringsAsFactors = FALSE, comment.char = "")
names(eggnog) <- gsub("^#", "", names(eggnog)); names(eggnog) <- gsub("^X\\.", "", names(eggnog))
symbol_map <- eggnog %>% dplyr::select(query, Preferred_name) %>% 
  dplyr::filter(Preferred_name != "" & Preferred_name != "-") %>% dplyr::distinct()

id_map <- mrna_map %>% dplyr::left_join(symbol_map, by = c("protein_id" = "query")) %>% 
  dplyr::filter(!is.na(Preferred_name)) %>% 
  dplyr::mutate(gene_id = paste0("gene-", locus_tag)) %>% 
  dplyr::select(gene_id, symbol = Preferred_name) %>% dplyr::distinct()

# ========== 2. 目标基因 ==========
TARGET_SYMBOLS <- c("BAG3", "BCL2", "CASP9", "DNAJA4", "DNAJB5", "HSPBP1", "NFKBIA")
target_map <- id_map %>% dplyr::filter(symbol %in% TARGET_SYMBOLS) %>% 
  dplyr::group_by(symbol) %>% dplyr::slice(1) %>% dplyr::ungroup()
cat("目标基因:\n"); print(target_map)

# ========== 3. 读取基因组（如果内存中没有） ==========
if (exists("genome_seqs")) {
  cat("使用内存中已有的基因组，", length(genome_seqs), "条序列\n")
} else {
  cat("读取基因组... ")
  read_fasta_basic <- function(fasta_file) {
    if (grepl("\\.gz$", fasta_file)) { con <- gzfile(fasta_file, "r") } else { con <- file(fasta_file, "r") }
    lines <- readLines(con); close(con)
    seqs <- list(); current_name <- NULL; current_seq <- character()
    for (line in lines) {
      if (grepl("^>", line)) {
        if (!is.null(current_name)) seqs[[current_name]] <- paste(current_seq, collapse = "")
        current_name <- sub("^>(\\S+).*", "\\1", line)
        current_seq <- character()
      } else { current_seq <- c(current_seq, line) }
    }
    if (!is.null(current_name)) seqs[[current_name]] <- paste(current_seq, collapse = "")
    return(seqs)
  }
  genome_seqs <- read_fasta_basic(GENOME_FASTA)
  cat(length(genome_seqs), "条序列\n")
}

# ========== 4. 从 GTF 提取外显子并拼接 cDNA ==========
gtf_lines <- readLines(GTF_FILE)
cdna_list <- list()

for (i in seq_len(nrow(target_map))) {
  gid <- target_map$gene_id[i]
  sym <- target_map$symbol[i]
  
  pattern <- paste0('transcript_id "', gid, '"')
  exon_lines <- gtf_lines[grepl(pattern, gtf_lines) & grepl("\texon\t", gtf_lines)]
  
  if (length(exon_lines) == 0) next
  
  exons <- do.call(rbind, lapply(exon_lines, function(line) {
    parts <- strsplit(line, "\t")[[1]]
    data.frame(chr = parts[1], start = as.numeric(parts[4]), end = as.numeric(parts[5]), 
               strand = parts[7], stringsAsFactors = FALSE)
  }))
  
  exons <- exons[order(exons$start), ]
  if (exons$strand[1] == "-") exons <- exons[order(exons$start, decreasing = TRUE), ]
  
  chr <- exons$chr[1]
  if (!(chr %in% names(genome_seqs))) next
  seq_chr <- genome_seqs[[chr]]
  
  exon_seqs <- sapply(seq_len(nrow(exons)), function(j) substr(seq_chr, exons$start[j], exons$end[j]))
  if (exons$strand[1] == "-") {
    exon_seqs <- sapply(exon_seqs, function(s) chartr("ATCGatcg", "TAGCtagc", paste(rev(strsplit(s, "")[[1]]), collapse = "")))
  }
  
  full_cdna <- paste(exon_seqs, collapse = "")
  cdna_list[[sym]] <- full_cdna
  
  cat(sym, ": cDNA 长度 =", nchar(full_cdna), "bp (", nrow(exons), "个外显子 )\n")
}

# ========== 5. 保存为 FASTA ==========
dir.create("results/RNAseq_analysis/Results/CDS_sequences", showWarnings = FALSE, recursive = TRUE)

dna_fasta_lines <- character()
for (s in names(cdna_list)) {
  dna_fasta_lines <- c(dna_fasta_lines, 
                       paste0(">", s, " | gene_id=", target_map$gene_id[target_map$symbol == s], " | cDNA"),
                       cdna_list[[s]])
}
writeLines(dna_fasta_lines, "results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA.fasta")
cat("\n✅ 已保存 cDNA 序列到 target_genes_cDNA.fasta\n")

# 汇总
summary_df <- data.frame(
  symbol = names(cdna_list),
  cDNA_bp = sapply(cdna_list, nchar),
  stringsAsFactors = FALSE
)
print(summary_df)
write.csv(summary_df, "results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA_summary.csv", row.names = FALSE)




setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 标准密码子表
CODON_TABLE <- list(
  TTT="F", TTC="F", TTA="L", TTG="L",
  CTT="L", CTC="L", CTA="L", CTG="L",
  ATT="I", ATC="I", ATA="I", ATG="M",
  GTT="V", GTC="V", GTA="V", GTG="V",
  TCT="S", TCC="S", TCA="S", TCG="S",
  CCT="P", CCC="P", CCA="P", CCG="P",
  ACT="T", ACC="T", ACA="T", ACG="T",
  GCT="A", GCC="A", GCA="A", GCG="A",
  TAT="Y", TAC="Y", TAA="*", TAG="*",
  CAT="H", CAC="H", CAA="Q", CAG="Q",
  AAT="N", AAC="N", AAA="K", AAG="K",
  GAT="D", GAC="D", GAA="E", GAG="E",
  TGT="C", TGC="C", TGA="*", TGG="W",
  CGT="R", CGC="R", CGA="R", CGG="R",
  AGT="S", AGC="S", AGA="R", AGG="R",
  GGT="G", GGC="G", GGA="G", GGG="G"
)

translate_orf <- function(dna) {
  dna <- toupper(dna)
  # 找第一个 ATG
  atg_pos <- regexpr("ATG", dna)
  if (atg_pos < 0) return(NULL)
  
  orf <- substr(dna, atg_pos, nchar(dna))
  # 截断到 3 的倍数
  len <- nchar(orf) - (nchar(orf) %% 3)
  orf <- substr(orf, 1, len)
  
  # 翻译
  codons <- sapply(seq(1, len, 3), function(i) substr(orf, i, i+2))
  aa <- sapply(codons, function(c) ifelse(is.null(CODON_TABLE[[c]]), "X", CODON_TABLE[[c]]))
  prot <- paste(aa, collapse = "")
  
  # 截断到第一个 STOP
  stop_pos <- regexpr("[*]", prot)
  if (stop_pos > 0) prot <- substr(prot, 1, stop_pos - 1)
  return(prot)
}

# 读取 cDNA FASTA
lines <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA.fasta")
prot_lines <- character()
current_header <- NULL
current_seq <- character()

for (line in lines) {
  if (grepl("^>", line)) {
    if (!is.null(current_header)) {
      dna <- paste(current_seq, collapse = "")
      prot <- translate_orf(dna)
      if (!is.null(prot)) {
        # 修改 header：把 cDNA 改成 protein，加上长度
        new_header <- gsub("\\| cDNA", paste0("| protein | length=", nchar(prot), "aa"), current_header)
        prot_lines <- c(prot_lines, new_header, prot)
        cat(gsub("^>", "", current_header), "-> 蛋白长度:", nchar(prot), "aa\n")
      } else {
        cat(gsub("^>", "", current_header), "-> 未找到 ATG\n")
      }
    }
    current_header <- line
    current_seq <- character()
  } else {
    current_seq <- c(current_seq, line)
  }
}

# 最后一个
if (!is.null(current_header)) {
  dna <- paste(current_seq, collapse = "")
  prot <- translate_orf(dna)
  if (!is.null(prot)) {
    new_header <- gsub("\\| cDNA", paste0("| protein | length=", nchar(prot), "aa"), current_header)
    prot_lines <- c(prot_lines, new_header, prot)
    cat(gsub("^>", "", current_header), "-> 蛋白长度:", nchar(prot), "aa\n")
  }
}

# 保存
writeLines(prot_lines, "results/RNAseq_analysis/Results/CDS_sequences/target_genes_protein.fasta")
cat("\n✅ 已保存到 target_genes_protein.fasta\n")





setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取 FASTA
lines <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA.fasta")

records <- list()
current_name <- NULL
current_seq <- character()

for (line in lines) {
  if (grepl("^>", line)) {
    if (!is.null(current_name)) {
      records[[current_name]] <- paste(current_seq, collapse = "")
    }
    current_name <- sub("^>(\\S+).*", "\\1", line)
    current_seq <- character()
  } else {
    current_seq <- c(current_seq, line)
  }
}
if (!is.null(current_name)) {
  records[[current_name]] <- paste(current_seq, collapse = "")
}

# 生成表格
df <- data.frame(
  gene = names(records),
  length_bp = nchar(unlist(records)),
  sequence = unlist(records),
  stringsAsFactors = FALSE
)

# 保存
write.csv(df, "results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA_with_sequence.csv", row.names = FALSE)
cat("✅ 已保存到 target_genes_cDNA_with_sequence.csv\n")

# 直接在控制台打印（方便你复制）
cat("\n========== 基因序列汇总 ==========\n")
for (i in seq_len(nrow(df))) {
  cat("\n>", df$gene[i], " | 长度:", df$length_bp[i], "bp\n")
  cat(df$sequence[i], "\n")
}




library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 1. DNAJB5 外显子边界确认（用于手动设计）==========
gtf_lines <- readLines("annotation/MCOR1.1.gtf")
pattern <- paste0('transcript_id "gene-MCOR_7799"')
exon_lines <- gtf_lines[grepl(pattern, gtf_lines) & grepl("\texon\t", gtf_lines)]

exons <- do.call(rbind, lapply(exon_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  data.frame(chr = parts[1], start = as.numeric(parts[4]), end = as.numeric(parts[5]), 
             strand = parts[7], length = as.numeric(parts[5]) - as.numeric(parts[4]) + 1,
             stringsAsFactors = FALSE)
}))
exons <- exons[order(exons$start), ]
if (exons$strand[1] == "-") exons <- exons[order(exons$start, decreasing = TRUE), ]

cat("========== DNAJB5 外显子信息 ==========\n")
print(exons)
cat("\nE1 长度:", exons$length[1], "bp\n")
cat("E2 长度:", exons$length[2], "bp\n")
cat("cDNA 总长度:", sum(exons$length), "bp\n")
cat("E1-E2 边界位置:", exons$length[1], "|", exons$length[1]+1, "\n\n")

# 读取 DNAJB5 cDNA（从之前生成的 FASTA）
cdna_lines <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA.fasta")
dnaJB5_cdna <- NULL
for (i in seq_along(cdna_lines)) {
  if (grepl("^>DNAJB5", cdna_lines[i])) {
    dnaJB5_cdna <- cdna_lines[i+1]
    break
  }
}
cat("DNAJB5 cDNA 序列（前 100 bp + 边界区域 + 后 100 bp）:\n")
cat("5'端:", substr(dnaJB5_cdna, 1, 100), "\n")
cat("边界:", substr(dnaJB5_cdna, max(1, exons$length[1]-50), min(nchar(dnaJB5_cdna), exons$length[1]+50)), "\n")
cat("3'端:", substr(dnaJB5_cdna, nchar(dnaJB5_cdna)-99, nchar(dnaJB5_cdna)), "\n\n")

# ========== 2. 排查 9 个新基因在 eggNOG 中的匹配 ==========
NEW_SYMBOLS <- c("CALR", "SMP", "Pif", "BAX", "BCL2", "Caspase3", "CHS", "CHI", "AMPK-alpha")

# 读取 eggNOG
eggnog <- read.delim("annotation/Galaxy4-[eggNOG Mapper on dataset 2_ annotations].tabular", 
                     stringsAsFactors = FALSE, comment.char = "")
names(eggnog) <- gsub("^#", "", names(eggnog)); names(eggnog) <- gsub("^X\\.", "", names(eggnog))

cat("========== 新基因 eggNOG 匹配排查 ==========\n")
for (s in NEW_SYMBOLS) {
  # 精确匹配 Preferred_name
  exact <- eggnog %>% dplyr::filter(Preferred_name == s) %>% dplyr::select(query, Preferred_name, Description) %>% dplyr::distinct()
  # 模糊匹配（Description 和 Preferred_name）
  fuzzy <- eggnog %>% dplyr::filter(grepl(s, Preferred_name, fixed = TRUE) | grepl(s, Description, fixed = TRUE)) %>% 
    dplyr::select(query, Preferred_name, Description) %>% dplyr::distinct()
  
  cat("\n---", s, "---\n")
  if (nrow(exact) > 0) {
    cat("✅ 精确匹配到", nrow(exact), "条:\n")
    print(head(exact, 3))
  } else if (nrow(fuzzy) > 0) {
    cat("⚠️ 无精确匹配，模糊匹配到", nrow(fuzzy), "条（前3）:\n")
    print(head(fuzzy, 3))
  } else {
    cat("❌ 未找到任何匹配\n")
  }
}



setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取 DNAJB5 cDNA
lines <- readLines("results/RNAseq_analysis/Results/CDS_sequences/target_genes_cDNA.fasta")
dnaJB5 <- NULL
for (i in seq_along(lines)) {
  if (grepl("^>DNAJB5", lines[i])) {
    dnaJB5 <- toupper(lines[i+1])
    break
  }
}

e1_len <- 825  # E1-E2 边界

calc_tm <- function(seq) {
  seq <- toupper(seq); nA <- nchar(gsub("[^A]", "", seq)); nT <- nchar(gsub("[^T]", "", seq))
  nG <- nchar(gsub("[^G]", "", seq)); nC <- nchar(gsub("[^C]", "", seq)); len <- nchar(seq)
  if (len < 14) (nA + nT) * 2 + (nG + nC) * 4 else 64.9 + 41 * (nG + nC - 16.4) / len
}
calc_gc <- function(seq) { seq <- toupper(seq); round(nchar(gsub("[^GC]", "", seq)) / nchar(seq) * 100, 1) }

# 扩大搜索：E1 深入 400bp，E2 全长
up <- substr(dnaJB5, max(1, e1_len-399), e1_len)    # 400 bp
down <- substr(dnaJB5, e1_len+1, nchar(dnaJB5))      # 234 bp

cat("DNAJB5 搜索区域: 上游", nchar(up), "bp, 下游", nchar(down), "bp\n\n")

candidates <- data.frame()

for (f_len in 25:17) {
  if (nchar(up) < f_len) next
  for (f_pos in (nchar(up) - f_len + 1):1) {
    f_seq <- substr(up, f_pos, f_pos + f_len - 1)
    f_tm <- calc_tm(f_seq); f_gc <- calc_gc(f_seq)
    if (f_tm < 50 || f_tm > 70) next
    if (f_gc < 25 || f_gc > 75) next
    if (grepl("AAAAA|TTTTT|CCCCC|GGGGG", f_seq)) next
    
    for (r_len in 25:17) {
      if (nchar(down) < r_len) next
      for (r_pos in 1:(nchar(down) - r_len + 1)) {
        r_seq <- substr(down, r_pos, r_pos + r_len - 1)
        r_seq_rc <- chartr("ATCGatcg", "TAGCtagc", paste(rev(strsplit(r_seq, "")[[1]]), collapse = ""))
        r_tm <- calc_tm(r_seq_rc); r_gc <- calc_gc(r_seq_rc)
        if (r_tm < 50 || f_tm > 70) next
        if (r_gc < 25 || r_gc > 75) next
        if (grepl("AAAAA|TTTTT|CCCCC|GGGGG", r_seq_rc)) next
        
        product_len <- (e1_len + 1 + r_pos + r_len - 2) - (max(1, e1_len-399) + f_pos - 1) + 1
        if (product_len < 80 || product_len > 250) next
        
        candidates <- rbind(candidates, data.frame(
          F_primer = f_seq, F_Tm = round(f_tm,1), F_GC = f_gc,
          R_primer = r_seq_rc, R_Tm = round(r_tm,1), R_GC = r_gc,
          product_bp = product_len,
          stringsAsFactors = FALSE
        ))
        if (nrow(candidates) >= 5) break
      }
      if (nrow(candidates) >= 5) break
    }
    if (nrow(candidates) >= 5) break
  }
  if (nrow(candidates) >= 5) break
}

if (nrow(candidates) > 0) {
  cat("找到", nrow(candidates), "对候选引物:\n")
  print(candidates)
  write.csv(candidates, "results/RNAseq_analysis/Results/DNAJB5_candidates.csv", row.names = FALSE)
} else {
  cat("仍未找到。建议改用普通引物（不跨内含子），在 E1 内部设计。\n")
}












library(tidyverse)
library(Biostrings)
library(rtracklayer)  # ← 加这行

setwd("/Volumes/Expansion/NCBI/PRJNA934294")
out_dir <- "results"

# 重新读取之前保存的稳定性结果
stab <- read.csv(file.path(out_dir, "reference_gene_stability_final.csv"), stringsAsFactors = FALSE)

# 把表达量合并进来（排除 RPS23）
match_df <- data.frame(
  Gene = c("EF1A", "TUBA", "ACTB", "RPS23"),
  Mean_TPM = c(16540.27, 2488.66, 4366.65, 89.19)
)

stab <- stab %>% left_join(match_df, by = "Gene")

# 修正判断：推荐 + 表达量 > 500 + 排除 18S/RPS23
usable <- stab[stab$Judgement == "✅ 推荐" & stab$Mean_TPM > 500 & stab$Gene != "r18S", ]
if (nrow(usable) < 2) {
  usable2 <- stab[stab$Judgement == "⚠️ 可用" & stab$Mean_TPM > 500 & stab$Gene != "r18S", ]
  usable <- rbind(usable, usable2)
}
top2 <- usable$Gene[1:2]
cat("\n>>> 最终推荐双内参：", paste(top2, collapse = " + "), "\n")

# ========== 路径 ==========
gtf_file    <- "annotation/MCOR1.1.gtf"
genome_file <- "genome/GCA_011752425.2_MCOR1.1_genomic.fna"

# 读取基因组
genome <- readDNAStringSet(genome_file)
names(genome) <- sapply(strsplit(names(genome), " "), `[`, 1)

# 读取 GTF
gtf <- import(gtf_file, format = "GTF")

# 提取 CDS 函数
extract_cdna <- function(gene_id, gtf, genome) {
  gr <- gtf[mcols(gtf)$gene_id == gene_id]
  if (length(gr) == 0) {
    alt <- ifelse(grepl("^gene-", gene_id), sub("^gene-", "", gene_id), paste0("gene-", gene_id))
    gr <- gtf[mcols(gtf)$gene_id == alt]
  }
  exons <- gr[mcols(gr)$type %in% c("exon", "CDS")]
  if (length(exons) == 0) return(NULL)
  
  txs <- split(exons, mcols(exons)$transcript_id)
  tx_seqs <- lapply(txs, function(tx) {
    tx <- sort(tx)
    seqs <- getSeq(genome, tx)
    s <- unlist(seqs)
    if (as.character(strand(tx)[1]) == "-") s <- reverseComplement(s)
    s
  })
  lens <- sapply(tx_seqs, length)
  return(tx_seqs[[which.max(lens)]])
}

# 引物设计函数（和目标基因一致）
calc_tm <- function(seq) {
  seq <- toupper(seq)
  nA <- nchar(gsub("[^A]", "", seq)); nT <- nchar(gsub("[^T]", "", seq))
  nG <- nchar(gsub("[^G]", "", seq)); nC <- nchar(gsub("[^C]", "", seq))
  len <- nchar(seq)
  if (len <= 13) (nA+nT)*2 + (nG+nC)*4 else 64.9 + 41*(nG+nC - 16.4)/len
}
calc_gc <- function(seq) {
  seq <- toupper(seq)
  round(nchar(gsub("[^GC]", "", seq)) / nchar(seq) * 100, 1)
}

design_primers <- function(dna, gene_name, product_min = 100, product_max = 200,
                           primer_len = 20, tm_target = 60, tm_tol = 2) {
  seq <- as.character(dna)
  L <- nchar(seq)
  res <- data.frame()
  
  for (f_start in seq(1, L - product_min, by = 10)) {
    for (p_size in seq(product_min, min(product_max, L - f_start + 1), by = 5)) {
      f_end <- f_start + primer_len - 1
      r_start <- f_start + p_size - primer_len
      r_end <- r_start + primer_len - 1
      if (f_end >= r_start) next
      
      f_seq <- substr(seq, f_start, f_end)
      r_seq <- substr(seq, r_start, r_end)
      r_rc <- as.character(reverseComplement(DNAString(r_seq)))
      
      f_tm <- calc_tm(f_seq); r_tm <- calc_tm(r_seq)
      f_gc <- calc_gc(f_seq); r_gc <- calc_gc(r_seq)
      
      ok <- abs(f_tm - r_tm) <= tm_tol &&
        f_tm >= 58 && f_tm <= 62 && r_tm >= 58 && r_tm <= 62 &&
        f_gc >= 40 && f_gc <= 60 && r_gc >= 40 && r_gc <= 60 &&
        !grepl("AAAA|TTTT|GGGG|CCCC", f_seq) &&
        !grepl("AAAA|TTTT|GGGG|CCCC", r_rc)
      
      if (ok) {
        res <- rbind(res, data.frame(
          Gene = gene_name, F_primer = f_seq, R_primer = r_rc,
          Product_size = p_size, F_Tm = round(f_tm,1), R_Tm = round(r_tm,1),
          F_GC = round(f_gc,1), R_GC = round(r_gc,1)
        ))
      }
    }
  }
  res
}

# 为 top2 设计引物
all_primers <- data.frame()
for (g in top2) {
  gid <- stab$Gene_ID[stab$Gene == g]
  cat("\n提取", g, "(", gid, ") ...\n")
  cds <- extract_cdna(gid, gtf, genome)
  if (!is.null(cds)) {
    cat("  cDNA 长度:", length(cds), "bp\n")
    p <- design_primers(cds, g)
    cat("  候选引物:", nrow(p), "对\n")
    if (nrow(p) > 0) {
      p <- p %>% arrange(abs(F_Tm - 60) + abs(R_Tm - 60)) %>% head(3)
      all_primers <- rbind(all_primers, p)
      writeXStringSet(DNAStringSet(setNames(cds, paste0(g, "_CDS"))),
                      file.path(out_dir, paste0(g, "_CDS.fasta")))
    }
  } else {
    cat("  ❌ 未找到 CDS\n")
  }
}

cat("\n=== 推荐引物（", paste(top2, collapse = " + "), "）===\n")
print(all_primers)
write.csv(all_primers, file.path(out_dir, "reference_gene_primers_final.csv"), row.names = FALSE)




library(dplyr)
library(Biostrings)
library(rtracklayer)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# ========== 1. 读取基因组和 GTF ==========
genome <- readDNAStringSet("genome/GCA_011752425.2_MCOR1.1_genomic.fna")
names(genome) <- sapply(strsplit(names(genome), " "), `[`, 1)
gtf <- import("annotation/MCOR1.1.gtf", format = "GTF")

hspa5_ids <- c("gene-MCOR_22078", "gene-MCOR_36097", "gene-MCOR_56583", 
               "gene-MCOR_56587", "gene-MCOR_56591", "gene-MCOR_56598", "gene-MCOR_56599")

# ========== 2. 提取所有 isoform（修复：用 subseq 代替 getSeq）==========
extract_isoforms <- function(gene_ids, gtf, genome) {
  results <- data.frame()
  gids_all <- mcols(gtf)$gene_id
  
  for (gid in gene_ids) {
    idx <- which(!is.na(gids_all) & gids_all == gid)
    gr <- gtf[idx]
    
    if (length(gr) == 0) {
      alt <- ifelse(grepl("^gene-", gid), sub("^gene-", "", gid), paste0("gene-", gid))
      idx <- which(!is.na(gids_all) & gids_all == alt)
      gr <- gtf[idx]
    }
    if (length(gr) == 0) next
    
    txs <- split(gr, mcols(gr)$transcript_id)
    for (tx_name in names(txs)) {
      tx <- txs[[tx_name]]
      exons <- tx[mcols(tx)$type == "exon"]
      if (length(exons) == 0) next
      exons <- sort(exons)
      
      # 修复：用 subseq 直接提取，不依赖 getSeq
      chr_vec <- as.character(seqnames(exons))
      seqs <- lapply(seq_along(exons), function(i) {
        subseq(genome[[chr_vec[i]]], start = start(exons)[i], width = width(exons)[i])
      })
      
      cdna <- unlist(DNAStringSet(seqs))
      if (as.character(strand(exons)[1]) == "-") cdna <- reverseComplement(cdna)
      
      results <- rbind(results, data.frame(
        gene_id = gid, transcript_id = tx_name,
        n_exons = length(exons), cdna_length = length(cdna),
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results)
}

cat("=== HSPA5 所有 isoform 评估 ===\n")
isoforms <- extract_isoforms(hspa5_ids, gtf, genome)
print(isoforms)

if (nrow(isoforms) == 0) stop("未找到任何 HSPA5 isoform")

best <- isoforms[which.max(isoforms$cdna_length), ]
cat("\n>>> 用于引物设计的 isoform（最长）:\n")
print(best)

# ========== 3. 提取该 isoform 的 cDNA 和外显子边界（同样用 subseq）==========
best_gid <- best$gene_id[1]
best_tx  <- best$transcript_id[1]

gids_all <- mcols(gtf)$gene_id
idx <- which(!is.na(gids_all) & gids_all == best_gid)
gr <- gtf[idx]
if (length(gr) == 0) {
  alt <- ifelse(grepl("^gene-", best_gid), sub("^gene-", "", best_gid), paste0("gene-", best_gid))
  idx <- which(!is.na(gids_all) & gids_all == alt)
  gr <- gtf[idx]
}

tx <- gr[mcols(gr)$transcript_id == best_tx]
exons <- tx[mcols(tx)$type == "exon"]
exons <- sort(exons)

chr_vec <- as.character(seqnames(exons))
exon_seqs <- lapply(seq_along(exons), function(i) {
  subseq(genome[[chr_vec[i]]], start = start(exons)[i], width = width(exons)[i])
})

if (as.character(strand(exons)[1]) == "-") {
  exon_seqs <- lapply(exon_seqs, reverseComplement)
}

cdna <- unlist(DNAStringSet(exon_seqs))
cdna_seq <- as.character(cdna)
exon_lengths <- sapply(exon_seqs, length)
cum_lengths <- cumsum(exon_lengths)
junctions <- cum_lengths[-length(cum_lengths)]

cat("\ncDNA 总长度:", length(cdna), "bp | 外显子数:", length(exons), 
    "| Junctions:", paste(junctions, collapse = ", "), "\n")

# ========== 4. 引物设计（与之前完全一致）==========
calc_tm <- function(seq) {
  seq <- toupper(seq)
  nA <- nchar(gsub("[^A]", "", seq)); nT <- nchar(gsub("[^T]", "", seq))
  nG <- nchar(gsub("[^G]", "", seq)); nC <- nchar(gsub("[^C]", "", seq))
  len <- nchar(seq)
  if (len < 14) (nA + nT) * 2 + (nG + nC) * 4 else 64.9 + 41 * (nG + nC - 16.4) / len
}
calc_gc <- function(seq) { 
  seq <- toupper(seq); round(nchar(gsub("[^GC]", "", seq)) / nchar(seq) * 100, 1) 
}
check_primer <- function(seq) {
  seq <- toupper(seq); issues <- c()
  if (grepl("AAAAA|TTTTT|CCCCC|GGGGG", seq)) issues <- c(issues, "polyN")
  if (calc_gc(seq) < 35 || calc_gc(seq) > 65) issues <- c(issues, "GC%")
  if (nchar(seq) < 18 || nchar(seq) > 25) issues <- c(issues, "length")
  paste(issues, collapse = ";")
}

design_primers_relaxed <- function(cdna_seq, junctions, gene_symbol, gene_id, n_candidates = 5) {
  candidates <- data.frame()
  
  for (j_idx in seq_along(junctions)) {
    junc <- junctions[j_idx]
    up   <- substr(cdna_seq, max(1, junc - 499), junc)
    down <- substr(cdna_seq, junc + 1, min(nchar(cdna_seq), junc + 500))
    
    for (f_len in 22:18) {
      if (nchar(up) < f_len) next
      for (f_pos in (nchar(up) - f_len + 1):1) {
        f_seq <- substr(up, f_pos, f_pos + f_len - 1)
        if (check_primer(f_seq) != "") next
        f_tm <- calc_tm(f_seq); f_gc <- calc_gc(f_seq)
        if (f_tm < 55 || f_tm > 65) next
        
        for (r_len in 22:18) {
          if (nchar(down) < r_len) next
          for (r_pos in 1:(nchar(down) - r_len + 1)) {
            r_seq <- substr(down, r_pos, r_pos + r_len - 1)
            r_rc <- chartr("ATCGatcg", "TAGCtagc", paste(rev(strsplit(r_seq, "")[[1]]), collapse = ""))
            if (check_primer(r_rc) != "") next
            r_tm <- calc_tm(r_rc); r_gc <- calc_gc(r_rc)
            if (r_tm < 55 || r_tm > 65) next
            
            f_start <- max(1, junc - 499) + f_pos - 1
            r_end   <- (junc + 1) + r_pos + r_len - 2
            product_len <- r_end - f_start + 1
            if (product_len < 60 || product_len > 250) next
            
            f_3at <- substr(f_seq, nchar(f_seq), nchar(f_seq)) %in% c("A", "T")
            r_3at <- substr(r_rc, nchar(r_rc), nchar(r_rc)) %in% c("A", "T")
            score <- as.integer(f_3at) + as.integer(r_3at)
            
            candidates <- rbind(candidates, data.frame(
              symbol = gene_symbol, gene_id = gene_id,
              junction = paste0("E", j_idx, "-E", j_idx + 1),
              F_primer = f_seq, F_Tm = round(f_tm, 1), F_GC = f_gc,
              R_primer = r_rc, R_Tm = round(r_tm, 1), R_GC = r_gc,
              product_bp = product_len, score = score,
              stringsAsFactors = FALSE
            ))
            if (nrow(candidates) >= n_candidates * length(junctions)) break
          }
          if (nrow(candidates) >= n_candidates * length(junctions)) break
        }
        if (nrow(candidates) >= n_candidates * length(junctions)) break
      }
      if (nrow(candidates) >= n_candidates * length(junctions)) break
    }
  }
  if (nrow(candidates) == 0) return(NULL)
  candidates[order(candidates$score), ]
}

cat("\n=== 设计 HSPA5 跨内含子引物 ===\n")
hspa5_primers <- design_primers_relaxed(cdna_seq, junctions, "HSPA5", best_gid, n_candidates = 5)

if (!is.null(hspa5_primers) && nrow(hspa5_primers) > 0) {
  print(hspa5_primers[, c("symbol", "junction", "F_primer", "F_Tm", "R_primer", "R_Tm", "product_bp", "score")])
  
  out_file <- "results/RNAseq_analysis/Results/qPCR_primers_HSPA5_relaxed.csv"
  write.csv(hspa5_primers, out_file, row.names = FALSE)
  cat("\n✅ 已保存:", out_file, "\n")
  
  cat("\n>>> 最优推荐（score 越低越好）:\n")
  print(hspa5_primers %>% slice_min(score, n = 1, with_ties = FALSE))
} else {
  cat("❌ 未找到合适引物\n")
}