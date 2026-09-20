library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取注释
gene_annot <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors=FALSE)

# 加载 HighConf
load("results/RNAseq_analysis/Results/Phase1_DEG_analysis.RData")

# 函数：HighConf → symbol 文本
make_string_input <- function(hc_obj, comp_name) {
  
  symbols <- gene_annot %>%
    filter(gene_id %in% hc_obj$all) %>%
    filter(Preferred_name != "" & !is.na(Preferred_name)) %>%
    distinct(Preferred_name) %>%
    pull(Preferred_name)
  
  # 保存：一行一个 symbol，纯文本
  out_file <- paste0("results/PPI/", comp_name, "_STRING_input.txt")
  writeLines(symbols, out_file)
  
  cat(comp_name, ":\n")
  cat("  HighConf genes:", length(hc_obj$all), "\n")
  cat("  Unique symbols:", length(symbols), "\n")
  cat("  Saved to:", out_file, "\n")
  cat("  First 10:", paste(head(symbols, 10), collapse=", "), "\n\n")
  
  return(symbols)
}

# 生成 T32 + T24
sym_T32 <- make_string_input(highconf_T32, "T32")
sym_T24 <- make_string_input(highconf_T24, "T24")

cat("========================================\n")
cat("文件已生成。下一步：\n")
cat("========================================\n")


setwd("/Volumes/Expansion/NCBI/PRJNA934294")

cat("=== T32 ===\n")
t32_edges <- read.delim("results/PPI/T32_string_edges.tsv", stringsAsFactors=FALSE)
cat("Rows:", nrow(t32_edges), "| Cols:", ncol(t32_edges), "\n")
cat("Colnames:", paste(colnames(t32_edges), collapse=", "), "\n")
cat("First 3 rows:\n")
print(head(t32_edges, 3))

cat("\n=== T24 ===\n")
t24_edges <- read.delim("results/PPI/T24_string_edges.tsv", stringsAsFactors=FALSE)
cat("Rows:", nrow(t24_edges), "| Cols:", ncol(t24_edges), "\n")
cat("Colnames:", paste(colnames(t24_edges), collapse=", "), "\n")
cat("First 3 rows:\n")
print(head(t24_edges, 3))


library(dplyr)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

# 读取注释和 DESeq2 结果
gene_annot <- read.csv("results/RNAseq_analysis/Results/gene_annotation_table.csv", stringsAsFactors=FALSE)
deseq_t32 <- read.csv("results/RNAseq_analysis/Results/DESeq2_T32_vs_Ctrl.csv", stringsAsFactors=FALSE)
deseq_t24 <- read.csv("results/RNAseq_analysis/Results/DESeq2_T24_vs_Ctrl.csv", stringsAsFactors=FALSE)

# 读取边文件，提取节点
t32_edges <- read.delim("results/PPI/T32_string_edges.tsv", stringsAsFactors=FALSE)
t24_edges <- read.delim("results/PPI/T24_string_edges.tsv", stringsAsFactors=FALSE)

t32_nodes <- unique(c(t32_edges[[1]], t32_edges[[2]]))
t24_nodes <- unique(c(t24_edges[[1]], t24_edges[[2]]))

# 生成节点属性
make_attrs <- function(nodes, deseq, label) {
  attrs <- gene_annot %>%
    filter(Preferred_name %in% nodes, Preferred_name != "", !is.na(Preferred_name)) %>%
    select(gene_id, Preferred_name) %>%
    left_join(deseq %>% select(gene_id=X, log2FoldChange, padj), by="gene_id") %>%
    group_by(Preferred_name) %>%
    summarise(log2FC = mean(log2FoldChange, na.rm=TRUE), padj = min(padj, na.rm=TRUE), .groups="drop") %>%
    filter(!is.na(log2FC))
  
  write.csv(attrs, paste0("results/PPI/", label, "_node_attrs.csv"), row.names=FALSE)
  cat(label, "nodes in network:", length(nodes), "| with log2FC:", nrow(attrs), "\n")
}

make_attrs(t32_nodes, deseq_t32, "T32")
make_attrs(t24_nodes, deseq_t24, "T24")

cat("\n文件已保存：\n")
cat("  results/PPI/T32_node_attrs.csv\n")
cat("  results/PPI/T24_node_attrs.csv\n")


library(igraph)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

edges_t24 <- read.delim("results/PPI/T24_string_edges.tsv default node.csv", stringsAsFactors=FALSE)
g_t24 <- graph_from_data_frame(edges_t24[, c("X.node1", "node2")], directed=FALSE)
comm_t24 <- cluster_walktrap(g_t24, steps=4)

mod_sizes <- sizes(comm_t24)
top3_ids <- order(mod_sizes, decreasing=TRUE)[1:3]

all_nodes <- V(g_t24)$name
module_label <- rep("Other", length(all_nodes))
for(i in 1:3){
  mod_nodes <- names(membership(comm_t24))[membership(comm_t24) == top3_ids[i]]
  module_label[all_nodes %in% mod_nodes] <- paste0("Module_", i)
}

mod_df <- data.frame(symbol = all_nodes, module = module_label, stringsAsFactors = FALSE)
write.csv(mod_df, "results/PPI/T24_all_nodes_module.csv", row.names = FALSE)
cat("Saved. Total nodes:", nrow(mod_df), "\n")
print(table(mod_df$module))



library(igraph)

setwd("/Volumes/Expansion/NCBI/PRJNA934294")

edges_t24 <- read.delim("results/PPI/T24_string_edges.tsv", stringsAsFactors=FALSE)
g_t24 <- graph_from_data_frame(edges_t24[, c("X.node1", "node2")], directed=FALSE)
comm_t24 <- cluster_walktrap(g_t24, steps=4)

mod_sizes <- sizes(comm_t24)
top3_ids <- order(mod_sizes, decreasing=TRUE)[1:3]

all_nodes <- V(g_t24)$name
module_label <- rep("Other", length(all_nodes))
for(i in 1:3){
  mod_nodes <- names(membership(comm_t24))[membership(comm_t24) == top3_ids[i]]
  module_label[all_nodes %in% mod_nodes] <- paste0("Module_", i)
}

mod_df <- data.frame(symbol = all_nodes, module = module_label, stringsAsFactors = FALSE)
write.csv(mod_df, "results/PPI/T24_all_nodes_module.csv", row.names = FALSE)
cat("Saved. Total nodes:", nrow(mod_df), "\n")
print(table(mod_df$module))