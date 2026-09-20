############################################################
# RNA-seq Differential Expression Analysis
# Mytilus coruscus temperature stress experiment
# DESeq2 workflow
############################################################


############################
# 1. Load packages
############################

library(DESeq2)
library(dplyr)
library(tidyverse)


############################
# 2. Working directory
############################

base_dir <- "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis"


############################
# 3. Read count matrix
############################


counts_file <- 
  "/Volumes/Expansion/NCBI/PRJNA934294/results/counts/counts_s0.txt"


raw_counts <- read.delim(
  counts_file,
  comment.char="#",
  check.names=FALSE
)


# extract count matrix

count_matrix <- raw_counts[,7:ncol(raw_counts)]


rownames(count_matrix) <- raw_counts$Geneid



############################
# 4. Rename samples
############################


colnames(count_matrix) <- c(
  "T32_3",
  "T32_2",
  "T32_1",
  "T24_3",
  "T24_2",
  "T24_1",
  "Ctrl_3",
  "Ctrl_2",
  "Ctrl_1"
)



############################
# 5. Sample information
############################


sample_info <- data.frame(
  
  row.names = colnames(count_matrix),
  
  condition=c(
    "T32",
    "T32",
    "T32",
    "T24",
    "T24",
    "T24",
    "Ctrl",
    "Ctrl",
    "Ctrl"
  )
  
)


sample_info$condition <- factor(
  sample_info$condition,
  levels=c(
    "Ctrl",
    "T24",
    "T32"
  )
)



############################
# 6. DESeq2 object
############################


dds <- DESeqDataSetFromMatrix(
  
  countData=count_matrix,
  
  colData=sample_info,
  
  design=~condition
  
)



############################
# 7. Filter low expression genes
############################


keep <- rowSums(
  counts(dds)>=10
)>=3


dds <- dds[keep,]



############################
# 8. Differential expression
############################


dds <- DESeq(dds)



############################
# 9. Extract comparisons
############################


res_T24 <- results(
  dds,
  contrast=c(
    "condition",
    "T24",
    "Ctrl"
  )
)


res_T32 <- results(
  dds,
  contrast=c(
    "condition",
    "T32",
    "Ctrl"
  )
)



res_T32vsT24 <- results(
  dds,
  contrast=c(
    "condition",
    "T32",
    "T24"
  )
)



############################
# 10. Export DEG tables
############################


write.csv(
  as.data.frame(res_T24),
  file.path(
    base_dir,
    "Results",
    "DESeq2_T24_vs_Ctrl.csv"
  )
)


write.csv(
  as.data.frame(res_T32),
  file.path(
    base_dir,
    "Results",
    "DESeq2_T32_vs_Ctrl.csv"
  )
)


write.csv(
  as.data.frame(res_T32vsT24),
  file.path(
    base_dir,
    "Results",
    "DESeq2_T32_vs_T24.csv"
  )
)



############################
# 11. Save R objects
############################


save(
  dds,
  res_T24,
  res_T32,
  res_T32vsT24,
  count_matrix,
  sample_info,
  
  file=file.path(
    base_dir,
    "Results",
    "DESeq2_complete_objects.RData"
  )
)


############################################################
# END
############################################################
source(
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Scripts/01_DESeq2/DESeq2_analysis.R"
)
sum(!is.na(res_T24$padj) & res_T24$padj < 0.05)

sum(!is.na(res_T32$padj) & res_T32$padj < 0.05)

sum(!is.na(res_T32vsT24$padj) & res_T32vsT24$padj < 0.05)
get_highconf
head(des_T24_up)
length(des_T24_up)
length(ed_T24_up)
length(lim_up)
ls(pattern="lim")
ls(pattern="up")
head(highconf_T24)
length(highconf_T24)
