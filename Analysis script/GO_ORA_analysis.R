gene_list_T32 <- as.character(highconf_T32$up)

# GO
go_enrich_T32 <- enricher(gene = gene_list_T32,
                          TERM2GENE = go_term2gene,
                          universe = universe_go,
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05)

# KEGG
kegg_enrich_T32 <- enricher(gene = gene_list_T32,
                            TERM2GENE = kegg_term2gene,
                            universe = universe_kegg,
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)

# 添加 GO 名称（如果 GO.db 已加载）
go_enrich_T32@result$Description <- go_term_names$name[match(go_enrich_T32@result$ID, go_term_names$term)]

cat("T32 vs Ctrl GO 条目数:", nrow(go_enrich_T32), "\n")
cat("T32 vs Ctrl KEGG 条目数:", nrow(kegg_enrich_T32), "\n")

# 查看前 15 条
head(go_enrich_T32, 15)


############################################################
# GO over-representation analysis
# RNA-seq High-confidence DEG
# Mytilus coruscus
############################################################


library(enrichit)
library(GO.db)
library(AnnotationDbi)
library(dplyr)


############################################################
# Load DEG list
############################################################


base_dir <-
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis"


deg_T24 <- read.table(
  file.path(
    base_dir,
    "Results/HighConf_DEG_T24vsCtrl.txt"
  ),
  header=FALSE,
  stringsAsFactors=FALSE
)


deg_T32 <- read.table(
  file.path(
    base_dir,
    "Results/HighConf_DEG_T32vsCtrl.txt"
  ),
  header=FALSE,
  stringsAsFactors=FALSE
)


genes_T24 <- deg_T24$V1
genes_T32 <- deg_T32$V1



############################################################
# GO ORA objects
############################################################


# Existing objects saved in Phase1_DEG_analysis.RData
# go_enrich
# go_enrich_T32


load(
  file.path(
    base_dir,
    "Results/Phase1_DEG_analysis.RData"
  )
)


############################################################
# Fix missing GO descriptions
############################################################


fix_go_desc <- function(enrich_obj){
  
  df <- enrich_obj@result
  
  
  missing_idx <- which(
    is.na(df$Description) |
      df$Description==""
  )
  
  
  if(length(missing_idx)>0){
    
    miss_ids <- df$ID[missing_idx]
    
    
    valid <- intersect(
      miss_ids,
      keys(GO.db,keytype="GOID")
    )
    
    
    if(length(valid)>0){
      
      terms <- AnnotationDbi::select(
        GO.db,
        keys=valid,
        columns="TERM",
        keytype="GOID"
      )
      
      
      colnames(terms)[1]="ID"
      
      
      df <- df %>%
        left_join(
          terms,
          by="ID"
        ) %>%
        mutate(
          Description=
            ifelse(
              is.na(TERM),
              Description,
              TERM
            )
        ) %>%
        select(-TERM)
      
    }
    
  }
  
  
  df
  
}


go_T24_fixed <- fix_go_desc(go_enrich)

go_T32_fixed <- fix_go_desc(go_enrich_T32)



############################################################
# Export GO tables
############################################################


write.csv(
  go_T24_fixed,
  file.path(
    base_dir,
    "Results/GO_T24_ORA_result.csv"
  ),
  row.names=FALSE
)


write.csv(
  go_T32_fixed,
  file.path(
    base_dir,
    "Results/GO_T32_ORA_result.csv"
  ),
  row.names=FALSE
)
list.files(
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Scripts/02_GO_KEGG"
)