############################################################
# High confidence DEG identification
# Consensus between DESeq2, edgeR and limma
############################################################


library(dplyr)


base_dir <- 
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis"



############################################################
# Function
############################################################


get_highconf <- function(
    des_up,
    des_down,
    ed_up,
    ed_down,
    lim_up,
    lim_down
){
  
  up_all <- c(
    des_up,
    ed_up,
    lim_up
  )
  
  
  down_all <- c(
    des_down,
    ed_down,
    lim_down
  )
  
  
  up_high <- names(
    table(up_all)[table(up_all)>=2]
  )
  
  
  down_high <- names(
    table(down_all)[table(down_all)>=2]
  )
  
  
  list(
    up=up_high,
    down=down_high,
    all=c(
      up_high,
      down_high
    )
  )
}



############################################################
# T24 vs Control
############################################################


highconf_T24 <- get_highconf(
  
  des_T24_up,
  des_T24_down,
  
  ed_T24_up,
  ed_T24_down,
  
  lim_T24_up,
  lim_T24_down
  
)



############################################################
# T32 vs Control
############################################################


highconf_T32 <- get_highconf(
  
  des_T32_up,
  des_T32_down,
  
  ed_T32_up,
  ed_T32_down,
  
  lim_T32_up,
  lim_T32_down
  
)



############################################################
# T32 vs T24
############################################################


highconf_T32v24 <- get_highconf(
  
  des_T32v24_up,
  des_T32v24_down,
  
  ed_T32v24_up,
  ed_T32v24_down,
  
  lim_T32v24_up,
  lim_T32v24_down
  
)



############################################################
# Export
############################################################

write.table(
  highconf_T24$all,
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/HighConf_DEG_T24vsCtrl.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  highconf_T32$all,
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/HighConf_DEG_T32vsCtrl.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  highconf_T32v24$all,
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Results/HighConf_DEG_T32vsT24.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


############################################################
# Check
############################################################


cat(
  "T24 HighConf:",
  length(highconf_T24$all),
  "\n"
)

cat(
  "T32 HighConf:",
  length(highconf_T32$all),
  "\n"
)

cat(
  "T32vsT24 HighConf:",
  length(highconf_T32v24$all),
  "\n"
)
list.files(
  "/Volumes/Expansion/NCBI/PRJNA934294/results/RNAseq_analysis/Scripts/01_DESeq2"
)