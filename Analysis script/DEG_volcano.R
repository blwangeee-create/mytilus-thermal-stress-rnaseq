############################################################
## Volcano plot function
## High-confidence DEG + top25 labels
## Gene name priority:
## Preferred_name > Description > gene ID
############################################################


make_volcano <- function(res, title){
  
  
  df <- as.data.frame(res)
  
  df$gene <- rownames(df)
  
  
  ## remove NA
  df <- df[
    !is.na(df$padj),
  ]
  
  
  ## DEG classification
  
  df$significance <- "NS"
  
  
  df$significance[
    df$padj < 0.05 &
      df$log2FoldChange > 1
  ] <- "Up"
  
  
  df$significance[
    df$padj < 0.05 &
      df$log2FoldChange < -1
  ] <- "Down"
  
  
  
  ########################################################
  ## Select top25 significant genes
  ########################################################
  
  
  top <- df %>%
    
    filter(
      significance!="NS"
    ) %>%
    
    arrange(
      padj
    ) %>%
    
    head(25)
  
  
  
  ########################################################
  ## Add annotation
  ########################################################
  
  
  top <-
    
    left_join(
      
      top,
      
      annotation_unique %>%
        dplyr::select(
          gene_id,
          Preferred_name,
          Description
        ),
      
      by=c(
        "gene"="gene_id"
      )
      
    )
  
  
  
  ########################################################
  ## Generate label
  ########################################################
  
  
  top$label <- NA
  
  
  ## 1. Preferred_name first
  
  top$label[
    !is.na(top$Preferred_name) &
      top$Preferred_name!="-" &
      top$Preferred_name!=""
  ] <-
    top$Preferred_name[
      !is.na(top$Preferred_name) &
        top$Preferred_name!="-" &
        top$Preferred_name!=""
    ]
  
  
  ## 2. Description second
  
  idx <- is.na(top$label) &
    !is.na(top$Description) &
    top$Description!="-" &
    top$Description!=""
  
  
  top$label[idx] <-
    top$Description[idx]
  
  
  
  ## 3. gene ID last
  
  top$label[
    is.na(top$label)
  ] <-
    top$gene[
      is.na(top$label)
    ]
  
  
  
  ########################################################
  ## Label color type
  ########################################################
  
  
  top$label_type <- "Gene"
  
  
  top$label_type[
    top$label ==
      top$Description
  ] <- "Description"
  
  
  
  
  ########################################################
  ## Volcano plot
  ########################################################
  
  
  p <- ggplot(
    
    df,
    
    aes(
      x=log2FoldChange,
      y=-log10(padj),
      color=significance
    )
    
  ) +
    
    
    geom_point(
      size=1.5,
      alpha=0.6
    ) +
    
    
    scale_color_manual(
      values=c(
        "Down"="#0B8BEE",
        "NS"="grey70",
        "Up"="#FFCD00"
      )
    )+
    
    ########################################################
  ## Add labels
  ########################################################
  
  
  ggrepel::geom_text_repel(
    
    data=top,
    
    aes(
      label=label,
      color=label_type
    ),
    
    size=3,
    
    max.overlaps=50,
    
    box.padding=0.5,
    
    point.padding=0.3,
    
    segment.color="grey60",
    
    show.legend=FALSE
    
  ) +
    
    
    scale_color_manual(
      
      values=c(
        
        "Down"="#0B8BEE",
        "NS"="grey70",
        "Up"="#FFCD00",
        "Gene"="black",
        "Description"="grey50"
        
      ),
      
      breaks=c(
        "Down",
        "NS",
        "Up"
      )
      
    ) +
    
    
    theme_classic() +
    
    
    labs(
      
      title=title,
      
      x="log2 Fold Change",
      
      y="-log10 adjusted P value"
      
    )
  
  
  
  return(p)
  
}
p_T24 <- make_volcano(
  res_T24,
  "T24 vs Control"
)


p_T32 <- make_volcano(
  res_T32,
  "T32 vs Control"
)


p_T32v24 <- make_volcano(
  res_T32vsT24,
  "T32 vs T24"
)
p_T24
p_T32
p_T32v24
