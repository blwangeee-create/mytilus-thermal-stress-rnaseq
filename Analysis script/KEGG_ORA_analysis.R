library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)

# 辅助函数：用 GO.db 补全 Description
fix_go_desc <- function(enrich_obj) {
  df <- enrich_obj@result
  missing_idx <- which(is.na(df$Description) | df$Description == "")
  if (length(missing_idx) > 0) {
    # 获取缺失 ID
    miss_ids <- df$ID[missing_idx]
    # 从 GO.db 提取 TERM
    valid <- intersect(miss_ids, keys(GO.db, keytype = "GOID"))
    if (length(valid) > 0) {
      terms <- AnnotationDbi::select(GO.db, keys = valid, columns = "TERM", keytype = "GOID")
      colnames(terms)[1] <- "ID"
      df <- df %>% left_join(terms, by = "ID") %>%
        mutate(Description = ifelse(is.na(TERM), Description, TERM)) %>%
        dplyr::select(-TERM)
    }
  }
  df
}

# 修复两个 ORA 结果的 Description
go_T24_fixed <- fix_go_desc(go_enrich)
go_T32_fixed <- fix_go_desc(go_enrich_T32)

# 取前20条显著通路，并过滤掉修复后仍为 NA 的极少数条目
top20_T24 <- go_T24_fixed %>%
  filter(pvalue < 0.05, !is.na(Description), Description != "") %>%
  arrange(pvalue) %>%
  head(20) %>%
  mutate(group = "24°C vs 13°C")

top20_T32 <- go_T32_fixed %>%
  filter(pvalue < 0.05, !is.na(Description), Description != "") %>%
  arrange(pvalue) %>%
  head(20) %>%
  mutate(group = "32°C vs 13°C")

# 合并并排序
go_plot <- bind_rows(top20_T24, top20_T32) %>%
  mutate(logP = -log10(pvalue))

# 让 y 轴从上到下显著性递减
go_plot <- go_plot %>%
  arrange(desc(logP)) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# 画棒棒糖图
p_go <- ggplot(go_plot, aes(x = logP, y = Description, color = group)) +
  geom_segment(aes(xend = 0, yend = Description), linewidth = 0.8) +
  geom_point(aes(size = Count), alpha = 0.9) +
  scale_color_manual(values = c("24°C vs 13°C" = "#FFD100", "32°C vs 13°C" = "#00A1E4")) +
  scale_size_continuous(range = c(3, 8)) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, max(go_plot$logP) * 1.1)
  ) +
  labs(x = "-log10(p-value)", y = NULL, title = "GO Enrichment: Temperature Effects") +
  theme_bw(base_family = "Arial") +
  theme(
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

ggsave("GO_lollipop_compare_fixed.pdf", p_go, width = 14, height = 10)
ggsave("GO_lollipop_compare_fixed.png", p_go, width = 14, height = 10, dpi = 300)