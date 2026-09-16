## Alpha-diversity metrics (Chao1, Shannon, Simpson) across sample groups.

source(here::here("scripts", "01-preload.R"))

diversity_16s <- t(estimateR(t_otu_df))
shannon <- diversity(t_otu_df, index = "shannon")
simpson <- diversity(t_otu_df, index = "simpson")
diversity_16s <- cbind(
  sample_id = rownames(diversity_16s),
  reads_number = rowSums(t_otu_df),
  diversity_16s,
  shannon = shannon,
  simpson = simpson
)
diversity_16s <- as.data.frame(diversity_16s)
write.table(diversity_16s, here(data_out, "16S_diversity.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

numcol <- c(
  "reads_number", "S.obs",
  "S.chao1", "se.chao1",
  "shannon", "simpson"
)

diversity_16s[numcol] <- lapply(diversity_16s[numcol],
                                function(x) (as.numeric(as.character(x))))

diversity_16s <- diversity_16s %>%
  mutate(
    Sample_Type = classify_group(sample_id),
    Ecosystem = classify_ecosystem(sample_id),
    Plant = case_when(
      grepl("S1_rAo|S1_bulk", sample_id) ~ "A.olchonensis",
      grepl("S2_rVo|S2_bulk", sample_id) ~ "V.olchonensis"
    )
  )

diversity_16s$Sample_Type <- factor(
  diversity_16s$Sample_Type, levels = group_levels
)
diversity_16s$Ecosystem <- factor(
  diversity_16s$Ecosystem, levels = ecosystem_levels
)
diversity_16s$Plant <- factor(diversity_16s$Plant,
                              levels = c("A.olchonensis", "V.olchonensis"))

chao1_boxplot <- ggplot(
  diversity_16s, aes(x = Sample_Type, y = S.chao1, fill = Ecosystem)
) +
  geom_boxplot(outlier.shape = 16, alpha = 0.8) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2),
    size = 1.5, alpha = 0.6
  ) +
  stat_summary(
    position = position_dodge(0.75),
    fun = mean, geom = "point", shape = 23, size = 3, fill = "white",
  ) +
  scale_fill_manual(values = ecosystem_colors) +
  labs(
    x = "",
    y = "Chao1 Index (Species Richness)",
    fill = "Sample Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    legend.position = "top"
  ) +
  theme_bw()

ggsave(
  here(images_out, "Chao1_Boxplot.pdf"),
  plot = chao1_boxplot, width = 8, height = 6
)

shannon_boxplot <- ggplot(
  diversity_16s, aes(x = Sample_Type, y = shannon, fill = Ecosystem)
) +
  geom_boxplot(outlier.shape = 16, alpha = 0.8) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2),
    size = 1.5, alpha = 0.6
  ) +
  stat_summary(
    position = position_dodge(0.75),
    fun = mean, geom = "point", shape = 23, size = 3, fill = "white"
  ) +
  scale_fill_manual(values = ecosystem_colors) +
  labs(
    x = "",
    y = "Shannon Index (Diversity)",
    fill = "Sample Type"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    legend.position = "top"
  )

ggsave(
  here(images_out, "Shannon_Boxplot.pdf"),
  plot = shannon_boxplot, width = 8, height = 6
)

cat("\n--- Statistical comparisons for Shannon index (Wilcoxon test) ---\n")
pairwise_shannon <- pairwise.wilcox.test(
  diversity_16s$shannon, diversity_16s$Sample_Type, p.adjust.method = "BH"
)
print(pairwise_shannon)

cat("\n--- Statistical comparisons for Chao1 index (Wilcoxon test) ---\n")
pairwise_chao1 <- pairwise.wilcox.test(
  diversity_16s$S.chao1, diversity_16s$Sample_Type, p.adjust.method = "BH"
)
print(pairwise_chao1)
