## Indicator species (IndVal) analysis at Family level

source(here::here("scripts", "01-preload.R"))

data_count <- normalize_otu(otu_df)
data_family <- aggregate_by_taxon(
  data_count, paste(tax_df$Order, tax_df$Family, sep = "_")
)

otu_table_for_indicator <- as.data.frame(t(data_family))

sample_names_indicator <- rownames(otu_table_for_indicator)
groups <- classify_group(sample_names_indicator)

indval_family <- multipatt(otu_table_for_indicator,
                           cluster = groups,
                           func = "r.g",
                           duleg = FALSE,
                           control = how(nperm = indicator_nperm))

indval_result <- indval_family$sign
indval_result$p.value <- indval_family$sign$p.value
indval_result$Family <- rownames(indval_result)

signif_indval_result <- indval_result %>%
  filter(p.value <= p_value_threshold) %>%
  arrange(p.value)

signif_indval_result <- signif_indval_result %>%
  mutate(Indicator_Group = case_when(
    s.S1_bulk == 1 & s.S1_rAo == 0 & s.S2_bulk == 0 & s.S2_rVo == 0 ~ "S1_bulk",
    s.S1_bulk == 0 & s.S1_rAo == 1 & s.S2_bulk == 0 & s.S2_rVo == 0 ~ "S1_rAo",
    s.S1_bulk == 0 & s.S1_rAo == 0 & s.S2_bulk == 1 & s.S2_rVo == 0 ~ "S2_bulk",
    s.S1_bulk == 0 & s.S1_rAo == 0 & s.S2_bulk == 0 & s.S2_rVo == 1 ~ "S2_rVo",
    s.S1_bulk == 1 & s.S1_rAo == 1 & s.S2_bulk == 0 & s.S2_rVo == 0 ~ "S1_bulk+S1_rAo",  # nolint
    s.S1_bulk == 0 & s.S1_rAo == 0 & s.S2_bulk == 1 & s.S2_rVo == 1 ~ "S2_bulk+S2_rVo",  # nolint
    TRUE ~ "Other"
  ))

strong_indicators <- signif_indval_result %>%
  filter(stat > indicator_stat_threshold) %>%
  arrange(desc(stat))

strong_barplot <- ggplot(
  strong_indicators %>% head(20),
  aes(x = stat, y = reorder(Family, stat), fill = Indicator_Group)
) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Indicator group") +
  labs(title = "Most impactful indicator families (stat > 0.5)",
       subtitle = "Top 20 indicator taxa with largest association",
       x = "Value of indicator status (stat)",
       y = "Bacterial Family") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")

ggsave(
  here(images_out, "Barplot_Strong_Indicators.pdf"),
  strong_barplot, width = 10, height = 8
)

heatmap_families <- strong_indicators %>%
  head(indicator_heatmap_top_n) %>%
  pull(Family)

heatmap_data <- otu_table_for_indicator[, heatmap_families, drop = FALSE]
heatmap_data_scaled <- as.data.frame(scale(heatmap_data))
heatmap_data_scaled$Group <- groups

heatmap_data_avg <- heatmap_data_scaled %>%
  group_by(Group) %>%
  summarise(across(all_of(heatmap_families), mean)) %>%
  tibble::column_to_rownames("Group")

heatmap_data_long <- heatmap_data_avg %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Group") %>%
  pivot_longer(cols = -Group, names_to = "Family", values_to = "Z_score")

heatmap_plot <- ggplot(
  heatmap_data_long, aes(x = Group, y = Family, fill = Z_score)
) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Z-score\n(relative abundance)") +
  labs(x = "Samples Group", y = "Family") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))

ggsave(
  here(images_out, "Heatmap_Strong_Indicators.pdf"),
  heatmap_plot, width = 10, height = 8
)

write.csv(
  signif_indval_result,
  here(data_out, "All_Indicator_Results.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  strong_indicators,
  here(data_out, "Strong_Indicators_Results.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

cat("Analysis done.\n")
cat(paste(
  "Total count of indicator species (p <", p_value_threshold, "):",
  nrow(signif_indval_result),
  "\n"
))
cat(paste(
  "Count of strong indicators (stat >", indicator_stat_threshold, "):",
  nrow(strong_indicators),
  "\n"
))
cat("Distribution by group:\n")
print(table(strong_indicators$Indicator_Group))

## ---- Boxplots for the top indicator families of each group ----

top_indicators_list <- list()

for (group in group_levels) {
  group_indicators <- strong_indicators %>%
    filter(Indicator_Group == group) %>%
    arrange(desc(stat)) %>%
    head(indicator_top_n_per_group)

  if (nrow(group_indicators) > 0) {
    top_indicators_list[[group]] <- group_indicators$Family
  }
}

all_top_families <- unique(unlist(top_indicators_list))
cat(paste(
  "\nTotal unique top indicators for plots:",
  length(all_top_families),
  "\n"
))

plot_data <- as.data.frame(otu_table_for_indicator) * 100
plot_data <- plot_data[, all_top_families, drop = FALSE]
plot_data$Sample_Group <- groups
plot_data$Sample_ID <- rownames(plot_data)

plot_data_long <- plot_data %>%
  pivot_longer(cols = -c(Sample_Group, Sample_ID),
               names_to = "Family",
               values_to = "Relative_Abundance")

create_family_boxplot <- function(family_name) {
  family_data <- plot_data_long %>%
    filter(Family == family_name)

  comparisons_list <- list(c("S1_bulk", "S1_rAo"),
                           c("S2_bulk", "S2_rVo"),
                           c("S1_rAo", "S2_rVo"))

  ggplot(
    family_data,
    aes(x = Sample_Group, y = Relative_Abundance, fill = Sample_Group)
  ) +
    geom_boxplot(
      outlier.shape = 16,
      outlier.size = 1.5,
      alpha = 0.8,
      width = 0.6
    ) +
    geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
    stat_compare_means(comparisons = comparisons_list,
                       method = "wilcox.test",
                       label = "p.signif",
                       hide.ns = TRUE,
                       tip.length = 0.01,
                       exact = FALSE) +
    scale_fill_manual(values = c("S1_bulk" = "#ff9999",
                                 "S1_rAo" = "#99ccff",
                                 "S2_bulk" = "#99ff99",
                                 "S2_rVo" = "#cc99ff")) +
    labs(title = family_name,
         x = "",
         y = "Rel. abundance (%)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.title.y = element_text(size = 9),
          legend.position = "none",
          plot.title = element_text(face = "italic", size = 10),
          panel.grid.minor = element_blank())
}

create_row_with_title <- function(plot_list, title_text) {
  row_plot <- ggarrange(plotlist = plot_list, ncol = length(plot_list))
  annotate_figure(row_plot,
                  top = text_grob(title_text,
                                  face = "bold",
                                  size = 12,
                                  hjust = 0,
                                  x = 0))
}

group_row_titles <- c(
  S1_bulk = "A. S1_bulk (Psammozem soil) indicator families",
  S1_rAo  = "B. S1_rAo (A. olchonensis rhizosphere) indicator families",
  S2_bulk = "C. S2_bulk (Litozem soil) indicator families",
  S2_rVo  = "D. S2_rVo (V. olchonensis rhizosphere) indicator families"
)

group_rows <- lapply(group_levels, function(group) {
  families <- top_indicators_list[[group]]
  if (length(families) == 0) return(NULL)
  plots <- lapply(setNames(families, families), create_family_boxplot)
  create_row_with_title(plots, group_row_titles[[group]])
})
group_rows <- Filter(Negate(is.null), group_rows)

combined_plot <- ggarrange(plotlist = group_rows, nrow = length(group_rows), ncol = 1)

ggsave(
  here(images_out, "Indicator_Families_Boxplots.pdf"),
  combined_plot, width = 12, height = 14
)

summary_stats <- plot_data_long %>%
  group_by(Family, Sample_Group) %>%
  summarise(
    Mean = round(mean(Relative_Abundance, na.rm = TRUE), 3),
    Median = round(median(Relative_Abundance, na.rm = TRUE), 3),
    SD = round(sd(Relative_Abundance, na.rm = TRUE), 3),
    Min = round(min(Relative_Abundance, na.rm = TRUE), 3),
    Max = round(max(Relative_Abundance, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(Family, Sample_Group)

write.csv(
  summary_stats,
  here(data_out, "Top_Indicators_Boxplot_Statistics.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

cat("TOP INDICATORS FOR PLOTS:\n")
for (group in group_levels) {
  if (!is.null(top_indicators_list[[group]])) {
    cat(paste(
      group, ":", paste(top_indicators_list[[group]], collapse = ", "), "\n"
    ))
  }
}

## ---- Scaled (0-1 per family) heatmap of the top indicator families ----

top_families_by_mean <- summary_stats %>%
  filter(Sample_Group %in% group_levels) %>%
  group_by(Family) %>%
  summarise(max_mean = max(Mean)) %>%
  arrange(desc(max_mean)) %>%
  head(indicator_heatmap_top_n)

stats_top <- summary_stats %>% filter(Family %in% top_families_by_mean$Family)

mean_matrix <- stats_top %>%
  select(Family, Sample_Group, Mean) %>%
  pivot_wider(names_from = Sample_Group, values_from = Mean, values_fill = 0)

mean_matrix_mat <- as.matrix(mean_matrix[, -1])
rownames(mean_matrix_mat) <- mean_matrix$Family

order_families <- order(apply(mean_matrix_mat, 1, max), decreasing = TRUE)
mean_matrix_sorted <- mean_matrix_mat[order_families, ]

mean_matrix_scaled <- t(apply(mean_matrix_sorted, 1, function(x) {
  if (max(x) == min(x)) return(rep(0, length(x)))
  (x - min(x)) / (max(x) - min(x))
}))

mean_df <- as.data.frame(mean_matrix_scaled)
mean_df$Family <- rownames(mean_df)
mean_long <- melt(mean_df,
                  id.vars = "Family",
                  variable.name = "Sample_Group",
                  value.name = "Scaled_Abundance")

mean_long$Family <- factor(
  mean_long$Family,
  levels = rev(rownames(mean_matrix_scaled))
)

# Original values for annotations
mean_original <- as.data.frame(mean_matrix_sorted)
mean_original$Family <- rownames(mean_original)
mean_original_long <- melt(mean_original,
                           id.vars = "Family",
                           variable.name = "Sample_Group",
                           value.name = "Original")

mean_long <- merge(
  mean_long, mean_original_long,
  by = c("Family", "Sample_Group")
)

scaled_heatmap_plot <- ggplot(
  mean_long, aes(x = Sample_Group, y = Family, fill = Scaled_Abundance)
) +
  geom_tile(color = "gray90", linewidth = 0.3) +
  geom_text(aes(label = round(Original, 3)), size = 2.5, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "red",
    limits = c(0, 1),
    name = "Relative\nabundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(title = paste("Top", indicator_heatmap_top_n, "Indicator Families"),
       subtitle = "Color intensity: white = minimum for this family, red = maximum for this family")  # nolint

ggsave(
  here(images_out, "Heatmap_Top15_Families_Scaled.pdf"),
  plot = scaled_heatmap_plot, width = 7, height = 10
)
