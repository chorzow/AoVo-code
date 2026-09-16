## Taxonomic composition (Phylum and Family level) analysis

source(here::here("scripts", "01-preload.R"))

data_count <- normalize_otu(otu_df)

## ---- helpers ----

# Collapse an aggregated taxon table to its top-N most abundant rows plus
# an "Other" row
collapse_top_n <- function(data_agg, n = top_n_taxa) {
  total <- rowSums(data_agg)
  total <- total / sum(total)
  data_agg <- cbind(data_agg, total = total)
  data_agg <- data_agg[order(data_agg$total, decreasing = TRUE), ]

  data_cum <- data_agg[seq_len(n), ]
  data_cum <- rbind(data_cum, colSums(data_agg[(n + 1):nrow(data_agg), ]))
  rownames(data_cum)[nrow(data_cum)] <- "Other"
  data_cum <- data_cum[, -ncol(data_cum)]  # remove total

  list(full = data_agg, cum = data_cum)
}

# Melt a per-sample x taxon table into long format for per-sample analyses.
to_long_pcm <- function(data_cum) {
  t_data_cum <- as.data.frame(t(data_cum))
  t_data_cum <- cbind(Sample = rownames(t_data_cum), t_data_cum)
  pcm <- melt(t_data_cum, id = "Sample")
  pcm$Sample <- factor(pcm$Sample, levels = unique(pcm$Sample))
  pcm
}

# write.table wrapper shared by every taxon-level table this script exports.
save_taxon_table <- function(df, id_col, filename) {
  out <- cbind(setNames(list(rownames(df)), id_col), df)
  write.table(out, here(data_out, filename),
              quote = FALSE, sep = "\t", row.names = FALSE)
}

# Grouped (by sample group) relative-abundance stacked bar chart.
plot_grouped_abundance <- function(data_cum, fill_label, palette = NULL) {
  sample_groups <- classify_group(colnames(data_cum))
  grouped <- t(apply(data_cum, 1, function(x) tapply(x, sample_groups, sum)))
  grouped_percent <- apply(grouped, 2, function(x) x / sum(x) * 100)

  plot_data <- melt(grouped_percent)
  colnames(plot_data) <- c(fill_label, "Sample_Group", "Abundance")
  plot_data$Sample_Group <- factor(plot_data$Sample_Group, levels = group_levels)

  fill_scale <- if (is.null(palette)) {
    scale_fill_brewer(palette = "Paired")
  } else {
    scale_fill_manual(values = palette)
  }

  p <- ggplot(
    plot_data, aes(x = Sample_Group, y = Abundance, fill = .data[[fill_label]])
  ) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
    fill_scale +
    labs(x = "Sample", y = "relative abundance (%)", fill = fill_label) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
      plot.caption = element_text(face = "italic", color = "red", hjust = 0.5),
      legend.position = "right"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    geom_hline(
      yintercept = c(25, 50, 75), linetype = "dashed", color = "gray", alpha = 0.3
    )

  list(plot = p, grouped = grouped, grouped_percent = grouped_percent)
}

# Coefficient of variation analysis + faceted boxplot of the most variable
# taxa within each sample group.
analyze_variability <- function(pcm, rank_label, plot_filename, table_prefix,
                                 print_n) {
  pcm$Group <- classify_group(pcm$Sample)
  pcm$Ecosystem <- classify_ecosystem(pcm$Sample)
  pcm$Plant <- classify_plant(pcm$Sample)

  cv_analysis <- pcm %>%
    group_by(Group, variable) %>%
    summarise(
      mean_abundance = mean(value),
      sd_abundance = sd(value),
      cv = sd(value) / mean(value) * 100,
      n = n(),
      .groups = "drop"
    )

  top_variable <- cv_analysis %>%
    group_by(Group) %>%
    arrange(desc(cv)) %>%
    slice_head(n = 10) %>%
    ungroup()

  cat("Top 10 most variable", rank_label, ", grouped:\n")
  print(top_variable, n = print_n)

  overall_variable <- cv_analysis %>%
    group_by(variable) %>%
    summarise(mean_cv = mean(cv)) %>%
    arrange(desc(mean_cv)) %>%
    slice_head(n = 8)

  variable_plot <- pcm %>%
    filter(variable %in% overall_variable$variable) %>%
    ggplot(aes(x = Group, y = value, fill = Ecosystem)) +
    geom_boxplot(outlier.shape = 16, alpha = 0.8) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.2), size = 1.5, alpha = 0.6
    ) +
    facet_wrap(~ variable, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = ecosystem_colors) +
    labs(x = "Sample", y = "abundance (%)", fill = "Biotope") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 8))

  ggsave(here(images_out, plot_filename), plot = variable_plot, width = 14, height = 10)

  write.table(cv_analysis,
    here(data_out, paste0(table_prefix, "_Variability_Analysis.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(top_variable,
    here(data_out, paste0("Top_Variable_", table_prefix, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  list(cv_analysis = cv_analysis, top_variable = top_variable, plot = variable_plot)
}

## ---- Phylum level ----

phylum_agg <- collapse_top_n(aggregate_by_taxon(data_count, tax_df$Phylum))
data_phylum <- phylum_agg$full
data_phylum_cum <- phylum_agg$cum

save_taxon_table(data_phylum, "Phylum", "16S_PhylumAOVO.tsv")
save_taxon_table(data_phylum_cum, "Phylum", "16S_Phylum_dominatAOVO.tsv")

pcm_phylum <- to_long_pcm(data_phylum_cum)

phylum_abundance <- plot_grouped_abundance(data_phylum_cum, "Phylum")
write.table(
  cbind(Phylum = rownames(phylum_abundance$grouped_percent), phylum_abundance$grouped_percent),
  here(data_out, "16S_Phylum_grouped_summary.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

phylum_plot <- phylum_abundance$plot
ggsave(
  here(images_out, "Dominant_Phyla_Grouped.pdf"),
  plot = phylum_plot, width = 10, height = 7
)

write.table(
  cbind(Phylum = rownames(phylum_abundance$grouped), phylum_abundance$grouped),
  here(data_out, "16S_Phylum_grouped_absolute.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

phylum_variability <- analyze_variability(
  pcm = pcm_phylum,
  rank_label = "phyla",
  plot_filename = "Phylum_Variability_Boxplot.pdf",
  table_prefix = "Phylum",
  print_n = 5
)

data_long <- data_phylum_cum %>%
  rownames_to_column(var = "Phylum") %>%
  pivot_longer(cols = -Phylum, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Group = factor(classify_group(Sample), levels = group_levels)) %>%
  filter(Abundance > 0)

cat("Groups found:", toString(unique(data_long$Group)), "\n")
cat("Number of observations in each group:\n")
print(table(data_long$Group))

# ---- Statistical analysis (Kruskal-Wallis + Dunn), Phylum level only ----

results_kw <- list()
results_dunn <- list()

phyla <- unique(data_long$Phylum)
for (phylum in phyla) {
  df <- data_long %>% filter(Phylum == phylum)

  kw <- kruskal.test(Abundance ~ Group, data = df)
  results_kw[[phylum]] <- data.frame(
    Phylum = phylum,
    Chi_squared = kw$statistic,
    p_value = kw$p.value,
    n = nrow(df)
  )

  if (kw$p.value < p_value_threshold) {
    dunn <- FSA::dunnTest(Abundance ~ Group, data = df, method = "bh")
    dunn_df <- dunn$res
    dunn_df$Phylum <- phylum
    results_dunn[[phylum]] <- dunn_df
  }
}

kw_results <- do.call(rbind, results_kw) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

dunn_results <- if (length(results_dunn) > 0) do.call(rbind, results_dunn) else NULL

cat("\n=== KRUSKAL-WALLIS RESULTS ===\n")
print(kw_results)
if (!is.null(dunn_results)) {
  cat("\n=== DUNN TEST RESULTS ===\n")
  print(dunn_results)
}

write.csv(kw_results, here(data_out, "kruskal_wallis_results.csv"), row.names = FALSE)
if (!is.null(dunn_results)) {
  write.csv(dunn_results, here(data_out, "dunn_test_results.csv"), row.names = FALSE)
}

## ---- Family level ----

data_family <- aggregate_by_taxon(
  data_count, paste(tax_df$Phylum, tax_df$Family, sep = "_")
)
family_agg <- collapse_top_n(data_family)
data_family <- family_agg$full
data_family_cum <- family_agg$cum

save_taxon_table(data_family, "Family", "16S_FamilyAOVO.tsv")
save_taxon_table(data_family_cum, "Family", "16S_Family_dominatAOVO.tsv")

pcm_family <- to_long_pcm(data_family_cum)

n_families <- nrow(data_family_cum)
family_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_families)

family_abundance <- plot_grouped_abundance(data_family_cum, "Family", palette = family_colors)
write.table(
  cbind(Family = rownames(family_abundance$grouped_percent), family_abundance$grouped_percent),
  here(data_out, "16S_Family_grouped_summary.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

family_plot <- family_abundance$plot +
  theme(legend.text = element_text(size = 8))
ggsave(
  here(images_out, "Dominant_Families_Grouped.pdf"),
  plot = family_plot, width = 12, height = 8
)

write.table(
  cbind(Family = rownames(family_abundance$grouped), family_abundance$grouped),
  here(data_out, "16S_Family_grouped_absolute.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

family_variability <- analyze_variability(
  pcm = pcm_family,
  rank_label = "families",
  plot_filename = "Family_Variability_Boxplot.pdf",
  table_prefix = "Family",
  print_n = 40
)

## ---- Combined Figure 3 ----

combined_plot <- (phylum_plot + labs(title = "A")) + (family_plot + labs(title = "B"))
ggsave(
  here(images_out, "Figure_3_Taxonomy_combined.pdf"),
  combined_plot, width = 14, height = 8
)

## ---- microshades composition plot (Phylum/Order stacked bars) ----

microshades_group_level    <- config$microshades$group_level
microshades_subgroup_level <- config$microshades$subgroup_level
microshades_selected_groups <- config$microshades$selected_groups

percent_fraction <- function(x) (x / sum(x) * 100)

get_relabund <- function(otu, tax, cur_samples, abund = 1) {
  ps <- phyloseq(otu, tax, cur_samples)
  cat("Transforming into percent abundance...\n")
  ps <- transform_sample_counts(ps, percent_fraction)
  cat("Filtering taxa...\n")
  ps_abund <- filter_taxa(ps, function(x) sum(x > abund) > 0, TRUE)
  cat("Transforming into relative percent...\n")
  transform_sample_counts(ps_abund, percent_fraction)
}

microshades_otu <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)
microshades_tax <- tax_table(as.matrix(tax_df))

abund <- get_relabund(microshades_otu, microshades_tax, sample_data(metadata), abund = 0)

theme_bars <- theme(
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.position = "right",
  axis.text = element_text(size = 12),
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 90, hjust = 1),
  axis.text.y = element_text(size = 12),
  strip.text = element_text(size = 12, face = "bold"),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12),
  strip.placement = "outside",
  panel.spacing.y = unit(1, "lines")
)

mdf <- prep_mdf(abund, subgroup_level = microshades_subgroup_level)
color_objs <- create_color_dfs(mdf,
  selected_groups = microshades_selected_groups,
  group_level = microshades_group_level,
  subgroup_level = microshades_subgroup_level,
  cvd = TRUE
)
mdf <- color_objs$mdf
cdf <- color_objs$cdf

reordered <- reorder_samples_by(mdf, cdf,
  sample_variable = "Sample",
  sample_ordering = sample_names
)

microshades_legend <- custom_legend(mdf, cdf,
  legend_key_size = 0.8, legend_text_size = 12,
  legend_orientation = "vertical",
  group_level = microshades_group_level, subgroup_level = microshades_subgroup_level
)

microshades_plot <- plot_microshades(reordered$mdf, reordered$cdf) +
  scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme_bars +
  theme(
    legend.position = "none",
    panel.spacing = unit(0, "lines"),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_blank()
  ) +
  facet_wrap(~Ecosystem, scales = "free_x", ncol = 1) +
  theme(axis.text.x = element_text(size = 8, angle = 65)) +
  theme(plot.margin = margin(6, 20, 6, 6)) +
  ylab("Relative abundance, %") +
  xlab("Sample type")

pdf(here(paste0(
  images_out, "/microshades_", str_to_lower(microshades_subgroup_level), ".pdf"
)))
print(microshades_plot)
dev.off()

pdf(here(paste0(
  images_out, "/microshades_", str_to_lower(microshades_subgroup_level), "_legend.pdf"
)), width = 10, height = 10)
print(microshades_legend)
dev.off()
