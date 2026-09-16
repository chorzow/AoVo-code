###rarecurve########

library(here)

config <- config::get()
data_folder <- config$input_folder
metadata_path <- config$metadata_file
all_otu_path <- config$otu_file
all_tax_path <- config$tax_file
kits_order <- config$kits_order
output_folder <- config$output_folder  # TODO: rename to output_dir here and in config
images_out <- config$images_output_folder


otu_df <- read.csv(
  here(data_folder, all_otu_path),
  header = TRUE, row.names = 1, check.names = FALSE, sep = "\t"
)
otu_matrix <- as.matrix(otu_df)
sample_names <- colnames(otu_matrix)
groups <- rep(c("S1_bulk", "S1_rAo", "S2_bulk", "S2_rVo"),
              times = c(10, 10, 10, 10))

colors <- brewer.pal(n = 4, name = "Set1")
names(colors) <- unique(groups)

rarecurve(t(otu_matrix),
          step = 100,
          col = colors[groups],
          lwd = 2,
          ylab = "Number of Observed OTUs",
          xlab = "Sequencing Depth (Number of Reads)",
          label = FALSE)

legend("bottomright",
       legend = names(colors),
       fill = colors,
       title = "Sample Group", # or "Habitat"
       bty = "n")

###Alpha-diversity metrics###
t_otu_df <- as.data.frame(t(otu_df))
diversity_16s <- t(estimateR(t_otu_df))
shannon <- diversity(data, index = "shannon")
simpson <- diversity(data, index = "simpson")
diversity_16s <- cbind(
  sample_id = rownames(diversity_16s),
  reads_number = rowSums(t_otu_df),
  diversity_16s,
  shannon = shannon,
  simpson = simpson
)
diversity_16s <- as.data.frame(diversity_16s)
write.table(diversity_16s, here(output_folder, "data/16S_diversity.tsv"),
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
    Sample_Type = case_when(
      grepl("S1_bulk", sample_id) ~ "S1_bulk",
      grepl("S1_rAo", sample_id)   ~ "S1_rAo",
      grepl("S2_bulk", sample_id) ~ "S2_bulk",
      grepl("S2_rVo", sample_id)   ~ "S2_rVo"
    ),
    Ecosystem = case_when(
      grepl("bulk", sample_id) ~ "Soil",
      grepl("rAo|rVo", sample_id) ~ "Rhizosphere"
    ),
    Plant = case_when(
      grepl("S1_rAo|S1_bulk", sample_id) ~ "A.olchonensis",
      grepl("S2_rVo|S2_bulk", sample_id) ~ "V.olchonensis"
    )
  )

diversity_16s$Sample_Type <- factor(diversity_16s$Sample_Type,
                                    levels = c("S1_bulk", "S1_rAo",
                                               "S2_bulk", "S2_rVo"))
diversity_16s$Ecosystem <- factor(diversity_16s$Ecosystem,
                                  levels = c("Soil", "Rhizosphere"))
diversity_16s$Plant <- factor(diversity_16s$Plant,
                              levels = c("A.olchonensis", "V.olchonensis"))
(chao1_boxplot <- ggplot(
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
  scale_fill_manual(
    values = c("Soil" = "#E69F00", "Rhizosphere" = "#56B4E9")
  ) +
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
)

(shannon_boxplot <- ggplot(
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
  scale_fill_manual(
    values = c("Soil" = "#E69F00", "Rhizosphere" = "#56B4E9")
  ) +
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

### beta-diversity ###
metadata <- data.frame(  # TODO: do that on the main metadata df
  Sample = sample_names,
  Group = ifelse(grepl("S1_bulk", sample_names), "S1_bulk",
                 ifelse(grepl("S1_rAo", sample_names), "S1_rAo",
                        ifelse(grepl("S2_bulk", sample_names), "S2_bulk",
                               "S2_rVo"))),
  Ecosystem = ifelse(grepl("bulk", sample_names), "Soil", "Rhizosphere"),
  Plant = ifelse(grepl("rAo", sample_names), "Astragalus", "Vicia")
)
rownames(metadata) <- sample_names


data_hell <- decostand(t_otu_df, method = "hellinger")
dist_bray <- vegdist(data_hell, method = "bray")

nmds <- metaMDS(dist_bray, k = 2, trymax = 100)

nmds_stress <- nmds$stress
cat("Stress value NMDS:", nmds_stress, "\n")

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)
nmds_points <- merge(nmds_scores, metadata, by = "Sample")

names(nmds_points)[names(nmds_points) == "NMDS1"] <- "NMDS1"
names(nmds_points)[names(nmds_points) == "NMDS2"] <- "NMDS2"

centroids_nmds <- nmds_points %>%
  group_by(Group) %>%
  summarise(
    Centroid_NMDS1 = mean(NMDS1),
    Centroid_NMDS2 = mean(NMDS2)
  )

# NMDS plot
(nmds_plot <- ggplot(
  nmds_points, aes(x = NMDS1, y = NMDS2, color = Group, shape = Ecosystem)
) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(
    data = centroids_nmds,
    aes(x = Centroid_NMDS1, y = Centroid_NMDS2, color = Group),
    size = 7, shape = 8, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Group), level = 0.95, linetype = 2, size = 0.5, alpha = 0.7
  ) +
  scale_color_manual(values = c("S1_bulk" = "#E41A1C",
                                "S1_rAo" = "#377EB8",
                                "S2_bulk" = "#4DAF4A",
                                "S2_rVo" = "#984EA3")
  ) +
  labs(
    title = "NMDS based on Bray-Curtis distance",
    subtitle = paste("Stress =", round(nmds_stress, 3)),
    x = "NMDS Axis 1",
    y = "NMDS Axis 2",
    color = "Sample group",
    shape = "Sample type"
  ) +
  theme_bw() +
  theme(legend.position = "right")
)

ggsave(
  here(images_out, "BetaDiversity_NMDS_Plot.png"),
  plot = nmds_plot, width = 10, height = 8, dpi = 300
)

# PCoA (Principal Coordinate Analysis)
pcoa <- cmdscale(dist_bray, k = 2, eig = TRUE)
pcoa_scores <- as.data.frame(pcoa$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$Sample <- rownames(pcoa_scores)
eigenvals <- pcoa$eig
var_explained <- eigenvals / sum(eigenvals) * 100

cat("\n--- PCoA Results ---\n")
cat("PCoA1 explains:", round(var_explained[1], 2), "% of variance\n")
cat("PCoA2 explains:", round(var_explained[2], 2), "% of variance\n")


pcoa_points <- merge(pcoa_scores, metadata, by = "Sample")


centroids_pcoa <- pcoa_points %>%
  group_by(Group) %>%
  summarise(
    Centroid_PCoA1 = mean(PCoA1),
    Centroid_PCoA2 = mean(PCoA2)
  )


(pcoa_plot <- ggplot(
  pcoa_points, aes(x = PCoA1, y = PCoA2, color = Group, shape = Ecosystem)
) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(
    data = centroids_pcoa,
    aes(x = Centroid_PCoA1, y = Centroid_PCoA2, color = Group),
    size = 7, shape = 8, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Group), level = 0.95, linetype = 2, size = 0.5, alpha = 0.7
  ) +
  scale_color_manual(values = c("S1_bulk" = "#E41A1C",
                                "S1_rAo" = "#377EB8",
                                "S2_bulk" = "#4DAF4A",
                                "S2_rVo" = "#984EA3")
  ) +
  labs(
    title = "PCoA based on Bray-Curtis distance",
    subtitle = paste0("PCoA1: ", round(var_explained[1], 1), "%, ",
                      "PCoA2: ", round(var_explained[2], 1), "%"),
    x = paste0("PCoA Axis 1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCoA Axis 2 (", round(var_explained[2], 1), "%)"),
    color = "Sample group",
    shape = "Sample type"
  ) +
  theme_bw() +
  theme(legend.position = "right")
)

ggsave(
  here(output_folder, "images", "BetaDiversity_PCoA_Plot.png"),
  plot = pcoa_plot, width = 10, height = 8, dpi = 300
)


# ANOSIM
anosim_result <- anosim(dist_bray, metadata$Group)
cat("\n--- ANOSIM Results ---\n")
print(anosim_result)

# PERMANOVA (simplified model using only one Group factor)
permanova_result <- adonis2(dist_bray ~ Group, data = metadata)
cat("\n--- PERMANOVA Results (Group) ---\n")
print(permanova_result)


write.csv(nmds_points,
  here(output_folder, "data", "NMDS_Coordinates.csv"), row.names = FALSE
)
write.csv(centroids_nmds,
  here(output_folder, "data", "NMDS_Centroids.csv"), row.names = FALSE
)
write.csv(pcoa_points,
  here(output_folder, "data", "PCoA_Coordinates.csv"), row.names = FALSE
)
write.csv(centroids_pcoa,
  here(output_folder, "data", "PCoA_Centroids.csv"), row.names = FALSE
)
write.csv(as.matrix(dist_bray),
  here(output_folder, "data", "BrayCurtis_Distance_Matrix.csv")
)
cat("\n--- All analyses completed ---\n")


#Taxonomy_Phylum #####
tax_df <- read.csv(
  here(data_folder, all_tax_path),
  header = TRUE, row.names = 1, check.names = FALSE, sep = "\t"
)

data_count <- otu_df
colsums <- colSums(data_count)
for (i in seq_len(ncol(data_count))) data_count[, i] <- data_count[, i] / colsums[i]  # TODO: rewrite
colSums(data_count)

data_phylum <- cbind(data_count, Phylum = paste(tax_df$Phylum))
data_phylum <- aggregate(. ~ Phylum, data_phylum, sum)
rownames(data_phylum) <- data_phylum[, 1]
data_phylum <- data_phylum[, -1]
total <- rowSums(data_phylum)
total <- total / sum(total)

data_phylum <- cbind(data_phylum, total = total)
data_phylum <- data_phylum[order(data_phylum$total, decreasing = TRUE), ]
data_phylum_cum <- data_phylum[1:10, ]  # extract top 10 dominant phyla
data_phylum_cum <- rbind(
  data_phylum_cum, colSums(data_phylum[11:nrow(data_phylum), ])
)

rownames(data_phylum_cum)[nrow(data_phylum_cum)] <- "Other"  # for unassigned
data_phylum_cum <- data_phylum_cum[, -ncol(data_phylum_cum)]  # remove total

t_data_phylum_cum <- as.data.frame(t(data_phylum_cum))
t_data_phylum_cum <- cbind(
  Sample = rownames(t_data_phylum_cum), t_data_phylum_cum
)

pcm <- melt(t_data_phylum_cum, id = "Sample")
pcm$Sample <- factor(pcm$Sample, levels = unique(pcm$Sample))

br <- ggplot(pcm, aes(x = Sample, y = value, fill = variable))
br +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Sample", y = "relative abundance")

write.table(
  cbind(Phylum = rownames(data_phylum), data_phylum),
  here(output_folder, "data", "16S_PhylumAOVO.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(
  cbind(Phylum = rownames(data_phylum_cum), data_phylum_cum),
  here(output_folder, "data", "16S_Phylum_dominatAOVO.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

data_phylum_cum <- data_phylum_cum[, -ncol(data_phylum_cum)]

sample_groups <- ifelse(
  grepl("S1_bulk", colnames(data_phylum_cum)), "S1_bulk",
  ifelse(grepl("rAo", colnames(data_phylum_cum)), "S1_rAo",
    ifelse(grepl("S2_bulk", colnames(data_phylum_cum)), "S2_bulk",
      "S2_rVo"
    )
  )
)

grouped_phyla <- t(apply(data_phylum_cum, 1, function(x) {
  tapply(x, sample_groups, sum)
}))

grouped_phyla_percent <- apply(grouped_phyla, 2, function(x) x / sum(x) * 100)

write.table(cbind(Phylum = rownames(grouped_phyla_percent),
                  grouped_phyla_percent),
            here(output_folder, "data", "16S_Phylum_grouped_summary.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)


plot_data <- melt(grouped_phyla_percent)
colnames(plot_data) <- c("Phylum", "Sample_Group", "Abundance")

plot_data$Sample_Group <- factor(plot_data$Sample_Group,
                                 levels = c("S1_bulk", "S1_rAo",
                                            "S2_bulk", "S2_rVo"))


(phylum_plot <- ggplot(
  plot_data, aes(x = Sample_Group, y = Abundance, fill = Phylum)
) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
  scale_fill_brewer(palette = "Paired") +
  labs(
    x = "Sample",
    y = "relative abundance (%)",
    fill = "Phylum"
  ) +
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
    yintercept = c(25, 50, 75), linetype = "dashed",
    color = "gray", alpha = 0.3
  )
)

ggsave(
  here(output_folder), "images", "Dominant_Phyla_Grouped.png",
  plot = phylum_plot, width = 10, height = 7, dpi = 300
)

write.table(cbind(Phylum = rownames(grouped_phyla),
                  grouped_phyla),
            here(output_folder, "data", "16S_Phylum_grouped_absolute.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)


pcm$Group <- ifelse(grepl("S1_bulk", pcm$Sample), "S1_bulk",
                    ifelse(grepl("rAo", pcm$Sample), "  S1_rAo",
                           ifelse(grepl("S2_bulk", pcm$Sample), "S2_bulk",
                                  "S2_rVo")))

pcm$Ecosystem <- ifelse(grepl("bulk", pcm$Sample), "Soil", "Rhizosphere")
pcm$Plant <- ifelse(grepl("Ao|Vo", pcm$Sample), "Astragalus", "Vicia")


cv_analysis <- pcm %>%
  group_by(Group, variable) %>%
  summarise(
    mean_abundance = mean(value),
    sd_abundance = sd(value),
    cv = sd(value) / mean(value) * 100,  # Coefficient of Variation in %
    n = n(),
    .groups = "drop"
  )


top_variable_phylum <- cv_analysis %>%
  group_by(Group) %>%
  arrange(desc(cv)) %>%
  slice_head(n = 10) %>%
  ungroup()

print("Top 10 most variable phyla, grouped:")
print(top_variable_phylum, n = 5)


overall_variable <- cv_analysis %>%
  group_by(variable) %>%
  summarise(mean_cv = mean(cv)) %>%
  arrange(desc(mean_cv)) %>%
  slice_head(n = 8)

variable_plot <- pcm %>%
  filter(variable %in% overall_variable$variable) %>%
  ggplot(aes(x = Group, y = value, fill = Ecosystem)) +
  geom_boxplot(outlier.shape = 16, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 1.5, alpha = 0.6) +
  facet_wrap(~ variable, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c("Soil" = "#E69F00", "Rhizosphere" = "#56B4E9")) +
  labs(x = "Sample",
       y = "abundance (%)",
       fill = "Biotope") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8))

ggsave(
  here(output_folder, "images", "Phylum_Variability_Boxplot.png"),
  plot = variable_plot, width = 14, height = 10, dpi = 300
)


write.table(cv_analysis,
  here(output_folder, "data", "Phylum_Variability_Analysis.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(top_variable_phylum,
  here(output_folder, "data", "Top_Variable_Phylum.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


sample_names <- colnames(data_phylum_cum)
sample_groups <- gsub("[0-9]+$", "", sample_names)
group_levels <- c("S1_bulk", "S1_rAo", "S2_bulk", "S2_rVo")

data_long <- data_phylum_cum %>%
  rownames_to_column(var = "Phylum") %>%
  pivot_longer(cols = -Phylum, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Group = gsub("[0-9]+$", "", Sample)) %>%
  mutate(Group = factor(Group, levels = group_levels)) %>%
  filter(Abundance > 0)

cat("Groups found:", toString(unique(data_long$Group)), "\n")
cat("Number of observations in each group:\n")
print(table(data_long$Group))


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


  if (kw$p.value < 0.05) {
    dunn <- FSA::dunnTest(Abundance ~ Group, data = df, method = "bh")
    dunn_df <- dunn$res
    dunn_df$Phylum <- phylum
    results_dunn[[phylum]] <- dunn_df
  }
}

kw_results <- do.call(rbind, results_kw) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

if (length(results_dunn) > 0) {
  dunn_results <- do.call(rbind, results_dunn)
} else {
  dunn_results <- NULL
}


cat("\n=== KRUSKAL-WALLIS RESULTS ===\n")
print(kw_results)

if (!is.null(dunn_results)) {
  cat("\n=== DUNN TEST RESULTS ===\n")
  print(dunn_results)
}


write.csv(kw_results, "kruskal_wallis_results.csv", row.names = FALSE)
if (!is.null(dunn_results)) {
  write.csv(dunn_results, "dunn_test_results.csv", row.names = FALSE)
}

#Taxonomy_Family####

data_count <- otu_df
colsums <- colSums(data_count)
for (i in seq_len(ncol(data_count))) data_count[, i] <- data_count[, i] / colsums[i]  # TODO: rewrite
colSums(data_count)

data_family <- cbind(data_count,
  Family = paste(tax_df$Phylum, tax_df$Family, sep = "_")
)

data_family <- aggregate(. ~ Family, data_family, sum)
rownames(data_family) <- data_family[, 1]
data_family <- data_family[, -1]
total <- rowSums(data_family)
total <- total / sum(total)

data_family <- cbind(data_family, total = total)
data_family <- data_family[order(data_family$total, decreasing = TRUE), ]
data_family_cum <- data_family[1:10, ]  # extracting top 10 dominating families
data_family_cum <- rbind(
  data_family_cum, colSums(data_family[11:nrow(data_family), ])
)
rownames(data_family_cum)[nrow(data_family_cum)] <- "Other"  # for unassigned
data_family_cum <- data_family_cum[, -ncol(data_family_cum)]  # remove total

t_data_family_cum <- as.data.frame(t(data_family_cum))
t_data_family_cum <- cbind(
  Sample = rownames(t_data_family_cum), t_data_family_cum
)
pcm <- melt(t_data_family_cum, id = "Sample")
pcm$Sample <- factor(pcm$Sample, levels = unique(pcm$Sample))

br <- ggplot(pcm, aes(x = Sample, y = value, fill = variable))
br +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Dominating families", x = "Samples", y = "Relative abundance")

write.table(cbind(Family = rownames(data_family), data_family),
  here(output_folder, "data", "16S_FamilyAOVO.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)
write.table(cbind(Family = rownames(data_family_cum), data_family_cum),
  here(output_folder, "data", "16S_Family_dominatAOVO.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)

data_family_cum <- data_family_cum[, -ncol(data_family_cum)]

sample_groups <- ifelse(
  grepl("S1_bulk", colnames(data_family_cum)), "S1_bulk",
  ifelse(grepl("rAo", colnames(data_family_cum)), "S1_rAo",
    ifelse(grepl("S2_bulk", colnames(data_family_cum)), "S2_bulk",
      "S2_rVo"
    )
  )
)


grouped_families <- t(apply(data_family_cum, 1, function(x) {
  tapply(x, sample_groups, sum)
}))


grouped_families_percent <- apply(
  grouped_families, 2, function(x) x / sum(x) * 100
)


write.table(cbind(Family = rownames(grouped_families_percent),
                  grouped_families_percent),
            here(output_folder, "data", "16S_Family_grouped_summary.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)


plot_data <- melt(grouped_families_percent)
colnames(plot_data) <- c("Family", "Sample_Group", "Abundance")

plot_data$Sample_Group <- factor(plot_data$Sample_Group,
                                 levels = c("S1_bulk", "S1_rAo",
                                            "S2_bulk", "S2_rVo"))


n_families <- length(unique(plot_data$Family))
family_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_families)

(family_plot <- ggplot(
  plot_data, aes(x = Sample_Group, y = Abundance, fill = Family)
) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.9) +
  scale_fill_manual(values = family_colors) +
  labs(
    x = "Sample",
    y = "relative abundance (%)",
    fill = "Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
    plot.caption = element_text(face = "italic", color = "red", hjust = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 8)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_hline(
    yintercept = c(25, 50, 75), linetype = "dashed",
    color = "gray", alpha = 0.3
  )
)


ggsave(
  here(output_folder, "images", "Dominant_Families_Grouped.png"),
  plot = family_plot, width = 12, height = 8, dpi = 300
)


write.table(
  cbind(Family = rownames(grouped_families), grouped_families),
  here(output_folder, "data", "16S_Family_grouped_absolute.tsv"),
  quote = FALSE, sep = "\t", row.names = FALSE
)


pcm$Group <- ifelse(grepl("S1_bulk", pcm$Sample), "S1_bulk",
                    ifelse(grepl("rAo", pcm$Sample), "  S1_rAo",
                           ifelse(grepl("S2_bulk", pcm$Sample), "S2_bulk",
                                  "S2_rVo")))

pcm$Ecosystem <- ifelse(grepl("soil", pcm$Sample), "Soil", "Rhizosphere")
pcm$Plant <- ifelse(grepl("Ao|soil1", pcm$Sample), "Astragalus", "Vicia")



cv_analysis <- pcm %>%
  group_by(Group, variable) %>%
  summarise(
    mean_abundance = mean(value),
    sd_abundance = sd(value),
    cv = sd(value) / mean(value) * 100, # Coefficient of Variation in %
    n = n(),
    .groups = "drop"
  )


top_variable_families <- cv_analysis %>%
  group_by(Group) %>%
  arrange(desc(cv)) %>%
  slice_head(n = 10) %>%
  ungroup()

print("Top 10 most variable families, grouped:")
print(top_variable_families, n = 40)



write.table(cv_analysis,
  here(output_folder, "data", "Family_Variability_Analysis.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(top_variable_families,
  here(output_folder, "data", "Top_Variable_Families.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


combined_plot <- (phylum_plot + labs(title = "A")) +
  (family_plot + labs(title = "B"))
combined_plot
ggsave(
  here(output_folder, "images", "Figure_3_Taxonomy_combined.png"),
  combined_plot, width = 14, height = 8, dpi = 300
)



## VENN DIAGRAMS FOR ALL FAMILIES (FULL MICROBIOME COMPARISON)#####

data_family_raw <- read.delim(
  here(output_folder, "data", "16S_FamilyAOVO.tsv"),
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

if ("total" %in% colnames(data_family_raw)) {
  data_family <- data_family_raw[, !colnames(data_family_raw) %in% "total"]
} else {
  data_family <- data_family_raw
}

data_family <- as.data.frame(lapply(data_family, as.numeric))
rownames(data_family) <- rownames(data_family_raw)

s1_bulk_samples <- colnames(data_family)[grepl("S1_bulk",
                                               colnames(data_family))]
s1_rao_samples <- colnames(data_family)[grepl("S1_rAo",
                                              colnames(data_family))]
s2_bulk_samples <- colnames(data_family)[grepl("S2_bulk",
                                               colnames(data_family))]
s2_rvo_samples <- colnames(data_family)[grepl("S2_rVo",
                                              colnames(data_family))]

get_families_for_group <- function(sample_names) {
  if (length(sample_names) == 0) return(character(0))
  presence_matrix <- data_family[, sample_names, drop = FALSE] > 0
  families_present <- rownames(data_family)[apply(presence_matrix, 1, any)]
  families_present
}

families_s1_bulk <- get_families_for_group(s1_bulk_samples)
families_s1_rao <- get_families_for_group(s1_rao_samples)
families_s2_bulk <- get_families_for_group(s2_bulk_samples)
families_s2_rvo <- get_families_for_group(s2_rvo_samples)

jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  intersection / union
}

jaccard_sites <- jaccard_index(families_s1_bulk, families_s2_bulk)
jaccard_s1 <- jaccard_index(families_s1_bulk, families_s1_rao)
jaccard_s2 <- jaccard_index(families_s2_bulk, families_s2_rvo)

create_fixed_venn <- function(set1, set2, labels, colors, title) {

  area1 <- length(set1)
  area2 <- length(set2)
  cross_area <- length(intersect(set1, set2))

  venn_plot <- draw.pairwise.venn(
    area1 = area1,
    area2 = area2,
    cross.area = cross_area,
    category = labels,
    fill = colors,
    alpha = 0.7,
    lty = "blank",
    cex = 1.2,
    cat.cex = 1.1,
    cat.pos = c(0, 0),
    cat.dist = 0.05,
    scaled = TRUE,
    fontfamily = "sans",
    cat.fontfamily = "sans",
    ind = FALSE
  )

  venn_grob <- gTree(children = venn_plot)

  result <- arrangeGrob(
    textGrob(title, gp = gpar(fontface = "bold", fontsize = 14), just = "top"),
    venn_grob,
    nrow = 2,
    heights = c(0.1, 0.9)
  )
  result
}

venn_sites_fixed <- create_fixed_venn(
  set1 = families_s1_bulk,
  set2 = families_s2_bulk,
  labels = c("S1_bulk", "S2_bulk"),
  colors = c("#FF9999", "#99FF99"),
  title = "A"
)

venn_s1_fixed <- create_fixed_venn(
  set1 = families_s1_bulk,
  set2 = families_s1_rao,
  labels = c("S1_bulk", "S1_rAo"),
  colors = c("#E41A1C", "#377EB8"),
  title = "B"
)

venn_s2_fixed <- create_fixed_venn(
  set1 = families_s2_bulk,
  set2 = families_s2_rvo,
  labels = c("S2_bulk", "S2_rVo"),
  colors = c("#4DAF4A", "#984EA3"),
  title = "C"
)

combined_fixed <- arrangeGrob(
  venn_sites_fixed,
  arrangeGrob(venn_s1_fixed, venn_s2_fixed, ncol = 2),
  nrow = 2,
  heights = c(1, 1)
)

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

output_path <- here(output_folder, "images", "Figure_Venn_diagrams_fixed.png")

ggsave(output_path,
       combined_fixed,
       width = 14,
       height = 10,
       dpi = 300)

grid.newpage()
grid.draw(combined_fixed)


###IndVal#####################

data_count <- otu_df
colsums <- colSums(data_count)
for (i in seq_len(ncol(data_count))) data_count[, i] <- data_count[, i] / colsums[i]
colSums(data_count)
data_family <- cbind(
  data_count, Family = paste(tax_df$Order, tax_df$Family, sep="_")
)
data_family <- aggregate(. ~ Family, data_family, sum)
rownames(data_family) <- data_family[, 1]
data_family <- data_family[, -1]

otu_table_for_indicator <- as.data.frame(t(data_family))

sample_names <- rownames(otu_table_for_indicator)
groups <- ifelse(
  grepl("S1_bulk", sample_names), "S1_bulk",
  ifelse(grepl("rAo", sample_names), "S1_rAo",
    ifelse(grepl("S2_bulk", sample_names), "S2_bulk",
      ifelse(grepl("rVo", sample_names),
        "S2_rVo", NA
      )
    )
  )
)

indval_family <- multipatt(otu_table_for_indicator,
                           cluster = groups,
                           func = "r.g",
                           duleg = FALSE,
                           control = how(nperm = 999))


indval_result <- indval_family$sign
indval_result$p.value <- indval_family$sign$p.value
indval_result$Family <- rownames(indval_result)


signif_indval_result <- indval_result %>%
  filter(p.value <= 0.05) %>%
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
  filter(stat > 0.5) %>%
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
  here(output_folder, "images", "Barplot_Strong_Indicators.png"),
  strong_barplot, width = 10, height = 8, dpi = 300
)


top_families <- strong_indicators %>% head(15) %>% pull(Family)

heatmap_data <- otu_table_for_indicator[, top_families, drop = FALSE]
heatmap_data_scaled <- as.data.frame(scale(heatmap_data))
heatmap_data_scaled$Group <- groups

heatmap_data_avg <- heatmap_data_scaled %>%
  group_by(Group) %>%
  summarise(across(all_of(top_families), mean)) %>%
  tibble::column_to_rownames("Group")

heatmap_data_long <- heatmap_data_avg %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Group") %>%
  pivot_longer(cols = -Group, names_to = "Family", values_to = "Z_score")

heatmap_plot <- ggplot(
  heatmap_data_long, aes(x = Group, y = Family, fill = Z_score)
) +
  geom_tile(color = "white", size = 0.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Z-score\n(relative abundance)") +
  labs(x = "Samples Group", y = "Family") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8))

ggsave(
  here(output_folder, "images", "Heatmap_Strong_Indicators.png"),
  heatmap_plot, width = 10, height = 8, dpi = 300
)


write.csv(
  signif_indval_result,
  here(output_folder, "data", "All_Indicator_Results.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)
write.csv(
  strong_indicators,
  here(output_folder, "data", "Strong_Indicators_Results.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

cat("Analysis done.\n")
cat(paste(
  "Total count of indicator species (p < 0.05):",
  nrow(signif_indval_result),
  "\n"
))
cat(paste(
  "Count of strong indicators (stat > 0.5):",
  nrow(strong_indicators),
  "\n"
))
cat("Distribution by group:\n")
print(table(strong_indicators$Indicator_Group))

# BOXPLOTS FOR TOP-3 INDICATOR FAMILIES FOR EACH GROUP
strong_indicators <- read.csv(
  here(output_folder, "data", "Strong_Indicators_Results.csv"),
  stringsAsFactors = FALSE
)

target_groups <- c("S1_bulk", "S1_rAo", "S2_bulk", "S2_rVo")


top_indicators_list <- list()

for (group in target_groups) {
  group_indicators <- strong_indicators %>%
    filter(Indicator_Group == group) %>%
    arrange(desc(stat)) %>%
    head(3)  # TODO: n indicators to config?

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

  p <- ggplot(
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

  p
}



plots_by_group <- list()
for (group in target_groups) {
  families <- top_indicators_list[[group]]
  if (length(families) > 0) {
    group_plots <- list()
    for (family in families) {
      group_plots[[family]] <- create_family_boxplot(family)
    }
    plots_by_group[[group]] <- group_plots
  }
}


create_row_with_title <- function(plot_list, title_text) {
  row_plot <- ggarrange(plotlist = plot_list, ncol = length(plot_list))
  row_plot <- annotate_figure(row_plot,
                              top = text_grob(title_text,
                                              face = "bold",
                                              size = 12,
                                              hjust = 0,
                                              x = 0))
  row_plot
}


row_s1_bulk <- create_row_with_title(
  plots_by_group[["S1_bulk"]],
  "A. S1_bulk (Psammozem soil) indicator families"
)

row_s1_rao <- create_row_with_title(
  plots_by_group[["S1_rAo"]],
  "B. S1_rAo (A. olchonensis rhizosphere) indicator families"
)

row_s2_bulk <- create_row_with_title(
  plots_by_group[["S2_bulk"]],
  "C. S2_bulk (Litozem soil) indicator families"
)

row_s2_rvo <- create_row_with_title(
  plots_by_group[["S2_rVo"]],
  "D. S2_rVo (V. olchonensis rhizosphere) indicator families"
)


combined_plot <- ggarrange(row_s1_bulk, row_s1_rao, row_s2_bulk, row_s2_rvo,
                           nrow = 4,
                           ncol = 1)

ggsave(
  here(output_folder, "images", "Indicator_Families_Boxplots.png"),
  combined_plot, width = 12, height = 14, dpi = 300
)

print(combined_plot)

# Save combined plot
ggsave(
  here(output_folder, "images", "Combined_Boxplots_Top_Indicators.png"),
  plot = combined_plot,
  width = n_cols * 5,  # TODO: n_cols missing
  height = n_rows * 4.5,  # TODO: n_rows missing
  dpi = 300
)

# Save each plot separately (optional)
for (family in all_top_families) {
  safe_name <- gsub("[^[:alnum:]]", "_", family)
  ggsave(
    here(output_folder, "images", paste0("Boxplot_", safe_name, ".png")),
    plot = all_plots[[family]],  # TODO: all_plots missing
    width = 6, height = 5, dpi = 300
  )
}


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
  here(output_folder, "data", "Top_Indicators_Boxplot_Statistics.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)


# cat("\nPLOTS SUCCESSFULLY CREATED\n")
# cat(paste("Boxplots created:", length(all_plots), "\n"))
# cat(paste("Combined plot saved: 'Combined_Boxplots_Top_Indicators.png'\n"))
# cat(paste("Statistics saved: 'Top_Indicators_Boxplot_Statistics.csv'\n\n"))

cat("TOP INDICATORS FOR PLOTS:\n")
for (group in target_groups) {
  if (!is.null(top_indicators_list[[group]])) {
    cat(paste(
      group, ":", paste(top_indicators_list[[group]], collapse = ", "), "\n"
    ))
  }
}

# Display combined plot
print(combined_plot)

# SCALING 0-1 FOR EACH FAMILY

stats <- read.csv(
  here(output_folder, "data", "Top_Indicators_Boxplot_Statistics.csv"),
  stringsAsFactors = FALSE
)


target_groups <- c("S1_bulk", "S1_rAo", "S2_bulk", "S2_rVo")
stats_filtered <- stats %>% filter(Sample_Group %in% target_groups)


top_families <- stats_filtered %>%
  group_by(Family) %>%
  summarise(max_mean = max(Mean)) %>%
  arrange(desc(max_mean)) %>%
  head(15)

stats_top <- stats_filtered %>% filter(Family %in% top_families$Family)


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

a <- ggplot(
  mean_long, aes(x = Sample_Group, y = Family, fill = Scaled_Abundance)
) +
  geom_tile(color = "gray90", size = 0.3) +
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
  labs(title = "Top 15 Indicator Families",
       subtitle = "Color intensity: white = minimum for this family, red = maximum for this family")  # nolint

ggsave(
  here(output_folder, "images", "Heatmap_Top15_Families_Scaled.png"),
  width = 7, height = 10, dpi = 300
)
