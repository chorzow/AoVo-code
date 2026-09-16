## Beta-diversity: NMDS, PCoA, ANOSIM and PERMANOVA on Bray-Curtis distance.

source(here::here("scripts", "01-preload.R"))

sample_meta <- data.frame(
  Sample = sample_names,
  Group = classify_group(sample_names),
  Plant = classify_plant(sample_names)
)


metadata <- merge(sample_meta, metadata, by = "Sample")
metadata <- metadata[match(sample_names, metadata$Sample), ]
rownames(metadata) <- metadata$Sample

data_hell <- decostand(t_otu_df, method = "hellinger")
dist_bray <- vegdist(data_hell, method = "bray")

nmds <- metaMDS(dist_bray, k = nmds_k, trymax = nmds_trymax)

nmds_stress <- nmds$stress
cat("Stress value NMDS:", nmds_stress, "\n")

nmds_scores <- as.data.frame(scores(nmds))
nmds_scores$Sample <- rownames(nmds_scores)
nmds_points <- merge(nmds_scores, metadata, by = "Sample")

centroids_nmds <- nmds_points %>%
  group_by(Group) %>%
  summarise(
    Centroid_NMDS1 = mean(NMDS1),
    Centroid_NMDS2 = mean(NMDS2)
  )

nmds_plot <- ggplot(
  nmds_points, aes(x = NMDS1, y = NMDS2, color = Group, shape = Ecosystem)
) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(
    data = centroids_nmds,
    aes(x = Centroid_NMDS1, y = Centroid_NMDS2, color = Group),
    size = 7, shape = 8, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Group), level = 0.95, linetype = 2, linewidth = 0.5, alpha = 0.7
  ) +
  scale_color_manual(values = group_colors) +
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

ggsave(
  here(images_out, "BetaDiversity_NMDS_Plot.pdf"),
  plot = nmds_plot, width = 10, height = 8
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

pcoa_plot <- ggplot(
  pcoa_points, aes(x = PCoA1, y = PCoA2, color = Group, shape = Ecosystem)
) +
  geom_point(size = 4, alpha = 0.8) +
  geom_point(
    data = centroids_pcoa,
    aes(x = Centroid_PCoA1, y = Centroid_PCoA2, color = Group),
    size = 7, shape = 8, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Group), level = 0.95, linetype = 2, linewidth = 0.5, alpha = 0.7
  ) +
  scale_color_manual(values = group_colors) +
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

ggsave(
  here(images_out, "BetaDiversity_PCoA_Plot.pdf"),
  plot = pcoa_plot, width = 10, height = 8
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
  here(data_out, "NMDS_Coordinates.csv"), row.names = FALSE
)
write.csv(centroids_nmds,
  here(data_out, "NMDS_Centroids.csv"), row.names = FALSE
)
write.csv(pcoa_points,
  here(data_out, "PCoA_Coordinates.csv"), row.names = FALSE
)
write.csv(centroids_pcoa,
  here(data_out, "PCoA_Centroids.csv"), row.names = FALSE
)
write.csv(as.matrix(dist_bray),
  here(data_out, "BrayCurtis_Distance_Matrix.csv")
)
cat("\n--- All analyses completed ---\n")
