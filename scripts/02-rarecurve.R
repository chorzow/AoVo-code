## Rarefaction curve of observed OTUs per sample group.

source(here::here("scripts", "01-preload.R"))

otu_matrix <- as.matrix(otu_df)
groups <- rep(group_levels, times = c(10, 10, 10, 10))
colors <- group_colors

pdf(here(images_out, "Rarecurve.pdf"))

rarecurve(t(otu_matrix),
          step = rarecurve_step,
          col = colors[groups],
          lwd = 2,
          ylab = "Number of Observed OTUs",
          xlab = "Sequencing Depth (Number of Reads)",
          label = FALSE)

legend("bottomright",
       legend = names(colors),
       fill = colors,
       title = "Sample Group",
       bty = "n")

dev.off()
