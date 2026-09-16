## Shared setup sourced by every numbered script: libraries, config, data
## loading, and helpers/variables reused across more than one script.

library(here)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(cowplot)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(FSA)
library(patchwork)
library(ggpubr)
library(indicspecies)
library(phyloseq)
library(microshades)


set.seed(123)

config <- config::get()

data_folder   <- config$input_folder
output_folder <- config$output_folder
images_out    <- config$images_output_folder
data_out      <- file.path(output_folder, "data")

metadata_file <- config$metadata_file
otu_file      <- config$otu_file
tax_file      <- config$tax_file

group_levels     <- config$sample_groups
ecosystem_levels <- config$ecosystems

rarecurve_step             <- config$analysis$rarecurve_step
nmds_k                     <- config$analysis$nmds_k
nmds_trymax                <- config$analysis$nmds_trymax
top_n_taxa                 <- config$analysis$top_n_taxa
p_value_threshold          <- config$analysis$p_value_threshold
indicator_nperm            <- config$analysis$indicator_nperm
indicator_stat_threshold   <- config$analysis$indicator_stat_threshold
indicator_top_n_per_group  <- config$analysis$indicator_top_n_per_group
indicator_heatmap_top_n    <- config$analysis$indicator_heatmap_top_n

dir.create(here(output_folder), showWarnings = FALSE, recursive = TRUE)
dir.create(here(data_out), showWarnings = FALSE, recursive = TRUE)
dir.create(here(images_out), showWarnings = FALSE, recursive = TRUE)

## ---- functions used in more than one script ----

# Assign each sample to one of the four sample groups from its ID.
classify_group <- function(sample_id) {
  case_when(
    grepl("S1_bulk", sample_id) ~ "S1_bulk",
    grepl("S1_rAo", sample_id)  ~ "S1_rAo",
    grepl("S2_bulk", sample_id) ~ "S2_bulk",
    grepl("S2_rVo", sample_id)  ~ "S2_rVo"
  )
}

# Assign each sample to its ecosystem (soil vs. rhizosphere) from its ID.
classify_ecosystem <- function(sample_id) {
  case_when(
    grepl("bulk", sample_id) ~ "Soil",
    grepl("rAo|rVo", sample_id) ~ "Rhizosphere"
  )
}

# Assign each sample to its host plant from its ID (bulk soil has no plant
# of its own, so it's labeled by the site's rhizosphere plant, matching the
# convention already used for beta-diversity's sample metadata).
classify_plant <- function(sample_id) {
  ifelse(grepl("rAo", sample_id), "Astragalus", "Vicia")
}

# Convert an OTU count table into per-sample relative abundance.
normalize_otu <- function(otu) {
  sweep(otu, 2, colSums(otu), "/")
}

# Sum OTU (relative) abundance up to an arbitrary taxonomic label vector.
aggregate_by_taxon <- function(data_count, labels) {
  data_agg <- cbind(data_count, Taxon = labels)
  data_agg <- aggregate(. ~ Taxon, data_agg, sum)
  rownames(data_agg) <- data_agg[, 1]
  data_agg[, -1]
}

## ---- OTU / taxonomy tables ----

otu_df <- read.csv(
  here(data_folder, otu_file),
  header = TRUE, row.names = 1, check.names = FALSE, sep = "\t"
)
tax_df <- read.csv(
  here(data_folder, tax_file),
  header = TRUE, row.names = 1, check.names = FALSE, sep = "\t"
)

sample_names <- colnames(otu_df)
t_otu_df <- as.data.frame(t(otu_df))

## ---- metadata (environmental metadata, keyed by the "Sample" column) ----

metadata <- read.csv(
  here(data_folder, metadata_file),
  header = TRUE, check.names = FALSE, sep = ","
)

metadata <- metadata[match(sample_names, metadata$Sample), ]
rownames(metadata) <- metadata$Sample
metadata$Ecosystem <- classify_ecosystem(metadata$Sample)

## ---- shared plotting palettes ----

group_colors <- setNames(brewer.pal(length(group_levels), "Set1"), group_levels)
ecosystem_colors <- setNames(c("#E69F00", "#56B4E9"), ecosystem_levels)
