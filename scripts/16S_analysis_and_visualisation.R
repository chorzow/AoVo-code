plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

readSampleData = function(sampledataFile, sep = "\t")
{
  metadata <- read.csv(sampledataFile, sep = sep)
  meta.df <- metadata %>% select(-"SampleID") %>% as.data.frame
  rownames(meta.df) <- metadata$SampleID
  return (sample_data(meta.df))
}

my_plot_bar = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL) 
{
  # This part is taken from https://github.com/joey711/phyloseq/issues/938#issuecomment-390420586
  p <- plot_bar(physeq, x = x, y = y, fill = fill)
  pd <- p$data %>% as_tibble %>% mutate(fill = as.character(fill))# %>%  replace_na(list(fill = "unknown"))
  fill_abun <- pd %>% group_by(fill[1]) %>% summarize(Abundance = sum(Abundance)) %>% arrange(Abundance)
  fill_levels <- fill_abun$fill
  pd0 <- pd %>% mutate(fill = factor(fill, fill_levels))
  p = ggplot(pd0, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  #p = p + theme_set(theme_bw())
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p = p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p = p + ggtitle(title)
  }
  return(p)
}



#library("dada2"); packageVersion("dada2")
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("magrittr") # For operator %>%
library("dplyr") # For "select" function
library("vegan"); packageVersion("vegan")
library("ranacapa")#ggrare
library("svglite")#save svg
#library("MicrobeR") #PCoA3D
library("ape")
library(tidyverse)
library(phyloseq)
library('microbiomeMarker'); packageVersion("microbiomeMarker")
library(knitr)
library(microbiome)


# Path to output directory.
Figures_outpath <- "C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\MetaRUS\\Collaborations\\SIFIBR_RAS\\16S_data_analysis\\97_clust_SILVA\\Figures\\"

# Make a Phyloseq object.
# Import frequency.
seqtab.nochim_OTU <- read.table(sep = "\t", file = "C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\MetaRUS\\Collaborations\\SIFIBR_RAS\\16S_data_analysis\\97_clust_SILVA\\all_OTU_frequency_97_clean.tsv", header = TRUE, check.names=FALSE, row.names=1)
seqtab.nochim_OTU <- t(seqtab.nochim_OTU)
# Import phylogeny.
taxa_OTU <- read.table(sep = "\t", file = "C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\MetaRUS\\Collaborations\\SIFIBR_RAS\\16S_data_analysis\\97_clust_SILVA\\all_OTUs_phylogeny_97_clean.tsv", header = TRUE, check.names=FALSE, row.names=1)
taxa_OTU <- as.matrix(taxa_OTU)
# Import metadata.
sample.data <- readSampleData("C:\\Users\\sutor\\OneDrive\\ThinkPad_working\\Sutor\\Science\\MetaRUS\\Collaborations\\SIFIBR_RAS\\16S_data_analysis\\97_clust_SILVA\\Metadata_clean_re_ordered.txt")
# Creating Phyloseq object.
ps <- phyloseq(otu_table(seqtab.nochim_OTU, taxa_are_rows = FALSE), sample_data(sample.data), tax_table(taxa_OTU))
dim(otu_table(ps))

table(sample_data(ps)$Location)
table(sample_data(ps)$Biome)
table(sample_data(ps)$Sample.type)
table(sample_data(ps)$Host.organism)


# Set custom pallets.
sponge_pallet_0 <- c("darkolivegreen3", "deepskyblue4")
sponge_pallet_1 <- c("darkgoldenrod2", "deepskyblue3", "darkkhaki", "indianred2", "deepskyblue4")
sponge_pallet_2 <- c("indianred2", "deepskyblue4")
sponge_pallet_3 <- c("cadetblue", "coral2", "darkseagreen", "darkgoldenrod2", "deepskyblue4", "darkolivegreen3", "bisque2", "cornflowerblue", "darkgoldenrod", "lightpink3")

# Rarefaction plot for all samples.
rare_gg_all_wes <- ggrare(ps, step = 100, color = "Sample.type", se = FALSE,  plot = TRUE) + theme_classic()
rare_gg_all_wes <- rare_gg_all_wes + scale_color_manual(values = sponge_pallet_1) + geom_line(size=1.1)
png(file=paste(Figures_outpath, "Rarefaction_curves_all_samples.png", sep=''), width=1000, height=867)
rare_gg_all_wes
dev.off()
svg(file=paste(Figures_outpath, "Rarefaction_curves_all_samples.svg", sep=''), width=1000, height=867)
rare_gg_all_wes
dev.off()

# Select leaf and control samples only.
ps_leaf_control <- subset_samples(ps, Sample.type %in% c("leaf", "laboratory water"))

# Rarefaction plot for leaf and control samples.
rare_gg_all_wes <- ggrare(ps_leaf_control, step = 100, color = "Sample.type", se = FALSE,  plot = TRUE) +
  theme_classic()
rare_gg_all_wes <- rare_gg_all_wes +scale_color_manual(values = sponge_pallet_1) +
  geom_line(size=1.1)
png(file=paste(Figures_outpath, "Rarefaction_curves_leaf_control_samples.png", sep=''), width=1000, height=1135)
rare_gg_all_wes
dev.off()
svg(file=paste(Figures_outpath, "Rarefaction_curves_leaf_control_samples.svg", sep=''), width=1000, height=1135)
rare_gg_all_wes
dev.off()

# Calculate bray-curtis distance matrix.
erie_bray <- phyloseq::distance(ps, method = "bray")

# Make a dataframe from the sample_data.
sampledf <- data.frame(sample_data(ps))

# Adonis test.
adonis_res <- adonis2(erie_bray ~ Sample.type, data = sampledf)
adonis_res

# Calculate Unifrac distancies.
random_tree <- rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
ps_unifrac <- merge_phyloseq(ps, sample.data, random_tree)

# Plot PCoA with unifrac distances.
ord_unifrac <- ordinate(ps_unifrac, "PCoA", "unifrac", weighted=TRUE)
PCoA_UW_all_samples_1_2 <- plot_ordination(ps_unifrac, ord_unifrac, color="Biome", shape="Location", axes=1:2) +
                           scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                           theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_all_samples_1_2.png", sep=''), width=650, height=500)
PCoA_UW_all_samples_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_all_samples_1_2.svg", sep=''), width=650, height=500)
PCoA_UW_all_samples_1_2
dev.off()

PCoA_UW_all_samples_2_3 <- plot_ordination(ps_unifrac, ord_unifrac, color="Biome", shape="Location", axes=2:3) +
                           scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                           theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_all_samples_2_3.png", sep=''), width=650, height=500)
PCoA_UW_all_samples_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_all_samples_2_3.svg", sep=''), width=650, height=500)
PCoA_UW_all_samples_2_3
dev.off()

PCoA_UW_all_samples_1_2_labels <- plot_ordination(ps_unifrac, ord_unifrac, color="Biome", shape="Location", label="F_reads_name", axes=1:2) +
                                  scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                  theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_all_samples_1_2_labels.png", sep=''), width=650, height=500)
PCoA_UW_all_samples_1_2_labels
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_all_samples_1_2_labels.svg", sep=''), width=650, height=500)
PCoA_UW_all_samples_1_2_labels
dev.off()



# Plot PCoA with Bray-Curtis distances.
ord_bc <- ordinate(ps, "PCoA", "bray")
PCoA_BC_all_samples_1_2 <- plot_ordination(ps, ord_bc, color="Biome", shape="Location", axes=1:2) +
                           scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                           theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_all_samples_1_2.png", sep=''), width=650, height=500)
PCoA_BC_all_samples_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_all_samples_1_2.svg", sep=''), width=650, height=500)
PCoA_BC_all_samples_1_2
dev.off()

PCoA_BC_all_samples_2_3 <- plot_ordination(ps, ord_bc, color="Biome", shape="Location", axes=2:3) +
                           scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                           theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_all_samples_2_3.png", sep=''), width=650, height=500)
PCoA_BC_all_samples_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_all_samples_2_3.svg", sep=''), width=650, height=500)
PCoA_BC_all_samples_2_3
dev.off()

PCoA_BC_all_samples_2_3_labels <- plot_ordination(ps, ord_bc, color="Biome", shape="Location", label="F_reads_name", axes=2:3) +
                                  scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                  theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_all_samples_2_3_labels.png", sep=''), width=650, height=500)
PCoA_BC_all_samples_2_3_labels
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_all_samples_2_3_labels.svg", sep=''), width=650, height=500)
PCoA_BC_all_samples_2_3_labels
dev.off()



## Select soil and rhizosphere samples only.
table(sample_data(ps)$Sample.type)
ps_soil_rhizo <- subset_samples(ps, Sample.type %in% c("soil", "rhizosphere"))
dim(otu_table(ps_soil_rhizo))

# Calculate bray-curtis distance matrix.
erie_bray_sr <- phyloseq::distance(ps_soil_rhizo, method = "bray")

# Make a dataframe from the sample_data.
sampledf_sr <- data.frame(sample_data(ps_soil_rhizo))

# Adonis test.
adonis_res_sr <- adonis2(erie_bray_sr ~ Sample.type, data = sampledf_sr)
adonis_res_sr

# Calculate Unifrac distancies.
random_tree_soil_rhizo <- rtree(ntaxa(ps_soil_rhizo), rooted=TRUE, tip.label=taxa_names(ps_soil_rhizo))
ps_soil_rhizo_unifrac <- merge_phyloseq(ps_soil_rhizo, sample.data, random_tree_soil_rhizo)

# Plot PCoA with unifrac distances.
ord_unifrac_soil_rhizo <- ordinate(ps_soil_rhizo_unifrac, "PCoA", "unifrac", weighted=TRUE)
PCoA_UW_soil_rhizosphere_1_2 <- plot_ordination(ps_soil_rhizo_unifrac, ord_unifrac_soil_rhizo, color="Location", shape="Sample.type", axes=1:2) +
                                scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_1_2.png", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_1_2.svg", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_1_2
dev.off()

PCoA_UW_soil_rhizosphere_2_3 <- plot_ordination(ps_soil_rhizo_unifrac, ord_unifrac_soil_rhizo, color="Location", shape="Sample.type", axes=2:3) +
  scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
  theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_2_3.png", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_2_3.svg", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_2_3
dev.off()

PCoA_UW_soil_rhizosphere_1_2_labels <- plot_ordination(ps_soil_rhizo_unifrac, ord_unifrac_soil_rhizo, color="Location", shape="Sample.type", label="F_reads_name", axes=1:2) +
                                       scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                       theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_1_2_labels.png", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_1_2_labels
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_1_2_labels.svg", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_1_2_labels
dev.off()


# Plot PCoA with Bray-Curtis distances.
ord_bc_soil_rhizo <- ordinate(ps_soil_rhizo, "PCoA", "bray")
PCoA_BC_soil_rhizosphere_1_2 <- plot_ordination(ps_soil_rhizo, ord_bc_soil_rhizo, color="Location", shape="Sample.type", axes=1:2) +
                                scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_1_2.png", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_1_2.svg", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_1_2
dev.off()

PCoA_BC_soil_rhizosphere_2_3 <- plot_ordination(ps_soil_rhizo, ord_bc_soil_rhizo, color="Location", shape="Sample.type", axes=2:3) +
                                scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_2_3.png", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_2_3.svg", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_2_3
dev.off()

PCoA_BC_soil_rhizosphere_2_3_labels <- plot_ordination(ps_soil_rhizo, ord_bc_soil_rhizo, color="Location", shape="Sample.type", label="F_reads_name", axes=2:3) +
                                       scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                       theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_2_3_labels.png", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_2_3_labels
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_2_3_labels.svg", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_2_3_labels
dev.off()



## Select soil and rhizosphere samples only and remove outlyer sample 14447.
ps_soil_rhizo_clean <- subset_samples(ps_soil_rhizo, F_reads_name!="14447_S158_L001_R1_001.fastq.gz")
dim(otu_table(ps_soil_rhizo_clean))

# Calculate bray-curtis distance matrix.
erie_bray_src <- phyloseq::distance(ps_soil_rhizo_clean, method = "bray")

# Make a dataframe from the sample_data.
sampledf_src <- data.frame(sample_data(ps_soil_rhizo_clean))

# Adonis test.
adonis_res_src <- adonis2(erie_bray_src ~ Sample.type, data = sampledf_src)
adonis_res_src

# Calculate Unifrac distancies.
random_tree_soil_rhizo_clean <- rtree(ntaxa(ps_soil_rhizo_clean), rooted=TRUE, tip.label=taxa_names(ps_soil_rhizo_clean))
ps_soil_rhizo_clean_unifrac <- merge_phyloseq(ps_soil_rhizo_clean, sample.data, random_tree_soil_rhizo_clean)

# Plot PCoA with unifrac distances.
ord_unifrac_soil_rhizo_clean <- ordinate(ps_soil_rhizo_clean_unifrac, "PCoA", "unifrac", weighted=TRUE)
PCoA_UW_soil_rhizosphere_clean_1_2 <- plot_ordination(ps_soil_rhizo_clean_unifrac, ord_unifrac_soil_rhizo_clean, color="Location", shape="Sample.type", axes=1:2) +
                                      scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                      theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_1_2.png", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_clean_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_1_2.svg", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_clean_1_2
dev.off()

PCoA_UW_soil_rhizosphere_clean_2_3 <- plot_ordination(ps_soil_rhizo_clean_unifrac, ord_unifrac_soil_rhizo_clean, color="Location", shape="Sample.type", axes=2:3) +
                                      scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                      theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_2_3.png", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_clean_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_2_3.svg", sep=''), width=650, height=500)
PCoA_UW_soil_rhizosphere_clean_2_3
dev.off()

# Make station number a discrete variable.
station_name_upd <- paste("Station", ps_soil_rhizo_clean_unifrac@sam_data[["Expedition.station.number"]])
sample_data(ps_soil_rhizo_clean_unifrac)$Expedition.station.number_upd <- station_name_upd

PCoA_UW_soil_rhizosphere_clean_stations_1_2 <- plot_ordination(ps_soil_rhizo_clean_unifrac, ord_unifrac_soil_rhizo_clean, color="Expedition.station.number_upd", shape="Sample.type", axes=1:2) +
                                               scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen", "darkgoldenrod2", "deepskyblue4", "darkolivegreen3")) + geom_point(size=5, alpha=0.75) +
                                               theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_stations_1_2.png", sep=''), width=730, height=500)
PCoA_UW_soil_rhizosphere_clean_stations_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_stations_1_2.svg", sep=''), width=730, height=500)
PCoA_UW_soil_rhizosphere_clean_stations_1_2
dev.off()

PCoA_UW_soil_rhizosphere_clean_stations_2_3 <- plot_ordination(ps_soil_rhizo_clean_unifrac, ord_unifrac_soil_rhizo_clean, color="Expedition.station.number_upd", shape="Sample.type", axes=2:3) +
                                               scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen", "darkgoldenrod2", "deepskyblue4", "darkolivegreen3")) + geom_point(size=5, alpha=0.75) +
                                               theme_bw()
png(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_stations_2_3.png", sep=''), width=730, height=500)
PCoA_UW_soil_rhizosphere_clean_stations_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_UW_soil_rhizosphere_clean_stations_2_3.svg", sep=''), width=730, height=500)
PCoA_UW_soil_rhizosphere_clean_stations_2_3
dev.off()

# Plot PCoA with Bray-Curtis distances.
ord_bc_soil_rhizo_clean <- ordinate(ps_soil_rhizo_clean, "PCoA", "bray")
PCoA_BC_soil_rhizosphere_clean_1_2 <- plot_ordination(ps_soil_rhizo_clean, ord_bc_soil_rhizo_clean, color="Location", shape="Sample.type", axes=1:2) +
                                      scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                      theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_1_2.png", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_clean_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_1_2.svg", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_clean_1_2
dev.off()

PCoA_BC_soil_rhizosphere_clean_2_3 <- plot_ordination(ps_soil_rhizo_clean, ord_bc_soil_rhizo_clean, color="Location", shape="Sample.type", axes=2:3) +
                                      scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen")) + geom_point(size=5, alpha=0.75) +
                                      theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_2_3.png", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_clean_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_2_3.svg", sep=''), width=650, height=500)
PCoA_BC_soil_rhizosphere_clean_2_3
dev.off()

# Make station number a discrete variable.
station_name_upd <- paste("Station", ps_soil_rhizo_clean@sam_data[["Expedition.station.number"]])
sample_data(ps_soil_rhizo_clean)$Expedition.station.number_upd <- station_name_upd

PCoA_BC_soil_rhizosphere_clean_stations_1_2 <- plot_ordination(ps_soil_rhizo_clean, ord_bc_soil_rhizo_clean, color="Expedition.station.number_upd", shape="Sample.type", axes=1:2) +
                                               scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen", "darkgoldenrod2", "deepskyblue4", "darkolivegreen3")) + geom_point(size=5, alpha=0.75) +
                                               theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_stations_1_2.png", sep=''), width=730, height=500)
PCoA_BC_soil_rhizosphere_clean_stations_1_2
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_stations_1_2.svg", sep=''), width=730, height=500)
PCoA_BC_soil_rhizosphere_clean_stations_1_2
dev.off()

PCoA_BC_soil_rhizosphere_clean_stations_2_3 <- plot_ordination(ps_soil_rhizo_clean, ord_bc_soil_rhizo_clean, color="Expedition.station.number_upd", shape="Sample.type", axes=2:3) +
                                               scale_colour_manual(values = c("cadetblue", "coral2", "darkseagreen", "darkgoldenrod2", "deepskyblue4", "darkolivegreen3")) + geom_point(size=5, alpha=0.75) +
                                               theme_bw()
png(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_stations_2_3.png", sep=''), width=730, height=500)
PCoA_BC_soil_rhizosphere_clean_stations_2_3
dev.off()
svg(file=paste(Figures_outpath, "PCoA_BC_soil_rhizosphere_clean_stations_2_3.svg", sep=''), width=730, height=500)
PCoA_BC_soil_rhizosphere_clean_stations_2_3
dev.off()



# Perform LefSe analysis.
# Compare composition of datasets.
# Compare general rhizosphere vs soil samples.
lef_out <- run_lefse(ps_soil_rhizo_clean, group = "Sample.type", norm = "CPM", taxa_rank = "none",
                     kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_all_samples.png", sep=''), width=500, height=450)
plot_ef_bar(lef_out)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_all_samples.svg", sep=''), width=500, height=450)
plot_ef_bar(lef_out)
dev.off()

lef_out_OTUs_list <- lef_out@marker_table$feature
lef_full_taxonomy_tt <- ps@tax_table[lef_out_OTUs_list, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df <- data.frame(as(tax_table(lef_full_taxonomy_tt), "matrix"))
lef_full_taxonomy_df$Full_taxonomy <- paste(lef_full_taxonomy_df$Order, lef_full_taxonomy_df$Family, lef_full_taxonomy_df$Genus, rownames(lef_full_taxonomy_df), sep='/')
lef_full_taxonomy_df$Full_taxonomy
lef_out@marker_table$feature <- lef_full_taxonomy_df$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_all_samples_full_names.png", sep=''), width=1000, height=450)
plot_ef_bar(lef_out)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_all_samples_full_names.svg", sep=''), width=1000, height=450)
plot_ef_bar(lef_out)
dev.off()

#Compare Astragalus olchonensis rhizosphere and soil at Station 1.
table(sample_data(ps_soil_rhizo_clean)$Host.organism)
table(sample_data(ps_soil_rhizo_clean)$Location)
table(sample_data(ps_soil_rhizo_clean)$Expedition.station.number_upd)

ps_soil_rhizo_clean_st_1 <- subset_samples(ps_soil_rhizo_clean, Expedition.station.number_upd=="Station 1")
table(sample_data(ps_soil_rhizo_clean_st_1)$Host.organism)

lef_out_st_1 <- run_lefse(ps_soil_rhizo_clean_st_1, group = "Host.organism", norm = "CPM", taxa_rank = "none",
                     kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_1_Astr_olchonensis.png", sep=''), width=500, height=450)
plot_ef_bar(lef_out_st_1)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_1_Astr_olchonensis.svg", sep=''), width=500, height=450)
plot_ef_bar(lef_out_st_1)
dev.off()

lef_out_OTUs_list_st_1 <- lef_out_st_1@marker_table$feature
lef_full_taxonomy_tt_st_1 <- ps@tax_table[lef_out_OTUs_list_st_1, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df_st_1 <- data.frame(as(tax_table(lef_full_taxonomy_tt_st_1), "matrix"))
lef_full_taxonomy_df_st_1$Full_taxonomy <- paste(lef_full_taxonomy_df_st_1$Order, lef_full_taxonomy_df_st_1$Family, lef_full_taxonomy_df_st_1$Genus, rownames(lef_full_taxonomy_df_st_1), sep='/')
lef_full_taxonomy_df_st_1$Full_taxonomy
lef_out_st_1@marker_table$feature <- lef_full_taxonomy_df_st_1$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_1_Astr_olchonensis_full_names.png", sep=''), width=1000, height=450)
plot_ef_bar(lef_out_st_1)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_1_Astr_olchonensis_full_names.svg", sep=''), width=1000, height=450)
plot_ef_bar(lef_out_st_1)
dev.off()

#Compare Vicia olchonensis rhizosphere and soil at Station 2.
ps_soil_rhizo_clean_st_2 <- subset_samples(ps_soil_rhizo_clean, Expedition.station.number_upd=="Station 2")
table(sample_data(ps_soil_rhizo_clean_st_2)$Host.organism)

lef_out_st_2 <- run_lefse(ps_soil_rhizo_clean_st_2, group = "Host.organism", norm = "CPM", taxa_rank = "none",
                          kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_2_Vicia_olchonensis.png", sep=''), width=500, height=1200)
plot_ef_bar(lef_out_st_2)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_2_Vicia_olchonensis.svg", sep=''), width=500, height=1200)
plot_ef_bar(lef_out_st_2)
dev.off()

lef_out_OTUs_list_st_2 <- lef_out_st_2@marker_table$feature
lef_full_taxonomy_tt_st_2 <- ps@tax_table[lef_out_OTUs_list_st_2, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df_st_2 <- data.frame(as(tax_table(lef_full_taxonomy_tt_st_2), "matrix"))
lef_full_taxonomy_df_st_2$Full_taxonomy <- paste(lef_full_taxonomy_df_st_2$Order, lef_full_taxonomy_df_st_2$Family, lef_full_taxonomy_df_st_2$Genus, rownames(lef_full_taxonomy_df_st_2), sep='/')
lef_full_taxonomy_df_st_2$Full_taxonomy
lef_out_st_2@marker_table$feature <- lef_full_taxonomy_df_st_2$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_2_Vicia_olchonensis_full_names.png", sep=''), width=1000, height=1200)
plot_ef_bar(lef_out_st_2)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_2_Vicia_olchonensis_full_names.svg", sep=''), width=1000, height=1200)
plot_ef_bar(lef_out_st_2)
dev.off()

#Compare Hedysarum zundukii rhizosphere and soil at Station 3.
ps_soil_rhizo_clean_st_3 <- subset_samples(ps_soil_rhizo_clean, Expedition.station.number_upd=="Station 3")
table(sample_data(ps_soil_rhizo_clean_st_3)$Host.organism)

lef_out_st_3 <- run_lefse(ps_soil_rhizo_clean_st_3, group = "Host.organism", norm = "CPM", taxa_rank = "none",
                          kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_3_Hedysarum_zundukii.png", sep=''), width=500, height=600)
plot_ef_bar(lef_out_st_3)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_3_Hedysarum_zundukii.svg", sep=''), width=500, height=600)
plot_ef_bar(lef_out_st_3)
dev.off()

lef_out_OTUs_list_st_3 <- lef_out_st_3@marker_table$feature
lef_full_taxonomy_tt_st_3 <- ps@tax_table[lef_out_OTUs_list_st_3, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df_st_3 <- data.frame(as(tax_table(lef_full_taxonomy_tt_st_3), "matrix"))
lef_full_taxonomy_df_st_3$Full_taxonomy <- paste(lef_full_taxonomy_df_st_3$Order, lef_full_taxonomy_df_st_3$Family, lef_full_taxonomy_df_st_3$Genus, rownames(lef_full_taxonomy_df_st_3), sep='/')
lef_full_taxonomy_df_st_3$Full_taxonomy
lef_out_st_3@marker_table$feature <- lef_full_taxonomy_df_st_3$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_3_Hedysarum_zundukii_full_names.png", sep=''), width=1000, height=600)
plot_ef_bar(lef_out_st_3)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_3_Hedysarum_zundukii_full_names.svg", sep=''), width=1000, height=600)
plot_ef_bar(lef_out_st_3)
dev.off()

#Compare Oxytropis triphylla, Astragalus chorinensis rhizosphere and soil at Station 4.
ps_soil_rhizo_clean_st_4 <- subset_samples(ps_soil_rhizo_clean, Expedition.station.number_upd=="Station 4")
table(sample_data(ps_soil_rhizo_clean_st_4)$Host.organism)

lef_out_st_4 <- run_lefse(ps_soil_rhizo_clean_st_4, group = "Host.organism", norm = "CPM", taxa_rank = "none",
                          kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_4_Astr_chorinensis_Oxy_triphylla.png", sep=''), width=500, height=800)
plot_ef_bar(lef_out_st_4)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_4_Astr_chorinensis_Oxy_triphylla.svg", sep=''), width=500, height=800)
plot_ef_bar(lef_out_st_4)
dev.off()

lef_out_OTUs_list_st_4 <- lef_out_st_4@marker_table$feature
lef_full_taxonomy_tt_st_4 <- ps@tax_table[lef_out_OTUs_list_st_4, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df_st_4 <- data.frame(as(tax_table(lef_full_taxonomy_tt_st_4), "matrix"))
lef_full_taxonomy_df_st_4$Full_taxonomy <- paste(lef_full_taxonomy_df_st_4$Order, lef_full_taxonomy_df_st_4$Family, lef_full_taxonomy_df_st_4$Genus, rownames(lef_full_taxonomy_df_st_4), sep='/')
lef_full_taxonomy_df_st_4$Full_taxonomy
lef_out_st_4@marker_table$feature <- lef_full_taxonomy_df_st_4$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_4_Astr_chorinensis_Oxy_triphylla_full_names.png", sep=''), width=1000, height=800)
plot_ef_bar(lef_out_st_4)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_4_Astr_chorinensis_Oxy_triphylla_full_names.svg", sep=''), width=1000, height=800)
plot_ef_bar(lef_out_st_4)
dev.off()

#Compare Oxytropis triphylla rhizosphere and soil at Station 5.
ps_soil_rhizo_clean_st_5 <- subset_samples(ps_soil_rhizo_clean, Expedition.station.number_upd=="Station 5")
table(sample_data(ps_soil_rhizo_clean_st_5)$Host.organism)

lef_out_st_5 <- run_lefse(ps_soil_rhizo_clean_st_5, group = "Host.organism", norm = "CPM", taxa_rank = "none",
                          kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_5_Oxy_triphylla.png", sep=''), width=500, height=600)
plot_ef_bar(lef_out_st_5)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_5_Oxy_triphylla.svg", sep=''), width=500, height=600)
plot_ef_bar(lef_out_st_5)
dev.off()

lef_out_OTUs_list_st_5 <- lef_out_st_5@marker_table$feature
lef_full_taxonomy_tt_st_5 <- ps@tax_table[lef_out_OTUs_list_st_5, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df_st_5 <- data.frame(as(tax_table(lef_full_taxonomy_tt_st_5), "matrix"))
lef_full_taxonomy_df_st_5$Full_taxonomy <- paste(lef_full_taxonomy_df_st_5$Order, lef_full_taxonomy_df_st_5$Family, lef_full_taxonomy_df_st_5$Genus, rownames(lef_full_taxonomy_df_st_5), sep='/')
lef_full_taxonomy_df_st_5$Full_taxonomy
lef_out_st_5@marker_table$feature <- lef_full_taxonomy_df_st_5$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_5_Oxy_triphylla_full_names.png", sep=''), width=1000, height=600)
plot_ef_bar(lef_out_st_5)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_5_Oxy_triphylla_full_names.svg", sep=''), width=1000, height=600)
plot_ef_bar(lef_out_st_5)
dev.off()

#Compare Oxytropis triphylla rhizosphere and soil at Station 6.
ps_soil_rhizo_clean_st_6 <- subset_samples(ps_soil_rhizo_clean, Expedition.station.number_upd=="Station 6")
table(sample_data(ps_soil_rhizo_clean_st_6)$Host.organism)

lef_out_st_6 <- run_lefse(ps_soil_rhizo_clean_st_6, group = "Host.organism", norm = "CPM", taxa_rank = "none",
                          kw_cutoff = 0.01, lda_cutoff = 2.5)

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_6_Oxy_triphylla_Oxy_microphylla.png", sep=''), width=500, height=2000)
plot_ef_bar(lef_out_st_6)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_6_Oxy_triphylla_Oxy_microphylla.svg", sep=''), width=500, height=2000)
plot_ef_bar(lef_out_st_6)
dev.off()

lef_out_OTUs_list_st_6 <- lef_out_st_6@marker_table$feature
lef_full_taxonomy_tt_st_6 <- ps@tax_table[lef_out_OTUs_list_st_6, c('Order', 'Family', 'Genus')]
lef_full_taxonomy_df_st_6 <- data.frame(as(tax_table(lef_full_taxonomy_tt_st_6), "matrix"))
lef_full_taxonomy_df_st_6$Full_taxonomy <- paste(lef_full_taxonomy_df_st_6$Order, lef_full_taxonomy_df_st_6$Family, lef_full_taxonomy_df_st_6$Genus, rownames(lef_full_taxonomy_df_st_6), sep='/')
lef_full_taxonomy_df_st_6$Full_taxonomy
lef_out_st_6@marker_table$feature <- lef_full_taxonomy_df_st_6$Full_taxonomy

png(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_6_Oxy_triphylla_Oxy_microphylla_full_names.png", sep=''), width=1000, height=2000)
plot_ef_bar(lef_out_st_6)
dev.off()
svg(file=paste(Figures_outpath, "LefSe_soil_rhizosphere_clean_st_6_Oxy_triphylla_Oxy_microphylla_full_names.svg", sep=''), width=1000, height=2000)
plot_ef_bar(lef_out_st_6)
dev.off()


# Plot alpha-diversity metrics.
plot_richness(ps_soil_rhizo_clean, x="Sample.type", measures=c("Chao1", "Shannon", "Simpson")) + 
  theme_set(theme_bw()) +
  geom_point(size=5)

# Plot Shannon index.
# Group by sample type.
p.shannon_st <- boxplot_alpha(ps_soil_rhizo_clean, 
                           index = c("shannon"),
                           x_var = "Sample.type",
                           fill.colors = c("darkgoldenrod2","deepskyblue3"))

png(file=paste(Figures_outpath, "Shannon_index_soil_rhizosphere_clean_all_samples.png", sep=''), width=500, height=450)
p.shannon_st
dev.off()
svg(file=paste(Figures_outpath, "Shannon_index_soil_rhizosphere_clean_all_samples.svg", sep=''), width=500, height=450)
p.shannon_st
dev.off()

# Group by location.
p.shannon_loc <- boxplot_alpha(ps_soil_rhizo_clean, 
                           index = c("shannon"),
                           x_var = "Location",
                           fill.colors = c("darkgoldenrod2", "deepskyblue3", "darkkhaki"))

png(file=paste(Figures_outpath, "Shannon_index_Location_clean_all_samples.png", sep=''), width=700, height=600)
p.shannon_loc
dev.off()
svg(file=paste(Figures_outpath, "Shannon_index_Location_clean_all_samples.svg", sep=''), width=700, height=600)
p.shannon_loc
dev.off()

# Group by species name.
p.shannon_ho <- boxplot_alpha(ps_soil_rhizo_clean, 
                              index = c("shannon"),
                              x_var = "Host.organism",
                              fill.colors = c("cadetblue", "coral2", "darkseagreen", "darkgoldenrod2", "deepskyblue4", "darkolivegreen3", "bisque2"))
p.shannon_ho <- p.shannon_ho + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png(file=paste(Figures_outpath, "Shannon_index_organism_clean_all_samples.png", sep=''), width=500, height=400)
p.shannon_ho
dev.off()
svg(file=paste(Figures_outpath, "Shannon_index_organism_clean_all_samples.svg", sep=''), width=500, height=400)
p.shannon_ho
dev.off()

# Plot Chao1 index.
p.chao1 <- boxplot_alpha(ps_soil_rhizo_clean, 
                           index = c("chao1"),
                           x_var = "Sample.type",
                           fill.colors = c("darkgoldenrod2","deepskyblue3"))

png(file=paste(Figures_outpath, "Chao1_index_soil_rhizosphere_clean_all_samples.png", sep=''), width=500, height=450)
p.chao1
dev.off()
svg(file=paste(Figures_outpath, "Chao1_index_soil_rhizosphere_clean_all_samples.svg", sep=''), width=500, height=450)
p.chao1
dev.off()














#ADD BARPLOTS.
ps
ps_soil_rhizo
ps_soil_rhizo_clean
ps_soil_rhizo_clean_st_1
ps_soil_rhizo_clean_st_2
ps_soil_rhizo_clean_st_3
ps_soil_rhizo_clean_st_4
ps_soil_rhizo_clean_st_5
ps_soil_rhizo_clean_st_6

#top 20 Family ps
genus.sum <- tapply(taxa_sums(ps), tax_table(ps)[, "Family"], sum, na.rm=TRUE)
top20genus <- names(sort(genus.sum, TRUE))[1:20]
Genustop20 <- prune_taxa((tax_table(ps)[, "Family"] %in% top20genus), ps)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new <- c("14390", "14396", "14397", "14398", "14399", "14922", "14923", "14924", "14925", "14926", "14400", "14401", "14406", "14407", "14408", "14409", "14927", "14928", "14929", "14930", "14410", "14415", "14416", "14417", "14418", "14931", "14932", "14933", "14934", "14935", "14420", "14421", "14425", "14426", "14427", "14428", "14429", "14936", "14937", "14938", "14430", "14431", "14432", "14433", "14434", "14435", "14436", "14437", "14438", "14439", "14440", "14441", "14442", "14445", "14446", "14447", "14448", "14939", "14940", "14941", "14450", "14451", "14453", "14455", "14456", "14458", "14459", "14942", "14943", "14944", "14652", "14653", "14654", "14655", "14656", "14657", "14658", "14659", "14660", "14661", "14662", "14663", "14664", "14666", "14667", "14668", "14669", "14670", "14945", "14946", "14672", "14673", "14675", "14676", "14678", "14679", "14947", "14948", "14949", "14950", "14682", "14683", "14684", "14688", "14689", "14951", "14952", "14953", "14954", "14955", "14692", "14693", "14698", "14699", "14956", "14957", "14958", "14959", "14960", "14961", "14702", "14703", "14709", "14962", "14963", "14964", "14965", "14966", "14967", "14968", "14712", "14713", "14719", "14720", "14721", "14969", "14970", "14971", "14722", "14731", "14732", "14733", "14734", "14735", "14736", "14737", "14738")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new)

ggsave( file = "All_samples_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 15, height = 10, dpi = 300)
ggsave( file = "All_samples_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 15, height = 10, dpi = 300)


#top 20 Family ps_soil_rhizo
genus.sum = tapply(taxa_sums(ps_soil_rhizo), tax_table(ps_soil_rhizo)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo)[, "Family"] %in% top20genus), ps_soil_rhizo)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new <- c("14390", "14396", "14397", "14398", "14399", "14922", "14923", "14924", "14925", "14926", "14400", "14401", "14406", "14407", "14408", "14409", "14927", "14928", "14929", "14930", "14410", "14415", "14416", "14417", "14418", "14931", "14932", "14933", "14934", "14935", "14420", "14421", "14425", "14426", "14427", "14428", "14429", "14936", "14937", "14938", "14430", "14431", "14432", "14433", "14434", "14435", "14436", "14437", "14438", "14439", "14440", "14441", "14442", "14445", "14446", "14447", "14448", "14939", "14940", "14941", "14450", "14451", "14453", "14455", "14456", "14458", "14459", "14942", "14943", "14944", "14652", "14653", "14654", "14655", "14656", "14657", "14658", "14659", "14660", "14661", "14662", "14663", "14664", "14666", "14667", "14668", "14669", "14670", "14945", "14946", "14672", "14673", "14675", "14676", "14678", "14679", "14947", "14948", "14949", "14950", "14682", "14683", "14684", "14688", "14689", "14951", "14952", "14953", "14954", "14955", "14692", "14693", "14698", "14699", "14956", "14957", "14958", "14959", "14960", "14961", "14702", "14703", "14709", "14962", "14963", "14964", "14965", "14966", "14967", "14968", "14712", "14713", "14719", "14720", "14721", "14969", "14970", "14971")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new)

ggsave( file = "ps_soil_rhizo_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 15, height = 10, dpi = 300)
ggsave( file = "ps_soil_rhizo_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 15, height = 10, dpi = 300)


#top 20 Family ps_soil_rhizo_clean
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean), tax_table(ps_soil_rhizo_clean)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean)[, "Family"] %in% top20genus), ps_soil_rhizo_clean)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new <- c("14390", "14396", "14397", "14398", "14399", "14922", "14923", "14924", "14925", "14926", "14400", "14401", "14406", "14407", "14408", "14409", "14927", "14928", "14929", "14930", "14410", "14415", "14416", "14417", "14418", "14931", "14932", "14933", "14934", "14935", "14420", "14421", "14425", "14426", "14427", "14428", "14429", "14936", "14937", "14938", "14430", "14431", "14432", "14433", "14434", "14435", "14436", "14437", "14438", "14439", "14440", "14441", "14442", "14445", "14446", "14448", "14939", "14940", "14941", "14450", "14451", "14453", "14455", "14456", "14458", "14459", "14942", "14943", "14944", "14652", "14653", "14654", "14655", "14656", "14657", "14658", "14659", "14660", "14661", "14662", "14663", "14664", "14666", "14667", "14668", "14669", "14670", "14945", "14946", "14672", "14673", "14675", "14676", "14678", "14679", "14947", "14948", "14949", "14950", "14682", "14683", "14684", "14688", "14689", "14951", "14952", "14953", "14954", "14955", "14692", "14693", "14698", "14699", "14956", "14957", "14958", "14959", "14960", "14961", "14702", "14703", "14709", "14962", "14963", "14964", "14965", "14966", "14967", "14968", "14712", "14713", "14719", "14720", "14721", "14969", "14970", "14971")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new)

ggsave( file = "ps_soil_rhizo_clean_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 15, height = 10, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 15, height = 10, dpi = 300)


#top 20 Family ps_soil_rhizo_clean_st_1
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean_st_1), tax_table(ps_soil_rhizo_clean_st_1)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean_st_1)[, "Family"] %in% top20genus), ps_soil_rhizo_clean_st_1)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new_st_1 <- c("14390", "14396", "14397", "14398", "14399", "14922", "14923", "14924", "14925", "14926", "14400", "14401", "14406", "14407", "14408", "14409", "14927", "14928", "14929", "14930")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new_st_1)

ggsave( file = "ps_soil_rhizo_clean_st_1_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_st_1_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)


#top 20 Family ps_soil_rhizo_clean_st_2
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean_st_2), tax_table(ps_soil_rhizo_clean_st_2)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean_st_2)[, "Family"] %in% top20genus), ps_soil_rhizo_clean_st_2)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new_st_2 <- c("14410", "14415", "14416", "14417", "14418", "14931", "14932", "14933", "14934", "14935", "14420", "14421", "14425", "14426", "14427", "14428", "14429", "14936", "14937", "14938")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new_st_2)

ggsave( file = "ps_soil_rhizo_clean_st_2_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_st_2_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)


#top 20 Family ps_soil_rhizo_clean_st_3
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean_st_3), tax_table(ps_soil_rhizo_clean_st_3)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean_st_3)[, "Family"] %in% top20genus), ps_soil_rhizo_clean_st_3)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new_st_3 <- c("14430", "14431", "14432", "14433", "14434", "14435", "14436", "14437", "14438", "14439", "14440", "14441", "14442", "14445", "14446", "14448", "14939", "14940", "14941")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new_st_3)

ggsave( file = "ps_soil_rhizo_clean_st_3_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_st_3_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)


#top 20 Family ps_soil_rhizo_clean_st_4
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean_st_4), tax_table(ps_soil_rhizo_clean_st_4)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean_st_4)[, "Family"] %in% top20genus), ps_soil_rhizo_clean_st_4)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new_st_4 <- c("14450", "14451", "14453", "14455", "14456", "14458", "14459", "14942", "14943", "14944", "14652", "14653", "14654", "14655", "14656", "14657", "14658", "14659", "14660", "14661", "14662", "14663", "14664", "14666", "14667", "14668", "14669", "14670", "14945", "14946")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new_st_4)

ggsave( file = "ps_soil_rhizo_clean_st_4_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_st_4_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)


#top 20 Family ps_soil_rhizo_clean_st_5
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean_st_5), tax_table(ps_soil_rhizo_clean_st_5)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean_st_5)[, "Family"] %in% top20genus), ps_soil_rhizo_clean_st_5)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new_st_5 <- c("14672", "14673", "14675", "14676", "14678", "14679", "14947", "14948", "14949", "14950", "14682", "14683", "14684", "14688", "14689", "14951", "14952", "14953", "14954", "14955")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new_st_5)

ggsave( file = "ps_soil_rhizo_clean_st_5_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_st_5_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)


#top 20 Family ps_soil_rhizo_clean_st_6
genus.sum = tapply(taxa_sums(ps_soil_rhizo_clean_st_6), tax_table(ps_soil_rhizo_clean_st_6)[, "Family"], sum, na.rm=TRUE)
top20genus = names(sort(genus.sum, TRUE))[1:20]
Genustop20 = prune_taxa((tax_table(ps_soil_rhizo_clean_st_6)[, "Family"] %in% top20genus), ps_soil_rhizo_clean_st_6)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Family", border_color = NA)
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                               "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
  theme(legend.position="bottom") 
bar_top20_genus <- bar_top20_genus + theme(panel.background = element_blank())
All_sample_order_new_st_6 <- c("14692", "14693", "14698", "14699", "14956", "14957", "14958", "14959", "14960", "14961", "14702", "14703", "14709", "14962", "14963", "14964", "14965", "14966", "14967", "14968", "14712", "14713", "14719", "14720", "14721", "14969", "14970", "14971")
bar_top20_genus$data$Sample <- factor(bar_top20_genus$data$Sample, levels = All_sample_order_new_st_6)

ggsave( file = "ps_soil_rhizo_clean_st_6_barplot_top_20_Family.svg", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)
ggsave( file = "ps_soil_rhizo_clean_st_6_barplot_top_20_Family.png", plot = bar_top20_genus, path = Figures_outpath, width = 10, height = 7, dpi = 300)






# Microbiome marker table
#dat<-marker_table(lef_out) %>% data.frame() %>% select(1:4)
#head(dat)
#write.table(dat, sep="\t", file = "C:/Users/rusla/OneDrive/Desktop/Projects/Hospital/part_1/Run_178/ASV/fetus_size/LefSe_none.tsv")


library(vegan)
library(microbiome)#главный тул здесь 
library(Rtsne)
library(ggplot2)
set.seed(423542)

method <- "tsne"
trans <- "hellinger"
distance <- "euclidean"

# Distance matrix for samples
ps_tsne <- microbiome::transform(ps, trans)

# Calculate sample similarities
dm <- vegdist(otu_table(ps_tsne), distance)

# Run TSNE
tsne_out <- Rtsne(dm, dims = 2) 
proj <- tsne_out$Y
rownames(proj) <- rownames(otu_table(ps_tsne))

p <- plot_landscape(proj, legend = T, size = 2) 
print(p)

#Make TSNE plot

library(tsnemicrobiota)
library(ggplot2)

tsne_res <- tsne_phyloseq(ps1, distance='wunifrac',
                          perplexity = 8, verbose=0, rng_seed = 3901)

plot_tsne_phyloseq(ps1, tsne_res, axes=1:2,
                   color = 'Group', shape="Group", title='t-SNE (Weighted UniFrac)') +
  scale_colour_manual(values = c("darkgoldenrod2", "deepskyblue3")) + geom_point(size=5, alpha=0.75) +
  theme(panel.background = element_blank())




#library("reticulate") - это если плот из python нужно построить в R
#https://datascienceplus.com/how-to-make-seaborn-pairplot-and-heatmap-in-r-write-python-in-r/






###top 10 class 17.11
class.sum = tapply(taxa_sums(ps), tax_table(ps)[, "Class"], sum, na.rm=TRUE)
class.sum 
top10class = names(sort(class.sum, TRUE))[1:10]
Classtop10 = prune_taxa((tax_table(ps)[, "Class"] %in% top10class), ps)
Classtop10 <- prune_samples(sample_sums(Classtop10)>1, Classtop10)
Classtop10_rel <- transform_sample_counts(Classtop10, function(x) x/sum(x))
#plot_bar_2 to remove bar border
bar_top10_class <- plot_bar_2(Classtop10_rel, fill="Class", border_color = NA)
bar_top10_class
bar_top10_class <- bar_top10_class + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkolivegreen3", "deepskyblue4", "darkgoldenrod2", "darkseagreen", "deepskyblue3", "darkkhaki", "bisque2", "cadetblue", "darkgoldenrod", "lightpink3")) +
  theme(legend.position="right") #"cadetblue" #"cornflowerblue"
bar_top10_class <- bar_top10_class + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank())
bar_top10_class + theme(panel.background = element_blank())

###top 15 class
class.sum = tapply(taxa_sums(ps), tax_table(ps)[, "Class"], sum, na.rm=TRUE)
class.sum 
top10class = names(sort(class.sum, TRUE))[1:15]
Classtop10 = prune_taxa((tax_table(ps)[, "Class"] %in% top10class), ps)
Classtop10 <- prune_samples(sample_sums(Classtop10)>1, Classtop10)
Classtop10_rel <- transform_sample_counts(Classtop10, function(x) x/sum(x))
#plot_bar_2 to remove bar border
bar_top10_class <- plot_bar_2(Classtop10_rel, fill="Class", border_color = NA)
bar_top10_class
bar_top10_class <- bar_top10_class + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkolivegreen3", "deepskyblue4", "darkgoldenrod2", "darkseagreen", "deepskyblue3", "darkkhaki", "bisque2", "cadetblue", "darkgoldenrod", "lightpink3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2")) +
  theme(legend.position="right") #"cadetblue" #"cornflowerblue"
bar_top10_class <- bar_top10_class + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank())
bar_top10_class + theme(panel.background = element_blank())
#Classtop10_rel <- transform_sample_counts(Classtop10, function(x) x/sum(x))
#my_plot_bar(Classtop10_rel, fill = "Class")
print(top10class)

###top 10 order
order.sum = tapply(taxa_sums(ps), tax_table(ps)[, "Order"], sum, na.rm=TRUE)
order.sum 
top10order = names(sort(order.sum, TRUE))[1:10]
Ordertop10 = prune_taxa((tax_table(ps)[, "Order"] %in% top10order), ps)
Ordertop10 <- prune_samples(sample_sums(Ordertop10)>1, Ordertop10)
Ordertop10_rel <- transform_sample_counts(Ordertop10, function(x) x/sum(x))
#plot_bar_2 to remove bar border
bar_top10_order <- plot_bar_2(Ordertop10_rel, fill="Order", border_color = NA)
bar_top10_order <- bar_top10_order + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("indianred2", "seagreen", "darkseagreen","paleturquoise1", "burlywood", "darkgoldenrod2", "deepskyblue4", "deepskyblue3", "darkkhaki", "lightpink3")) +
  theme(legend.position="right") #"cadetblue" #"cornflowerblue" #"darkolivegreen3" #"deepskyblue4"
bar_top10_order <- bar_top10_order + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank())
bar_top10_order + theme(panel.background = element_blank())

#top 20 genus
genus.sum = tapply(taxa_sums(ps), tax_table(ps)[, "Genus"], sum, na.rm=TRUE)
genus.sum
top20genus = names(sort(genus.sum, TRUE))[1:20]
top20genus
Genustop20 = prune_taxa((tax_table(ps)[, "Genus"] %in% top20genus), ps)
Genustop20 <- prune_samples(sample_sums(Genustop20)>1, Genustop20)
Genustop20_rel <- transform_sample_counts(Genustop20, function(x) x/sum(x))
print(top20genus)
bar_top20_genus <- plot_bar_2 (Genustop20_rel, fill="Genus", border_color = NA)
bar_top20_genus
bar_top20_genus <- bar_top20_genus + geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values = c("darkblue", "darkgoldenrod2", "darkorchid", "paleturquoise3", "darkolivegreen1","darkorange","royalblue2", "darksalmon", "khaki2", "indianred2", "darkslategrey", "darkkhaki", "seagreen", "burlywood","paleturquoise1", "plum1",
                                         "firebrick4", "darkolivegreen4", "deepskyblue3", "darkgoldenrod3")) +  #dodgerblue3 #broun1 #lightskyblue #cyan1 #deeppink2 #darkolivegreen4 #coral2 #darkgreen
                                           theme(legend.position="bottom") 
bar_top20_genus + theme(panel.background = element_blank())


#
plot_richness(ps, measures=c("Simpson", "Shannon"))
plot_richness(ps, x="Крупный_плод", color="Крупный_плод", measures=c("Chao1", "Shannon")) + geom_point(size=2, alpha=0.7) +
   theme(panel.background = element_blank())
net <- plot_net(ps, "bray", color = "diet_before", laymeth = "reingold.tilford") +
  scale_color_manual(values = sponge_pallet)
net

plot_richness(ps, x="Крупный_плод", measures=c("Chao1", "Shannon"))

#save data for process
ps@tax_table
ps@otu_table
seqtab.nochim.export <- t(ps@otu_table)
write.table(seqtab.nochim.export, sep="\t", file = "C:/Users/rusla/OneDrive/Desktop/Projects/Hospital/part_4/OTU_97/All_freq_Uterus_vagina.tsv")
write.table(ps@tax_table, sep="\t", file = "C:/Users/rusla/OneDrive/Desktop/Projects/Hospital/part_4/OTU_97/All_otu_Uterus_Vagina_taxa.tsv")
write.table(ps@sam_data, sep="\t", file = "C:/Users/rusla/OneDrive/Desktop/Projects/Hospital/part_4/OTU_97/metadata_Uterus_Vagina.tsv")



#PCoA 3D
#meta <- read.csv(sep = "\t", file = "C:/Users/rusla/OneDrive/Desktop/Projects/Hospital/part_3/metadata_defrost.txt", header = TRUE, check.names=FALSE, row.names=1, skipNul = TRUE)
#otu <- read.csv(sep = "\t", file = "C:/Users/rusla/OneDrive/Desktop/Projects/Hospital/part_3/OTU_97/all_OTU_frequency_sort.tsv", header = TRUE, check.names=FALSE, row.names=1, skipNul = TRUE)
#pcoa3d <- PCoA3D(METRIC="weightedunifrac", METADATA=meta, FEATURES=otu, COLOR="Status", SHAPE="GD", TREE=random_tree, AXIS=c(1,2,3), PALETTE=sponge_pallet_2)
#pcoa3d


# Homogeneity of dispersion test
#beta <- betadisper(erie_bray, sampledf$GD)
#permutest(beta)

#lef_out<-run_lefse(ps1, group = "Состояние", subgroup = "Крупный_плод", norm = "CPM", taxa_rank = "none",
#kw_cutoff = 0.05, lda_cutoff = 2)

#lef_out
#plot_ef_bar(lef_out)
