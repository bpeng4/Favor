#Load Package
rm(list = ls())
gc()
suppressMessages({
  library("tidyverse")
  library("phyloseq")
  library("rstatix")
  library("vegan")
  library("picante")
  library("ggpubr")
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)  
  library(tibble)
  library(poppr)  # For AMOVA if genetic analysis is needed
  library(viridis)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(agricolae)
  library(microbiomeSeq)
  library(microbiomeutilities)
  library(viridis)
})
no_of_cores = 16
setwd("/work/benson/bpeng4/Favor/Favor_16S")
load("intermediate/Favor_Forward_ps.rda")
ls()
ps

###Pre-arrangement on  the phyloseq file
################
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))
ps.clean

## Standardize Data by Rarefication
##
#Check Total Reads to find the Thresholds for Rarefication
OTU<-otu_table(ps.clean) |>as.data.frame()
OTU$Total<- rowSums(OTU)
#Rarefication
ps.rare<-rarefy_even_depth(ps.clean, sample.size = 13634,
                           rngseed = 111, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

ps.rare <-  subset_samples(ps.rare, Sample != "FBB0")

ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])

ps_S5S9 <- subset_samples(ps.rare, Microbiome %in% c("S5","S9"))

ps_S6S9S12 <- subset_samples(ps.rare, Microbiome %in% c("S1","S5","S6", "S9"))

ps_S6S9S12Rep1 <- subset_samples(ps_S6S9S12, Rep == "1")
ps_S6S9S12Rep2 <- subset_samples(ps_S6S9S12, Rep == "2")

## LEfSe analysis and plot

lefse_result <- run_lefse(
  ps = ps_S5S9,
  group = "Sample",            
  norm = "CPM",
  lda_cutoff = 2.0,
  kw_cutoff = 0.05,
  wilcoxon_cutoff = 0.05
)




# Cladogram plot
plot_cladogram(lefse_result, color = c('#6A3D9A', '#E41A1C'))



# LDA score barplot: use ggplot2 + marker_table
lda_tbl <- marker_table(lefse_result)

# LEfSe LDA Barplot
library(dplyr)
library(stringr)
library(ggplot2)

lda_tbl_df <- as_tibble(lda_tbl) %>%
  mutate(feature = str_trunc(feature, 200)) %>%
  mutate(genus = word(feature, 1)) %>%
  arrange(desc(ef_lda)) %>%
  mutate(genus = factor(genus, levels = rev(unique(genus))))

lda_tbl_df <- lda_tbl_df %>%
  mutate(last_level = str_extract(feature, "[^|]+$"))

# Define your colors vector (adjust length as needed)
my_colors <- c('#6A3D9A', '#FFFF99','#F781BF', '#CCEBC5',  '#FF7F00', '#E41A1C')

ggplot(lda_tbl_df, aes(x = last_level, y = ef_lda, fill = enrich_group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = my_colors) +
  labs(x = "Taxa", y = "LDA Score", title = "LEfSe LDA Barplot") +
  theme_minimal(base_size = 14)







