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
  library(poppr)  # For AMOVA if genetic analysis is needed
  library(viridis)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
  library(agricolae)
  library(microbiomeSeq)
  library(microbiomeutilities)
  library(viridis)
  library(RColorBrewer)
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

## Standardize Data by Calculating Relative Abundance
##
#Filter Out Taxa Exist in At Least 25% samples
ps.clean.p0 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 251}, prune=TRUE)
ps.clean.p0

#Change to relative abundance
ps.clean.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))
ps.sample.re<-subset_samples(ps.clean.re, !ps.clean.re@sam_data$Sample == "FBB0")
######################

###Plot Stack Bar chart for Baseline relative abundance
####################
ps.re.FBB<-subset_samples(ps.clean.re, ps.clean.re@sam_data$Sample == "FBB0")

#Set up mycolor
mycolor<-c(
  '#F781BF', '#CCEBC5', '#FFFF99', '#FF7F00', '#8DD3C7', '#7FC97F', '#999999',  
  '#CAB2D6', '#FDBF6F', '#D95F02', '#666666', '#A6761D', '#B2DF8A', '#377EB8', 
  '#6A3D9A', '#FB9A99', '#E41A1C', '#FFFFB3', '#FDB462', '#F0027F', '#D9D9D9', '#4DAF4A', 
  '#FFFF33', '#B15928', '#FDC086', '#80B1D3', '#A6CEE3', '#BC80BD', '#1B9E77', '#E31A1C',  
  '#FCCDE5', '#386CB0', '#984EA3', '#FFED6F', '#66A61E', '#E6AB02',  
  '#BF5B17', '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C','#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4','#a48cf4', 
  '#f45cf2', '#f66bad', '#4A192C', '#8D948D', '#FE0000', '#F4A900', '#49678D',
  '#922B3E', '#1E213D', '#E6D690', '#B44C43', '#47402E', '#2C5545', '#F5D033')

sample_order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
                  "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16",
                  "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24")

ps.re.FBB@sam_data$Microbiome <- factor(ps.re.FBB@sam_data$Microbiome, levels = sample_order)


#Plot for Microbiome Baseline Profile
phyloseq::plot_bar(ps.re.FBB, fill = "Family") + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Microbiome, scales = "free_x",ncol = 5) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, size = 8),
        legend.key.size = unit(17, "pt"))+
  scale_fill_manual(values=mycolor)+
  guides(fill = guide_legend(ncol = 1))

########################
####Alpha Diversity Analysis and Plotting
## alpha diversity should be calculated before filtering on abundance and prevalence
#Remove NCs
tree = phyloseq::phy_tree(ps)
samp = data.frame(phyloseq::otu_table(ps))

##Alpha Diversity Among the Samples
adiv <- data.frame(
  phyloseq::estimate_richness(ps, measures = c( "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps)), Sample_Name, Microbiome, Rep,Sample)) %>%
  dplyr::select(-se.chao1)

#Glance the median values
adiv %>%
  group_by( Sample ) %>%
  dplyr::summarise(median_Chao1 = median(Chao1),
                   median_Shannon = median(Shannon),
                   median_Simpson = median(Simpson),
                   median_InvSimpson = median(InvSimpson),
                   median_Fisher = median(Fisher),
                   median_PD = median(PD))

## Statistical Tests
# ANOVA + Duncan
model <- aov(Simpson ~ Sample + Microbiome + Rep, data = adiv)
summary(model)
duncan_result <- duncan.test(model, "Sample", console = TRUE)

write.csv(duncan_result[4], file = "/work/benson/bpeng4/Favor/Plots/AlphaDiversity_Duncan_Mean_SE.csv" )
write.csv(duncan_result[6], file = "/work/benson/bpeng4/Favor/Plots/AlphaDiversity_Post_Hoc_Comparison.csv" )


## Plotting
#Order Sample Name
sample_order <- c("FBB0", "A554", "PHG47", "PHG39", "B57", "HP301",
                  "LH60", "764", "PHG35", "NC262")
adiv0 <- adiv %>%
  filter(Sample %in% sample_order) %>%
  mutate(Sample = factor(Sample, levels = sample_order)) %>%
  arrange(Sample)

adiv0 <- adiv0[adiv0$Sample != "FBB0",]
  
adiv0 %>%
  pivot_longer(cols = c("Shannon", "Chao1", "Simpson", "PD"), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
  ggplot(aes(x = Sample, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Sample), height = 0, width = .2) +
  labs(
    title = "Alpha Diversity Measures",
    x = "", 
    y = "Alpha Diversity Measures"
  )+
  facet_wrap(~ metric, scales = "free", ncol = 1) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14),
        title = element_text(size = 18)) +
  stat_anova_test()



AdivPlot<-function(adiv, title_text){
  adiv %>%
    pivot_longer(cols = c("Shannon", "Chao1", "Simpson", "PD"), 
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
    ggplot(aes(x = Sample, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = Sample), height = 0, width = .2) +
    labs(
      title = title_text,
      x = "", 
      y = "Alpha Diversity Measures"
    )+
    labs(x = "", y = "Alpha Diversity Measures") +
    facet_wrap(~ metric, scales = "free", ncol = 1) +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 14),
          title = element_text(size = 18)) +
    stat_anova_test()
}

# Loop with dynamic titles
for (i in c("S1", "S2", "S3", "S4","S5", "S6", "S7", 
            "S9", "S12", "S14", "S15", "S16",
            "S17", "S18", "S19", "S20", "S21", "S22", "S24")) {
  
  subset_data <- adiv0[adiv0$Microbiome == i, ]
  plot_title <- paste(i, "Alpha Diversity by Seed Type")
  print(AdivPlot(subset_data, plot_title))
}
 


####Beta Diversity Analysis and Plotting
########################
ps.beta= ps.clean.re

###Principle Component Analysis
# Bray-Curtis
ordBC <- ordinate(ps.beta, "PCoA", "bray")
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")
smpID <- ps.beta@sam_data$Sample_Name

#Calculate the variance explained by the first two components
variance_BC <- ordBC$values$Relative_eig[1:2] * 100  # convert to percentage
names(variance_BC) <- c("PC1", "PC2")
variance_BC
variance_JC <- ordJC$values$Relative_eig[1:2] * 100
names(variance_JC) <- c("PC1", "PC2")
variance_JC
variance_UF <- ordUF$values$Relative_eig[1:2] * 100
names(variance_UF) <- c("PC1", "PC2")
variance_UF
variance_wUF <- ordwUF$values$Relative_eig[1:2] * 100
names(variance_wUF) <- c("PC1", "PC2")
variance_wUF

# add sample_data info
BC <- merge(data.frame(ordBC$vectors[,1:2], sample = smpID, method = 'BC'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(BC) <- BC$Row.names
BC$Row.names <- NULL  # Remove the extra column created by merge

Jaccard <- merge(data.frame(ordJC$vectors[,1:2], sample = smpID,method = 'Jaccard'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(Jaccard) <- Jaccard$Row.names
Jaccard$Row.names <- NULL  # Remove the extra column created by merge

unifrac <- merge(data.frame(ordUF$vectors[,1:2], sample = smpID,method = 'unifrac'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(unifrac) <- unifrac$Row.names
unifrac$Row.names <- NULL  # Remove the extra column created by merge

wunifrac <- merge(data.frame(ordwUF$vectors[,1:2], sample = smpID,method = 'wunifrac'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(wunifrac) <- wunifrac$Row.names
wunifrac$Row.names <- NULL  # Remove the extra column created by merge

# Keep first 2 vectors (latent variables, PCs) of each distance matrix
df <- rbind(BC,Jaccard,unifrac,wunifrac)

#Form a column named MicroSamp
df$MicroSamp <- paste(df$Microbiome, df$Sample, sep = "_")
# Calculate the mean for each level of MicroSamp
summary_df <- df %>%
  group_by(method,Microbiome,Sample, MicroSamp) %>%
  summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE)
  )
#Order Category. levels
summary_df$Sample <- factor(summary_df$Sample, 
                               levels = c("FBB0", "A554", "PHG47", "PHG39", "B57", "HP301",
                                           "LH60", "764", "PHG35", "NC262"))
#Plot without arrows
# Define color palette
yourcolor <- c(
  '#F781BF', '#FF7F00', '#8DD3C7', '#7FC97F', '#A6761D', '#6A3D9A', '#FB9A99',
  '#E41A1C', '#FDB462', '#F0027F', '#D9D9D9', '#FFFF33', 
  '#FDC086', '#80B1D3', '#BC80BD', '#1B9E77', '#E31A1C', '#FCCDE5',
  '#FFED6F', '#66A61E', '#E6AB02', '#BF5B17',
  '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C', '#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4',
  '#a48cf4', '#f45cf2', '#f66bad'
)

#Set up FBB0 as the starting points
summary_df0<-summary_df
fbb0_df <- summary_df0 %>%
  filter(Sample == "FBB0") %>%
  select(method, Microbiome, Axis.1, Axis.2) %>%
  mutate(Axis.1_start = Axis.1, Axis.2_start = Axis.2) %>%
  select(-Axis.1, -Axis.2)  # Remove old columns to avoid duplicates

summary_df0 <- summary_df0 %>%
  left_join(fbb0_df, by = c("Microbiome", "method"))  # Match on method & Microbiome

# Order the levels
sample_order <- c("S1", "S2", "S3", "S4","S5", "S6", "S7", 
                  "S9", "S12", "S14", "S15", "S16",
                  "S17", "S18", "S19", "S20", "S21", "S22", "S24")

summary_df0$Microbiome <- factor(summary_df0$Microbiome, levels = sample_order)


#Plot with arrows
# Add the explanation rate for each principle component
variance_labels <- c(
  BC = "Bray-Curtis     PC1 (30.77%) vs PC2 (12.40%)",
  Jaccard = "Jaccard     PC1 (16.24%) vs PC2 (6.20%)",
  unifrac = "Unifrac     PC1 (34.14%) vs PC2 (14.38%)",
  wunifrac = "Weighted Unifrac     PC1 (46.74%) vs PC2 (16.02%)"
)

ggplot(summary_df0, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") + 
  facet_wrap(~method, scales = 'free', 
             labeller = as_labeller(variance_labels)) +  # <-- Add custom labels
  scale_color_manual(values = yourcolor) +
  scale_shape_manual(values = c(19,2:10)) +
  labs(shape = "Seed Type", color = "Microbiome",
       x = "PC1", y = "PC2") +  # You can still label axes generically
  theme_minimal()
#############

#Beta Diversity Plot with arrows for each microbiome
####################
# Define plotting function
BetaPlot <- function(sub_df, micro_name, color_palette) {
  ggplot(sub_df, aes(Axis.1, Axis.2, color = Microbiome, shape = Sample.x)) + 
    geom_point(size = 3) + 
    geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                     xend = Axis.1, yend = Axis.2),
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = "black") + 
    facet_wrap(~method, scales = 'free', 
               labeller = as_labeller(variance_labels)) +  # <-- Add custom labels
    scale_color_manual(values = yourcolor) +
    scale_shape_manual(values = c(19,2:10)) +
    labs(shape = "Seed Type", color = "Microbiome",
         x = "PC1", y = "PC2") +  # You can still label axes generically
    theme_minimal()+
    scale_color_manual(values = yourcolor_named) +
    ggtitle(paste(micro_name, "Beta Diversity Shift from the Baseline"))
}

# Loop over each microbiome and generate plot
microbiomes <- c("S1", "S2", "S3", "S4","S5", "S6", "S7", 
                 "S9", "S12", "S14", "S15", "S16",
                 "S17", "S18", "S19", "S20", "S21", "S22", "S24")

yourcolor_named <- setNames(yourcolor[1:length(microbiomes)], microbiomes)


for (i in seq_along(microbiomes)) {
  micro_name <- microbiomes[i]
  sub_df <- summary_df0[summary_df0$Microbiome == micro_name, ]
  
  if (nrow(sub_df) > 0) {
    print(BetaPlot(sub_df, micro_name, yourcolor))
  } else {
    message(paste("No data for", micro_name, "- skipping plot."))
  }
}



####Heatmap
#Heatmap for Category
####################
#glommate to Genus Level
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])

#Plotting
heat.sample <- plot_taxa_heatmap(ps.genus,
                                 subset.top = 50,
                                 VariableA = "Sample",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10")

#############Calculate the average read counts for each sample
# Extract OTU (abundance) table and metadata from phyloseq
otu_table_df <- as.data.frame(otu_table(ps.genus))
meta_df <- as.data.frame(sample_data(ps.genus))

# Add sample information to the OTU table
otu_table_df$Sample <- meta_df$Sample[match(rownames(otu_table_df), rownames(meta_df))]

# Aggregate by sample (mean abundance per sample)
otu_table_avg <- otu_table_df %>%
  group_by(Sample) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Convert back to matrix and ensure row names are correct
otu_table_avg <- column_to_rownames(otu_table_avg, var = "Sample") %>%
  as.matrix()

# Subtract values in the "FBB0" row from every row
otu_table_avg <- sweep(otu_table_avg, 2, otu_table_avg["FBB0", ], FUN = "-")
otu_table_avg <- otu_table_avg[rownames(otu_table_avg) != "FBB0",]

# Convert back to a phyloseq OTU table
otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)

# Ensure tax_table matches the taxa in otu_table_avg_phylo
common_taxa <- intersect(rownames(tax_table(ps.genus)), colnames(otu_table_avg_phylo))
tax_table_filtered <- tax_table(ps.genus)[common_taxa, ]

# Create a new sample metadata table
meta_avg_phylo <- sample_data(data.frame(Sample = rownames(otu_table_avg_phylo),
                                         row.names = rownames(otu_table_avg_phylo)))

# Set factor levels for Sample in the desired order
meta_avg_phylo$Sample <- factor(meta_avg_phylo$Sample, 
                                  levels = c("A554", "PHG47", "PHG39", "B57", "HP301",
                                             "LH60", "764", "PHG35", "NC262"))  # Replace with actual category names

# Create a new phyloseq object with matching OTU and tax tables
ps.genus.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)

# Generate heatmap of top 50 taxa, now grouped by sample
heat.sample <- plot_taxa_heatmap(ps.genus.avg,
                                   subset.top = 50,
                                   VariableA = "Sample",
                                   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                   transformation = "log10")

# Plot the heatmap
heat.sample

####################

#Heatmap for Sample
##############
heatmapsample <- function(ps.genus) {

  # Extract OTU (abundance) table and metadata from phyloseq
  otu_table_df <- as.data.frame(otu_table(ps.genus))
  meta_df <- as.data.frame(sample_data(ps.genus))
  
  # Add sample information to the OTU table
  otu_table_df$Sample <- meta_df$Sample[match(rownames(otu_table_df), rownames(meta_df))]
  
  # Aggregate by sample (mean abundance per sample)
  otu_table_avg <- otu_table_df %>%
    group_by(Sample) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # Convert back to matrix and ensure row names are correct
  otu_table_avg <- column_to_rownames(otu_table_avg, var = "Sample") %>%
    as.matrix()
  
  # Subtract values in the "FBB0" row from every row
  otu_table_avg <- sweep(otu_table_avg, 2, otu_table_avg["FBB0", ], FUN = "-")
  otu_table_avg <- otu_table_avg[rownames(otu_table_avg) != "FBB0",]
  
  # Convert back to a phyloseq OTU table
  otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)
  
  # Ensure tax_table matches the taxa in otu_table_avg_phylo
  common_taxa <- intersect(rownames(tax_table(ps.genus)), colnames(otu_table_avg_phylo))
  tax_table_filtered <- tax_table(ps.genus)[common_taxa, ]
  
  # Create a new sample metadata table
  meta_avg_phylo <- sample_data(data.frame(Sample = rownames(otu_table_avg_phylo),
                                           row.names = rownames(otu_table_avg_phylo)))
  
  # Set factor levels for Sample in the desired order
  meta_avg_phylo$Sample <- factor(meta_avg_phylo$Sample, 
                                  levels = c("A554", "PHG47", "PHG39", "B57", "HP301",
                                             "LH60", "764", "PHG35", "NC262"))  # Replace with actual category names
  
  # Create a new phyloseq object with matching OTU and tax tables
  ps.genus.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
  
  # Step 2: Apply log2(x + 1)
  ps.genus.log2 <- transform_sample_counts(ps.genus.avg, function(x) log2(x + 1))
  
  # Step 3: Convert log2(x + 1) to 10^(log2(x + 1))
  ps.genus.log2.10 <- transform_sample_counts(ps.genus.log2, function(x) 10^x)
  
  # Plot heatmap of top 50 taxa
  heat.sample <- plot_taxa_heatmap(
    ps.genus.log2.10,
    subset.top = 50,
    VariableA = c("Sample"),
    heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    transformation = "log10"
  )
  return(heat.sample)
}

##############


###Permanova Calculation and Visualization
##############
#####Write the function 
ps.filter<-function(ps.genus){
  #Extract OTU/ASV table and metadata:
  otu_table <- as.data.frame(otu_table(ps.genus))
  metadata <- data.frame(sample_data(ps.genus))  # Force conversion
  taxa_table<- as.data.frame(tax_table(ps.genus))

  #Use the adonis function (PERMANOVA) from vegan to analyze treatment(Sample) effects
  metadata$Sample <- as.factor(metadata$Sample)  # Ensure it's a factor
  adonis_results <- adonis2(otu_table ~ Sample, data = metadata, permutations = 999, method = "bray")
  
  #Extract p-values
  p_value <- adonis_results$`Pr(>F)`[1]  # Extract p-value for treatment
  
  #test each taxon individually:
  p_values <- apply(otu_table, 2, function(x) {
    df <- data.frame(x = x, Sample = metadata$Sample)
    fit <- aov(x ~ Sample, data = df)
    summary(fit)[[1]][["Pr(>F)"]][1]  # Extract p-value
  })
  
  #Adjust for multiple testing (FDR correction):
  p_values_adj <- p.adjust(p_values, method = "fdr")
  
  #Mark taxa as significant if p < 0.05
  signif_taxa <- names(p_values_adj[p_values_adj < 0.05])

  #Filter for only the significant taxa
  ps.genus.filtered <- prune_taxa(taxa_names(ps.genus) %in% signif_taxa, ps.genus)
  
  return(ps.genus.filtered)
}
###############

###Heatmap for significant taxa calculated by Permanova for each Subject 
#glommate to Genus Level
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])
ps.genus.filtered <- ps.filter(ps.genus)
heatmapsample(ps.genus.filtered)

#Subset for each Subject
ps.genus.S1<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S1")
ps.genus.filtered.S1 <- ps.filter(ps.genus.S1)
heatmapsample(ps.genus.filtered.S1)

ps.genus.S2<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S2")
ps.genus.filtered.S2 <- ps.filter(ps.genus.S2)
#heatmapsample(ps.genus.filtered.S2)

ps.genus.S3<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S3")
ps.genus.filtered.S3 <- ps.filter(ps.genus.S3)
heatmapsample(ps.genus.filtered.S3)

ps.genus.S4<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S4")
ps.genus.filtered.S4 <- ps.filter(ps.genus.S4)
heatmapsample(ps.genus.filtered.S4)

ps.genus.S7<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S7")
ps.genus.filtered.S7 <- ps.filter(ps.genus.S7)
heatmapsample(ps.genus.filtered.S7)

ps.genus.S9<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S9")
ps.genus.filtered.S9 <- ps.filter(ps.genus.S9)
heatmapsample(ps.genus.filtered.S9)

ps.genus.S12<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S12")
ps.genus.filtered.S12 <- ps.filter(ps.genus.S12)
heatmapsample(ps.genus.filtered.S12)

ps.genus.S14<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S14")
ps.genus.filtered.S14 <- ps.filter(ps.genus.S14)
heatmapsample(ps.genus.filtered.S14)

ps.genus.S15<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S15")
ps.genus.filtered.S15 <- ps.filter(ps.genus.S15)
heatmapsample(ps.genus.filtered.S15)

ps.genus.S16<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S16")
ps.genus.filtered.S16 <- ps.filter(ps.genus.S16)
heatmapsample(ps.genus.filtered.S16)

ps.genus.S17<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S17")
ps.genus.filtered.S17 <- ps.filter(ps.genus.S17)
heatmapsample(ps.genus.filtered.S17)

ps.genus.S19<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S19")
ps.genus.filtered.S19 <- ps.filter(ps.genus.S19)
heatmapsample(ps.genus.filtered.S19)

ps.genus.S20<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S20")
ps.genus.filtered.S20 <- ps.filter(ps.genus.S20)
heatmapsample(ps.genus.filtered.S20)

ps.genus.S21<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S21")
ps.genus.filtered.S21 <- ps.filter(ps.genus.S21)
heatmapsample(ps.genus.filtered.S21)

ps.genus.S24<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S2")
ps.genus.filtered.S24 <- ps.filter(ps.genus.S24)
#heatmapsample(ps.genus.filtered.S24)

#ps.rare<-subset_samples(ps.rare, !ps.rare@sam_data$Category == "FBB16")

ps.family = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[5])
ps.family.filtered <- ps.filter(ps.family)
heatmapsample(ps.family.filtered)

ps.order = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[4])
ps.order.filtered <- ps.filter(ps.order)
heatmapsample(ps.order.filtered)

ps.class = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[3])
ps.class.filtered <- ps.filter(ps.class)
heatmapsample(ps.class.filtered)


ps.genus.sample <- subset_samples(ps.genus, !ps.genus@sam_data$Category %in% c("FBB0","FS"))
#Modify mannually by add "1" to each value in ps.genus otu table
otu_table(ps.genus.sample) <- otu_table(ps.genus.sample) + 1
ps.genus.sample.filtered <- ps.filter(ps.genus.sample)
heatmapsample(ps.genus.sample.filtered)



##############
####Prebiotic Potential Index Calculation
##Define Beneficial Bacteria and Harmful Bacteria
Beneficial<-c("Akkermansia", "Anaerostipes", "Barnesiella", "Bifidobacterium",
              "Blautia", "Butyricicoccus", "Catenibacterium", "Coprococcus",
              "Eubacterium", "Faecalibacterium", "Fusicatenibacter", "Lactobacillus",
              "Megasphaera", "Oscillibacter", "Parabacteroides", "Prevotella",
              "Roseburia", "Ruminococcus")
Harmful <- c("Bilophila", "Desulfovibrio", "Escherichia-Shigella","Fusobacterium", 
             "Haemophilus", "Klebsiella", "Paraprevotella", "Parasutterella", 
             "Streptococcus", "Sutterella", "Veillonella")

####Glomerate the phyloseq files to Beneficial and Harmful Categories
########
#Write a function to classified taxa to Beneficial and Harmful Categories
TaxGlom<-function(ps.rare){
  ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])
  
  genus_level <- rank_names(ps.rare)[6]  # Usually "Genus"
  
  #Modify the tax_table to include a new "Group" column
  tax_df <- as.data.frame(tax_table(ps.genus))
  tax_df$Group <- ifelse(tax_df[[genus_level]] %in% Beneficial, "Beneficial",
                         ifelse(tax_df[[genus_level]] %in% Harmful, "Harmful", "Other"))
  
  #Reassign the modified taxonomy table back to the phyloseq object
  tax_table(ps.genus) <- tax_table(as.matrix(tax_df))
  
  # Rename the new column as a taxonomic rank (so tax_glom can use it)
  colnames(tax_table(ps.genus))[ncol(tax_table(ps.genus))] <- "Group"
  
  #Glom based on your new "Group" column
  ps.grouped <- tax_glom(ps.genus, taxrank = "Group")
  
  #(optional): Remove "Other" if you only want Beneficial + Harmful
  ps.grouped <- subset_taxa(ps.grouped, Group %in% c("Beneficial", "Harmful"))
  
  #Extract OTU (ASV) table and tax_table
  otu_df <- as.data.frame(t(otu_table(ps.grouped)))
  tax_df <- as.data.frame(tax_table(ps.grouped))
  
  #Add Group column to OTU table using tax_table
  otu_df$Group <- tax_df$Group
  
  #Sum abundance across samples by Group
  grouped_abund <- otu_df %>%
    group_by(Group) %>%
    summarise(across(everything(), sum))
  
  #Convert back to phyloseq components (optional)
  otu_grouped <- as.matrix(grouped_abund[, -1])
  rownames(otu_grouped) <- grouped_abund$Group
  
  otu_ps <- otu_table(otu_grouped, taxa_are_rows = TRUE)
  otu_ps <- t(otu_ps)
  tax_grouped <- tax_table(matrix(data = grouped_abund$Group, ncol = 1,
                                  dimnames = list(grouped_abund$Group, "Group")))
  
  #Create new phyloseq object
  ps.grouped_final <- phyloseq(otu_ps, tax_grouped, sample_data(ps.grouped))
  
  return(ps.grouped_final)
}

#########
ps.grouped_final<-TaxGlom(ps.rare)

ps.rareS1<-subset_samples(ps.rare, Microbiome == "S1")
ps.rareS2<-subset_samples(ps.rare, Microbiome == "S2")
ps.rareS3<-subset_samples(ps.rare, Microbiome == "S3")
ps.rareS4<-subset_samples(ps.rare, Microbiome == "S4")
ps.rareS5<-subset_samples(ps.rare, Microbiome == "S5")
ps.rareS6<-subset_samples(ps.rare, Microbiome == "S6")
ps.rareS7<-subset_samples(ps.rare, Microbiome == "S7")
ps.rareS8<-subset_samples(ps.rare, Microbiome == "S8")
ps.rareS9<-subset_samples(ps.rare, Microbiome == "S9")
ps.rareS10<-subset_samples(ps.rare, Microbiome == "S10")
ps.rareS11<-subset_samples(ps.rare, Microbiome == "S11")
ps.rareS12<-subset_samples(ps.rare, Microbiome == "S12")
ps.rareS13<-subset_samples(ps.rare, Microbiome == "S13")
ps.rareS14<-subset_samples(ps.rare, Microbiome == "S14")
ps.rareS15<-subset_samples(ps.rare, Microbiome == "S15")
ps.rareS16<-subset_samples(ps.rare, Microbiome == "S16")

ps.S1.grouped_final<-TaxGlom(ps.rareS1)
ps.S2.grouped_final<-TaxGlom(ps.rareS2)
ps.S3.grouped_final<-TaxGlom(ps.rareS3)
ps.S4.grouped_final<-TaxGlom(ps.rareS4)
ps.S5.grouped_final<-TaxGlom(ps.rareS5)
ps.S6.grouped_final<-TaxGlom(ps.rareS6)
ps.S7.grouped_final<-TaxGlom(ps.rareS7)
ps.S8.grouped_final<-TaxGlom(ps.rareS8)
ps.S9.grouped_final<-TaxGlom(ps.rareS9)
ps.S10.grouped_final<-TaxGlom(ps.rareS10)
ps.S11.grouped_final<-TaxGlom(ps.rareS11)
ps.S12.grouped_final<-TaxGlom(ps.rareS12)
ps.S13.grouped_final<-TaxGlom(ps.rareS13)
ps.S14.grouped_final<-TaxGlom(ps.rareS14)
ps.S15.grouped_final<-TaxGlom(ps.rareS15)
ps.S16.grouped_final<-TaxGlom(ps.rareS16)


#Calculate the average read counts for each category
#############
# Extract OTU (abundance) table and metadata from phyloseq
otu_table_df <- as.data.frame(otu_table(ps.grouped_final))
meta_df <- as.data.frame(sample_data(ps.grouped_final))

# Add category information to the OTU table
otu_table_df$Category <- meta_df$Category[match(rownames(otu_table_df), rownames(meta_df))]

# Aggregate by category (mean abundance per category)
otu_table_avg <- otu_table_df %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Convert back to matrix and ensure row names are correct
otu_table_avg <- column_to_rownames(otu_table_avg, var = "Category") %>%
  as.matrix()

# Convert back to a phyloseq OTU table
otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)

# Ensure tax_table matches the taxa in otu_table_avg_phylo
common_taxa <- intersect(rownames(tax_table(ps.grouped_final)), colnames(otu_table_avg_phylo))
tax_table_filtered <- tax_table(ps.grouped_final)[common_taxa, ]

# Create a new sample metadata table
meta_avg_phylo <- sample_data(data.frame(Category = rownames(otu_table_avg_phylo),
                                         row.names = rownames(otu_table_avg_phylo)))

# Set factor levels for Category in the desired order
meta_avg_phylo$Category <- factor(meta_avg_phylo$Category, 
                                  levels = c("FBB0", "FBB16",                       
                                              "FS", "Postbiotic","Prebiotic control",           
                                              "Prebiotic fiber", "Synergistic effect (Fiber)",  
                                              "Vitamin","Vitamin (Synergistic effect)",
                                              "Vitamin (Synergistic)","Vitamin + fiber"
                                              ))  # Replace with actual category names

# Create a new phyloseq object with matching OTU and tax tables
ps.grouped.final.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
############
ps.grouped.final.avg

#####Calculate  PPIa for each Category level
OTU_Average<- ps.grouped.final.avg@otu_table |> as.data.frame()
OTU_Average$Total<-rowSums(OTU_Average)
OTU_Average$PPIa<-(OTU_Average$Beneficial- OTU_Average$Harmful)/OTU_Average$Total 

write.csv(OTU_Average, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_Category.csv")



#Write a funciton to calculate the average read counts for each sample
#############

SampleRead <- function(ps.grouped_final){
  # Extract OTU (abundance) table and metadata from phyloseq
  otu_table_df <- as.data.frame(otu_table(ps.grouped_final))
  meta_df <- as.data.frame(sample_data(ps.grouped_final))
  
  # Add sample information to the OTU table
  otu_table_df$Sample <- meta_df$Sample[match(rownames(otu_table_df), rownames(meta_df))]
  
  # Aggregate by sample (mean abundance per sample)
  otu_table_avg <- otu_table_df %>%
    group_by(Sample) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # Convert back to matrix and ensure row names are correct
  otu_table_avg <- column_to_rownames(otu_table_avg, var = "Sample") %>%
    as.matrix()
  
  # Convert back to a phyloseq OTU table
  otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)
  
  # Ensure tax_table matches the taxa in otu_table_avg_phylo
  common_taxa <- intersect(rownames(tax_table(ps.grouped_final)), colnames(otu_table_avg_phylo))
  tax_table_filtered <- tax_table(ps.grouped_final)[common_taxa, ]
  
  # Create a new sample metadata table
  meta_avg_phylo <- sample_data(data.frame(Sample = rownames(otu_table_avg_phylo),
                                           row.names = rownames(otu_table_avg_phylo)))
  
  # Set factor levels for Sample in the desired order
  meta_avg_phylo$Sample <- factor(meta_avg_phylo$Sample, 
                                  levels = c("beetroot powder", "BetaVia Pure WD",                     
                                             "carrot powder", "Cellulose",                           
                                             "dried apple pomace", "dried spinich powder",                
                                             "ginger powder", "Ground Beet Pulp","guar gum",                  
                                             "HiSmooth flax seed fiber", "Jerusalum articohoke powder (inulin)",
                                             "kappa carageenan ", "lemon powder",                       
                                             "locust bean gum", "Micronized cassava fiber",            
                                             "Mushroom chitosan", "Oat Bagasse hydrolysate-(B)",         
                                             "Oat Bagasse hydrolysate-(PBV)", "Oat Bagasse-Control",                 
                                             "Olive Flour", "pomegranate powder",                  
                                             "Psyllium Husk","xanthan gum","FBB16","FBB0"))  # Replace with actual category names
  
  # Create a new phyloseq object with matching OTU and tax tables
  ps.grouped.final.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
  
  OTU_Average<- ps.grouped.final.avg@otu_table |> as.data.frame()
  OTU_Average$Total<-rowSums(OTU_Average)
  OTU_Average$PPIa<-(OTU_Average$Beneficial- OTU_Average$Harmful)/OTU_Average$Total 
  
  return(OTU_Average)
}

############
OTU_Average <- SampleRead(ps.grouped_final)
OTU_Average_S1 <- SampleRead(ps.S1.grouped_final)
OTU_Average_S2 <- SampleRead(ps.S2.grouped_final)
OTU_Average_S3 <- SampleRead(ps.S3.grouped_final)
OTU_Average_S4 <- SampleRead(ps.S4.grouped_final)
OTU_Average_S5 <- SampleRead(ps.S5.grouped_final)
OTU_Average_S6 <- SampleRead(ps.S6.grouped_final)
OTU_Average_S7 <- SampleRead(ps.S7.grouped_final)
OTU_Average_S8 <- SampleRead(ps.S8.grouped_final)
OTU_Average_S9 <- SampleRead(ps.S9.grouped_final)
OTU_Average_S10 <- SampleRead(ps.S10.grouped_final)
OTU_Average_S11 <- SampleRead(ps.S11.grouped_final)
OTU_Average_S12 <- SampleRead(ps.S12.grouped_final)
OTU_Average_S13 <- SampleRead(ps.S13.grouped_final)
OTU_Average_S14 <- SampleRead(ps.S14.grouped_final)
OTU_Average_S15 <- SampleRead(ps.S15.grouped_final)
OTU_Average_S16 <- SampleRead(ps.S16.grouped_final)


#####Calculate  PPIa for each Category level


write.csv(OTU_Average, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_Sample.csv")
write.csv(OTU_Average_S1, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S1.csv")
write.csv(OTU_Average_S2, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S2.csv")
write.csv(OTU_Average_S3, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S3.csv")
write.csv(OTU_Average_S4, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S4.csv")
write.csv(OTU_Average_S5, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S5.csv")
write.csv(OTU_Average_S6, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S6.csv")
write.csv(OTU_Average_S7, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S7.csv")
write.csv(OTU_Average_S8, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S8.csv")
write.csv(OTU_Average_S9, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S9.csv")
write.csv(OTU_Average_S10, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S10.csv")
write.csv(OTU_Average_S11, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S11.csv")
write.csv(OTU_Average_S12, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S12.csv")
write.csv(OTU_Average_S13, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S13.csv")
write.csv(OTU_Average_S14, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S14.csv")
write.csv(OTU_Average_S15, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S15.csv")
write.csv(OTU_Average_S16, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S16.csv")



#########GMHI Calculation
ps.species = phyloseq::tax_glom(ps.clean.re, taxrank = rank_names(ps.clean.re)[7])
#Output otu tables
OutOTU<-function(ps.run, filename = "otu_otuput"){
  OTU<-otu_table(ps.run) |>as.data.frame()
  Taxa<-tax_table(ps.run) |>as.data.frame()
  Taxa$Name<- apply(Taxa[,6:7], 1, function(x) paste(x, collapse = "_"))
  colnames(OTU)<-Taxa$Name
  OTU<-t(OTU) |> as.data.frame()
  return(OTU)
}

OTU<-OutOTU(ps.species)

write.csv(OTU, file = "/work/benson/bpeng4/Givaudan/species_relative_abundances.csv")
