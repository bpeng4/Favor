library(ggplot2)
library(dplyr)
library(ggrepel)
library(broom)

#### Bifidobacterium
Bifi_Genus <- read.csv(file = "/Users/bopeng/Documents/Favor/Cache/Bifidobacterium_Genus.csv")

Bifi_Genus$Bifi_log2 <- log(Bifi_Genus$Bifidobacterium, base = 2)

#### Bacteroides
Bac_Genus <- read.csv(file = "/Users/bopeng/Documents/Favor/Cache/Bacteroides_Genus.csv")

Bac_Genus$Bac_log10 <- log10(Bac_Genus$Bacteroides)

#### Bif_Bac
Bac_Bif_Genus <- merge(Bifi_Genus, Bac_Genus, by=c("Sample_Name", "Seed_Type", "Microbiome"))

#### Remove Samples with "-Inf" Bifidobacterium values
Bac_Bif_Genus <- Bac_Bif_Genus[Bac_Bif_Genus$Bifi_log2 != "-Inf",]

#### Scatter Plot
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(broom)
library(purrr)  # Needed for map_dbl and map2_dbl

# Optional: for reproducible label layout
set.seed(123)

# Step 1: Calculate R² and Spearman ρ per Seed_Type
model_info <- Bac_Bif_Genus %>%
  group_by(Microbiome) %>%
  summarise(
    # Fit linear model
    mod = list(lm(Bifi_log2 ~ Bac_log10, data = pick(everything()))),
    # Spearman correlation
    spearman_rho = cor(Bac_log10, Bifi_log2, method = "spearman", use = "complete.obs"),
    # Mean x for label position
    x = mean(Bac_log10, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    # Extract R²
    r2 = map_dbl(mod, ~ glance(.x)$r.squared),
    # Predict y at mean x
    y = map2_dbl(mod, x, ~ predict(.x, newdata = data.frame(Bac_log10 = .y))),
    # Construct label
    label = paste0(Microbiome, "\nR² = ", round(r2, 2), "\nρ = ", round(spearman_rho, 2))
  )

# Step 2: Plot with labels
ggplot(Bac_Bif_Genus, aes(x = Bac_log10, y = Bifi_log2, color = Microbiome)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_label_repel(
    data = model_info,
    aes(x = x, y = y, label = label, fill = Microbiome),
    color = "black", show.legend = FALSE,
    size = 3, label.size = 0.3,
    max.overlaps = 100
  ) +
  theme_bw() +
  guides(color = 'none', fill = 'none') +
  labs(x = "Bacteroides (log10(CFU))", y = "Bifidobacterium (log2(CFU))")

