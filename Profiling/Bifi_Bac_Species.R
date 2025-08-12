library(ggplot2)
library(dplyr)
library(ggrepel)
library(broom)

#### Bifidobacterium
Bifi_Genus <- read.csv(file = "/Users/bopeng/Documents/Favor/Cache/Bif_Species.csv")

#### Bacteroides
Bac_Genus <- read.csv(file = "/Users/bopeng/Documents/Favor/Cache/Bac_Species.csv")

#### Bif_Bac
Bac_Bif_Genus <- merge(Bifi_Genus, Bac_Genus, by=c("Sample_Name", "Seed_Type", "Microbiome"))

#### Scatter Plot
# Load required libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(broom)
library(purrr)  # Needed for map_dbl and map2_dbl

# Step 1: Calculate R² and Spearman ρ per Seed_Type
model_info <- Bac_Bif_Genus %>%
  group_by(Microbiome) %>%
  summarise(
    # Fit linear model
    mod = list(lm(Bif_adolescentis_log2 ~ Bac_ovatus_log2, data = pick(everything()))),
    # Spearman correlation
    spearman_rho = cor(Bac_ovatus_log2, Bif_adolescentis_log2, method = "spearman", use = "complete.obs"),
    # Mean x for label position
    x = mean(Bac_ovatus_log2, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    # Extract R²
    r2 = map_dbl(mod, ~ glance(.x)$r.squared),
    # Predict y at mean x
    y = map2_dbl(mod, x, ~ predict(.x, newdata = data.frame(Bac_ovatus_log2 = .y))),
    # Construct label
    label = paste0(Microbiome, "\nR² = ", round(r2, 2), "\nρ = ", round(spearman_rho, 2))
  )

# Step 2: Plot with labels
ggplot(Bac_Bif_Genus, aes(x = Bac_ovatus_log2, y = Bif_adolescentis_log2, color = Microbiome)) +
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
  labs(x = "Bacteroides ovatus(log2(CFU))", y = "Bifidobacterium adolescentis(log2(CFU))")


### Write a function for correlation and looping
library(rlang)

B_B_Correlation <- function(data, bif_col, bac_col) {
  bif_sym <- sym(bif_col)
  bac_sym <- sym(bac_col)
  
  model_info <- data %>%
    group_by(Microbiome) %>%
    summarise(
      mod = list(lm(!!bif_sym ~ !!bac_sym, data = pick(everything()))),
      spearman_rho = cor(.data[[bac_col]], .data[[bif_col]], method = "spearman", use = "complete.obs"),
      x = mean(.data[[bac_col]], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      r2 = map_dbl(mod, ~ glance(.x)$r.squared),
      y = map2_dbl(mod, x, ~ predict(.x, newdata = setNames(data.frame(.y), bac_col))),
      label = paste0(Microbiome, "\nR² = ", round(r2, 2), "\nρ = ", round(spearman_rho, 2))
    )
  
  ggplot(data, aes(x = .data[[bac_col]], y = .data[[bif_col]], color = Microbiome)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
    geom_label_repel(
      data = model_info,
      aes(x = x, y = y, label = label, fill = Microbiome),
      color = "black", show.legend = FALSE,
      size = 3, label.size = 0.3, max.overlaps = 100
    ) +
    theme_bw() +
    guides(color = 'none', fill = 'none') +
    labs(x = paste0(bac_col, " (log2 CFU)"), y = paste0(bif_col, " (log2 CFU)"))
}

Bif_Species <- c("Bif_adolescentis_log2", "Bif_breve_log2", 
                 "Bif_catenulatum_log2",  "Bif_bifidum_log2"  )

Bac_Species <- c("Bac_cutis_log2","Bac_timonensis_log2","Bac_kribbi_log2",              
                 "Bac_faecichinchillae_log2",          
                 "Bac_caccae_log2","Bac_ovatus_log2",              
                 "Bac_acidifaciens_log2",       
                 "Bac_fragilis_log2","Bac_thetaiotaomicron_log2",    
                 "Bac_intestinalis_log2",             
                 "Bac_cellulosilyticus_log2","Bac_clarus_log2" ,             
                 "Bac_plebeius_log2","Bac_coprocola_log2",           
                 "Bac_stercoris_log2","Bac_eggerthii_log2",          
                 "Bac_salyersiae_log2","Bac_nordii_log2",              
                 "Bac_uniformis_log2",
                 "Bac_massiliensis_log2",        
                 "Bac_vulgatus_log2","Bac_dorei_log2")

for (i in Bif_Species) {
  for (j in Bac_Species) {
    print(paste("Plotting", i, "vs", j))
    print(B_B_Correlation(Bac_Bif_Genus, i, j))
  }
}










