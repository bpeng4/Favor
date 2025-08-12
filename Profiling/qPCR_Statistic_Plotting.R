library(agricolae)

## Prepare Standard Curve
SDcurve<-read.csv(file = "/Users/bopeng/Documents/Favor/Documents/Bifido Test2 04252025.csv")
sdcurve<-lm(formula = SDcurve$Logarithmic.concentration~SDcurve$CtSYBR.CT.)

plot(x=SDcurve$CtSYBR.CT.,y=SDcurve$Logarithmic.concentration, ylab="log(CFU)", xlab="CT")
curve(expr = sdcurve$coefficients[1]+sdcurve$coefficients[2]*x, 
      col= "blue", add = TRUE)
legend(x = 22, y = 7.5, legend = "log(CFU)= 9.8941 - 0.2967*CT", bty = "n")
legend(x = 22, y = 7, legend = "r = -0.9986166", bty = "n")
cor(SDcurve$CtSYBR.CT.,SDcurve$Logarithmic.concentration)

## Input Data
set1<- read.csv(file = "/Users/bopeng/Documents/Favor/Cache/FavorPericarp071725_All.csv")
#CT to CFU
set1$CT <- as.numeric(set1$CT)
set1 <- set1[!is.na(set1$CT),]
set1$LogCon<-sdcurve$coefficients[1]+sdcurve$coefficients[2]*set1$CT
set1$CFU<-10^(set1$LogCon)
set1$log2<-log2(set1$CFU)

## Statistical Tests
# ANOVA + Duncan
model <- aov(CFU ~ Seed_Name*Microbiome, data = set1)
summary(model)
duncan_result <- duncan.test(model, c("Seed_Name", "Microbiome"), console = TRUE)

write.csv(duncan_result[4], file = "/Users/bopeng/Documents/Favor/Plots/Favor_Pericarp_071525_Duncan_Mean_SE.csv" )
write.csv(duncan_result[6], file = "/Users/bopeng/Documents/Favor/Plots/Favor_Pericarp_071525_Post_Hoc_Comparison.csv" )

## Plotting Preparation
# Load required library
library(ggplot2)
library(dplyr)
library(ggbreak)
library(tidyr)


# Define colors
colors<-c(
  '#F781BF', '#7FC97F', '#FF7F00', '#8DD3C7', '#999999',  
  '#CAB2D6', '#E41A1C', '#386CB0', '#f7754f', '#ce9032',
  '#FFC300', '#B429DC', '#1CF7C2', '#9F4004', '#E9A971',
  '#1132F0', '#86940f',
  '#3ca2f4','#a48cf4')



#Order the Microbiomes
microbiome_order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                      "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", 
                      "S19","S20", "S21", "S22", "S23", "S24")

set1 <- set1 %>%
  filter(Microbiome %in% microbiome_order) %>%
  mutate(Microbiome = factor(Microbiome, levels = microbiome_order)) %>%
  arrange(Microbiome)

#First Set "S9 S15 S18 S21"
#Subset Microbiomes
set1sub <- set1[set1$Microbiome %in% c("S9", "S15", "S18", "S21"),]


model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, "Microbiome", console = TRUE)
duncan_result <- duncan.test(model, "Seed_Name", console = TRUE)


#Order the samples
sample_order <- c( "LH60", "PHG47", "NC262", "B57", "HP301", "PHG39", "A554", 
                   "764", "PHG35",  "FBB0")


set1sub <- set1sub %>%
  filter(Seed_Name %in% sample_order) %>%
  mutate(Seed_Name = factor(Seed_Name, levels = sample_order)) %>%
  arrange(Seed_Name)

#Duncan Test for the sub-microbiome
model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, c("Seed_Name", "Microbiome"), console = TRUE)

#Prepare grouping labels
symbols <- duncan_result[6] |> as.data.frame()
symbols$Sample <- rownames(symbols)
symbols <- symbols %>%
  separate(Sample, into = c("Seed", "Microbiome"), sep = ":")

microbiome_order <- c("S9","S15","S18", "S21")
symbols <- symbols %>%
  mutate(
    Seed = factor(Seed, levels = sample_order),
    Microbiome = factor(Microbiome, levels = microbiome_order)
  ) %>%
  arrange(Seed, Microbiome)


## Summarizing data and Plotting
bar_data <- set1sub %>%
  group_by(Seed_Name, Microbiome) %>%
  summarise(log2 = mean(log2), .groups = "drop") %>%
  mutate(group_label = symbols$groups.groups)

subcolors <- c('#386CB0', '#FFC300', '#9F4004','#86940f')

ggplot(bar_data, aes(x = Seed_Name, y = log2, fill = Microbiome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = subcolors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("S9 S15 S18 S21") +
  scale_y_break(c(1, 10), scales = 300) +
  labs(y = "Bifidobacteria (Log2(CFU))", x = "Sample Type") +
  geom_text(
    aes(label = group_label),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  )


#Second Set "S5 S14 S20 S24"
#Subset Microbiomes
set1sub <- set1[set1$Microbiome %in% c("S20","S5","S14", "S24"),]

model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, "Microbiome", console = TRUE)
duncan_result <- duncan.test(model, "Seed_Name", console = TRUE)



#Order the samples
sample_order <- c( "PHG47", "764", "NC262", "B57", "PHG35", "HP301",
                   "PHG39", "LH60", "A554", "FBB0")


set1sub <- set1sub %>%
  filter(Seed_Name %in% sample_order) %>%
  mutate(Seed_Name = factor(Seed_Name, levels = sample_order)) %>%
  arrange(Seed_Name)

#Duncan Test for the sub-microbiome
model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, c("Seed_Name", "Microbiome"), console = TRUE)

#Prepare grouping labels
symbols <- duncan_result[6] |> as.data.frame()
symbols$Sample <- rownames(symbols)
symbols <- symbols %>%
  separate(Sample, into = c("Seed", "Microbiome"), sep = ":")

microbiome_order <- c("S5","S14", "S20", "S24")
symbols <- symbols %>%
  mutate(
    Seed = factor(Seed, levels = sample_order),
    Microbiome = factor(Microbiome, levels = microbiome_order)
  ) %>%
  arrange(Seed, Microbiome)


## Summarizing data and Plotting
bar_data <- set1sub %>%
  group_by(Seed_Name, Microbiome) %>%
  summarise(log2 = mean(log2), .groups = "drop") %>%
  mutate(group_label = symbols$groups.groups)

subcolors <- c('#999999', '#ce9032', '#1132F0', '#a48cf4')

ggplot(bar_data, aes(x = Seed_Name, y = log2, fill = Microbiome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = subcolors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("S5 S14 S20 S24") +
  scale_y_break(c(1, 10), scales = 300) +
  labs(y = "Bifidobacteria (Log2(CFU))", x = "Sample Type") +
  geom_text(
    aes(label = group_label),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  )

# Third Set "S6","S12", "S16","S22"
#Subset Microbiomes
set1sub <- set1[set1$Microbiome %in% c("S6","S12", "S16","S22"),]


model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, "Microbiome", console = TRUE)
duncan_result <- duncan.test(model, "Seed_Name", console = TRUE)



#Order the samples
sample_order <- c( "HP301", "PHG47", "PHG35", "764", "NC262", "A554", "PHG39",
                   "B57", "LH60", "FBB0")


set1sub <- set1sub %>%
  filter(Seed_Name %in% sample_order) %>%
  mutate(Seed_Name = factor(Seed_Name, levels = sample_order)) %>%
  arrange(Seed_Name)

#Duncan Test for the sub-microbiome
model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, c("Seed_Name", "Microbiome"), console = TRUE)

#Prepare grouping labels
symbols <- duncan_result[6] |> as.data.frame()
symbols$Sample <- rownames(symbols)
symbols <- symbols %>%
  separate(Sample, into = c("Seed", "Microbiome"), sep = ":")

microbiome_order <- c("S6 S12 S16 S22")
symbols <- symbols %>%
  mutate(
    Seed = factor(Seed, levels = sample_order),
    Microbiome = factor(Microbiome, levels = microbiome_order)
  ) %>%
  arrange(Seed, Microbiome)


## Summarizing data and Plotting
bar_data <- set1sub %>%
  group_by(Seed_Name, Microbiome) %>%
  summarise(log2 = mean(log2), .groups = "drop") %>%
  mutate(group_label = symbols$groups.groups)

subcolors <- c('#CAB2D6', '#f7754f', '#B429DC', '#3ca2f4')

ggplot(bar_data, aes(x = Seed_Name, y = log2, fill = Microbiome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = subcolors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("S6 S12 S16 S22") +
  scale_y_break(c(1, 10), scales = 300) +
  labs(y = "Bifidobacteria (Log2(CFU))", x = "Sample Type") +
  geom_text(
    aes(label = group_label),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  )


#Fourth Set "S2","S19","S17"
#Subset Microbiomes
set1sub <-set1[set1$Microbiome %in% c("S2","S19","S17"),]

model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, "Microbiome", console = TRUE)
duncan_result <- duncan.test(model, "Seed_Name", console = TRUE)



#Order the samples
sample_order <- c( "PHG47", "764", "PHG35", "NC262", "LH60", "HP301",
                   "B57",  "A554", "PHG39", "FBB0")

set1sub <- set1sub %>%
  filter(Seed_Name %in% sample_order) %>%
  mutate(Seed_Name = factor(Seed_Name, levels = sample_order)) %>%
  arrange(Seed_Name)

#Duncan Test for the sub-microbiome
model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, c("Seed_Name", "Microbiome"), console = TRUE)

#Prepare grouping labels
symbols <- duncan_result[6] |> as.data.frame()
symbols$Sample <- rownames(symbols)
symbols <- symbols %>%
  separate(Sample, into = c("Seed", "Microbiome"), sep = ":")

microbiome_order <- c()
symbols <- symbols %>%
  mutate(
    Seed = factor(Seed, levels = sample_order),
    Microbiome = factor(Microbiome, levels = microbiome_order)
  ) %>%
  arrange(Seed, Microbiome)


## Summarizing data and Plotting
bar_data <- set1sub %>%
  group_by(Seed_Name, Microbiome) %>%
  summarise(log2 = mean(log2), .groups = "drop") %>%
  mutate(group_label = symbols$groups.groups)

subcolors <- c('#7FC97F', '#1CF7C2', '#E9A971')

ggplot(bar_data, aes(x = Seed_Name, y = log2, fill = Microbiome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = subcolors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("S2 S17 S19") +
  scale_y_break(c(1, 10), scales = 300) +
  labs(y = "Bifidobacteria (Log2(CFU))", x = "Sample Type") +
  geom_text(
    aes(label = group_label),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  )


# Fifth Set
#Subset Microbiomes
set1sub <- set1[set1$Microbiome %in% c("S1","S3", "S7"),]

model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, "Microbiome", console = TRUE)
duncan_result <- duncan.test(model, "Seed_Name", console = TRUE)



#Order the samples
sample_order <- c( "764", "PHG39","PHG47", "B57", "PHG35", "LH60", "A554", 
                   "NC262", "HP301", "FBB0")


set1sub <- set1sub %>%
  filter(Seed_Name %in% sample_order) %>%
  mutate(Seed_Name = factor(Seed_Name, levels = sample_order)) %>%
  arrange(Seed_Name)

#Duncan Test for the sub-microbiome
model <- aov(CFU ~ Seed_Name*Microbiome, data = set1sub)
summary(model)
duncan_result <- duncan.test(model, c("Seed_Name", "Microbiome"), console = TRUE)

#Prepare grouping labels
symbols <- duncan_result[6] |> as.data.frame()
symbols$Sample <- rownames(symbols)
symbols <- symbols %>%
  separate(Sample, into = c("Seed", "Microbiome"), sep = ":")

microbiome_order <- c("S1","S3", "S7")
symbols <- symbols %>%
  mutate(
    Seed = factor(Seed, levels = sample_order),
    Microbiome = factor(Microbiome, levels = microbiome_order)
  ) %>%
  arrange(Seed, Microbiome)


## Summarizing data and Plotting
bar_data <- set1sub %>%
  group_by(Seed_Name, Microbiome) %>%
  summarise(log2 = mean(log2), .groups = "drop") %>%
  mutate(group_label = symbols$groups.groups)

subcolors <- c('#F781BF', '#FF7F00', '#E41A1C')

ggplot(bar_data, aes(x = Seed_Name, y = log2, fill = Microbiome)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = subcolors) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    legend.title = element_blank()
  ) +
  ggtitle("S1 S3 S7") +
  scale_y_break(c(1, 10), scales = 300) +
  labs(y = "Bifidobacteria (Log2(CFU))", x = "Sample Type") +
  geom_text(
    aes(label = group_label),
    position = position_dodge(width = 0.8),
    vjust = -0.5,
    size = 3.5
  )



