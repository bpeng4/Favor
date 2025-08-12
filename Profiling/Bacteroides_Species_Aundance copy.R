library(GGally)

B_Species <- read.csv(file = "/Users/bopeng/Documents/Favor/Cache/Bacteroides_Species.csv")

B_Species_Long <- reshape(B_Species, direction = "long",
                          varying = c( "B.caccae", "B.ovatus", "B.acidifaciens",
                                       "B.fragilis", "B.thetaiotaomicron",
                                       "B.cellulosilyticus", "B.stercoris",
                                       "B.eggerthii","B.uniformis", "B.dorei",       
                                       "B.massiliensis","B.vulgatus"),
                          timevar = "Species", idvar = "Sample_Name")
B_Species_Long$RA <- B_Species_Long$B

B_Species_Long$Combine <- paste(B_Species_Long$Seed_Type, B_Species_Long$Reps,
                                 B_Species_Long$Species, sep = "-") 

B_Species_Wide <- reshape(B_Species_Long, direction = "wide",
                          idvar = "Combine", timevar = "Microbiome",
                          v.names = "RA")

#replace all NA values to 0
B_Species_Wide[is.na(B_Species_Wide)] <- 0

### Modify the microbiome names
# Assume your dataframe is named df
colnames(B_Species_Wide) <- gsub("^RA\\.", "", colnames(B_Species_Wide))

# Reorder the columns
# Get the column names from 7 to 25
cols_7_25 <- names(B_Species_Wide)[7:25]

# Wanted Order
microbiome_order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                      "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", 
                      "S19","S20", "S21", "S22", "S23", "S24")

# Reorder those columns using your custom vector
ordered_cols <- microbiome_order[microbiome_order %in% cols_7_25]

# Combine with the rest of the columns
B_Species_Wide <- B_Species_Wide[, c(names(B_Species_Wide)[1:6], ordered_cols )]


colors<-c(
  '#F781BF', '#7FC97F', '#FF7F00', '#8DD3C7', '#999999',  
  '#CAB2D6', '#E41A1C', '#386CB0', '#ce9032',
  '#FFC300', '#B429DC', '#1CF7C2', '#9F4004', '#E9A971',
  '#1132F0', '#86940f',
  '#3ca2f4','#a48cf4')

ggparcoord(data = B_Species_Wide,
           columns = 7:25,
           alphaLines = 0.2,
           boxplot = TRUE,
           groupColumn = "Species",
           scale = 'std') +
  scale_color_manual(values = colors)+
  facet_wrap(~ Species)




### Standarized by FBB0
#Set up FBB0 as the starting points
library(dplyr)
B_Species_Long0<-B_Species_Long

fbb0_df <- B_Species_Long0 %>%
  filter(Seed_Type == "FBB0") %>%
  select(Microbiome, Reps, Species, RA) %>%
  mutate(RA_start = RA) %>%
  select(-RA)  # Remove old columns to avoid duplicates

B_Species_Long0 <- B_Species_Long0 %>%
  left_join(fbb0_df, by = c("Microbiome", "Reps", "Species"))  

B_Species_Long0$RA_adjust <- B_Species_Long0$RA - B_Species_Long0$RA_start

B_Species_Long0 <- B_Species_Long0[B_Species_Long0$Seed_Type !="FBB0",]

B_Species_Wide0 <- reshape(B_Species_Long0, direction = "wide",
                          idvar = "Combine", timevar = "Microbiome",
                          v.names = "RA_adjust")

#replace all NA values to 0
B_Species_Wide0[is.na(B_Species_Wide0)] <- 0

### Modify the microbiome names
# Assume your dataframe is named df
colnames(B_Species_Wide0) <- gsub("^RA_adjust\\.", "", colnames(B_Species_Wide0))

# Reorder the columns
# Get the column names from 9 to 27
cols_9_27 <- names(B_Species_Wide0)[9:27]

# Wanted Order
microbiome_order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                      "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", 
                      "S19","S20", "S21", "S22", "S23", "S24")

# Reorder those columns using your custom vector
ordered_cols <- microbiome_order[microbiome_order %in% cols_9_27]

# Combine with the rest of the columns
B_Species_Wide0 <- B_Species_Wide0[, c(names(B_Species_Wide0)[1:8], ordered_cols )]

#Plotting
colors<-c(
  '#F781BF', '#7FC97F', '#FF7F00', '#8DD3C7', '#999999',  
  '#CAB2D6', '#E41A1C', '#386CB0', '#ce9032',
  '#FFC300', '#B429DC', '#1CF7C2', '#9F4004', '#E9A971',
  '#1132F0', '#86940f',
  '#3ca2f4','#a48cf4')


ggparcoord(data = B_Species_Wide0,
           columns = 9:27,
           alphaLines = 0.2,
           boxplot = TRUE,
           groupColumn = "Species",
           scale = 'std') +
  scale_color_manual(values = colors)+
  facet_wrap(~ Species, ncol =2)

# Write a function for plotting
Bac_Species <- function(B_Species_Wide, Seed) {
  ggparcoord(data = B_Species_Wide,
             columns = 9:27,
             alphaLines = 0.2,
             boxplot = TRUE,
             groupColumn = "Species",
             scale = 'std') +
    scale_color_manual(values = colors) +
    facet_wrap(~ Species, ncol = 2) +
    labs(x = "Microbiomes", y = "Bacteroides Species Standardized Read Counts") +
    ggtitle(paste(Seed , "-FBB0"))
}

# Subset by Seed Type and generate plots
Seed_Types <- unique(B_Species_Wide0$Seed_Type)

for (i in Seed_Types) {
  B_Species_Wide0_sub <- B_Species_Wide0[B_Species_Wide0$Seed_Type == i, ]
  
  p <- Bac_Species(B_Species_Wide0_sub, i)
  print(p)  # Ensure the plot is displayed
}



#Plot for Bacteroides celluloilyticus
B_cellulosilyticus <- B_Species_Wide0[B_Species_Wide0$Species == "cellulosilyticus",]

ggparcoord(data = B_cellulosilyticus,
           columns = 9:27,
           alphaLines = 0.2,
           boxplot = TRUE,
           groupColumn = "Seed_Type",
           scale = 'std') +
  scale_color_manual(values = colors)+
  facet_wrap(~ Seed_Type, ncol =2)+
  labs(x = "Microbiomes", y = "(Seed Types - FBB0) Standardized Read Counts") +
  ggtitle(paste("Bacteroides cellulosilyticus" ))
