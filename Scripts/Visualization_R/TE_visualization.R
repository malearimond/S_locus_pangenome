### Visualization of the S-locus and repeat content across S-haplotypes ###
## 02.04.2025 ##

## load packages ##

library(ggplot2)
library(ggforce)
library(tidyr)
library(ggpubr)

### compute S-locus length with samtools faidx across fasta file ###
Length_of_slocus <- read.delim("~/Path/to/Length_of_slocus.txt")

Length_of_slocus$length <- as.numeric(Length_of_slocus$length..kb.)

# histogram of S-locus length #
ggplot(Length_of_slocus, aes(x = length, fill = Mating_system)) +
  geom_histogram(binwidth = 50, color = "black", alpha = 0.7) +
  labs(title = "Length of S-haplotypes across SI and SC individuals",
       x = "Length [kb]",
       y = "Count",
       fill = "Mating System")+
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16), 
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),
    legend.position = "top" # Remove legend (if you want it visible, change this line)
  )+
  scale_fill_manual(values = c("SI" = "#0072B2", "SC" = "brown"))



#### Repeat and TE annotation across S-locus of SC and SI ####
## set shape settings for mating cluster ##
shapes_filter_mating <- c(
  'Poland' = 0, 'Czech' = 2, 'Pyrenees' = 3, 
  'LaPalma' = 4, 'Main' = 5, 'Pyrenees_OG' = 6, 'Main_OG' = 8, 
  'Czech_OG' = 9, 'LaPalma_OG' = 10, 'Poland_OG' = 11)
# Read in TE data, inferred from RepeatMasker 

## read in data ##
Repeat_annotation <- read.delim("~/path/to/Results/Repeat_annotation.txt")

## remove all NA files ##
Repeat_annotation <- Repeat_annotation[!is.na(Repeat_annotation$Mating_cluster), ]

Repeat_annotation$Mating_system <- as.factor(Repeat_annotation$Mating_system)
Repeat_annotation$Mating_cluster <- as.factor(Repeat_annotation$Mating_cluster)


## plot per SC and SI individuals ##

Repeat_per_mating <- ggplot(Repeat_annotation, aes(x = Mating_system, y = Total_repeats, 
                                                   fill = Mating_system)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +  # Boxplots for SI and SC
  geom_jitter(aes(shape = Mating_cluster), position = position_jitter(0.2), alpha = 1, size=8) + # Points get unique shapes
  labs(title = "Total coverage of repeats across S-haplotypes by mating system",
       x = "Mating system",
       y = "Total Repeats [%]",
       fill = "Mating system",
       shape = "Mating cluster") + # Legend for shapes
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("SC" = "brown", "SI" = "#0072B2")) + # Colors for boxplots
  scale_shape_manual(values = shapes_filter_mating) +  # Assign symbols to clusters
  guides(fill = "none") +  # Remove only fill legend
  theme(
    axis.line = element_line(color = "black"),  
    axis.ticks = element_line(color = "black"),  
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16), 
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank()
  )



## define order across mating cluster ##
Repeat_annotation$Mating_cluster_new <- factor(Repeat_annotation$Mating_cluster, 
                                           levels = c("Main", "Main_OG", "Pyrenees", "Pyrenees_OG", 
                                                      "LaPalma", "LaPalma_OG", "Czech", "Czech_OG",  
                                                      "Poland", "Poland_OG"))

Repeat_per_cluster <-ggplot(Repeat_annotation, aes(x = Mating_cluster_new, y = Total_repeats, fill = Mating_system)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.2), alpha = 0.7) + # Add points for visibility
  labs(title = "Total coverage of repeats across S-haplotypes by mating cluster",
       x = "Mating cluster",
       y = "Total Repeats [%]",
       fill = "Mating system") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_manual(values = c("SC" = "brown", "SI" = "#0072B2"))+
  theme(
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Add axis ticks
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    text = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16), 
    axis.ticks.length = unit(.25, "cm"),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),
    legend.position = "top" # Remove legend (if you want it visible, change this line)
  )


## circular plot of each TE family across S-locus ##

colnames(Repeat_annotation) <- c("Genome_ID" ,"Total_repeats", "SINE","SINE elements" ,"LINE", "LINE elements","LTR","LTR elements","DNA elements","DNA.TE elements" , "Unclassified", "Unclassified elements", "Mating_system","Mating_cluster" , "Mating_cluster_new") 

Repeat_annotation_for_circ <- Repeat_annotation[Repeat_annotation$Genome_ID %in% c("PL01_005" ,"CZ01_005" ,"ES27_024","ES03_014", "ES20_028", "GR01_002_hap1", "IT19_019_hap2", "IT20_003_hap2" ,"IT19_019_hap1", "IT19_002_hap1"),]
long_data <- pivot_longer(Repeat_annotation_for_circ, 
                          cols = c("SINE", "LINE", "LTR", "DNA elements", "Unclassified"), 
                          names_to = "TE family", values_to = "Coverage")
Circulae_TEs <-ggplot(long_data, aes(x = Mating_cluster, y = Coverage, fill = `TE family`)) +
  geom_col(position = "stack") + 
  coord_polar() +
  labs(title = "TE families across S-haplotypes",
       x = "", y = "Coverage [%]") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1, size = 14)) +
  scale_fill_manual(values = c("SINE" = "#E69F00", "LINE" = "#56B4E9", 
                               "LTR" = "#009E73", "DNA elements" = "#F0E442", 
                               "Unclassified" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 10)) + # Grid every 10%
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank()   # Remove y-axis ticks
  )
ggarrange(Repeat_per_mating, Circulae_TEs,
          ncol=1, nrow=2)

