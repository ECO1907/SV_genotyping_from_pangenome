###This script is used to calculate the LD density and to plot the results

#1. Load the necessary packages

library('tidyverse')
library('ggplot2')
library('dplyr')
library("RColorBrewer")

#2. Calculate the mean LD for 1000 bp bins for the 3 scenarios

##a. For SNPs-SNPs scenario

ld_snps <- read_table("/path/to/SNP_SNP.ld", col_types = cols()) %>%
  mutate(distance = abs(BP_B - BP_A))
ld_snps$distance <- as.numeric(ld_snps$distance)

ld_bins_snps <- ld_snps %>%
  group_by(cut(distance, breaks = seq(0, max(distance), by = 1000))) %>%
  summarize(mean_R2 = mean(R2), distance = mean(distance))

ld_snps_100k <- ld_bins_snps[ld_bins_snps$distance < 100000,] #Cut the maximum distance to 100 Kb

##b. For SVs-SVs scenario

ld_svs <- read_table("/path/to/SV_SV.ld", col_types = cols()) %>%
  mutate(distance = abs(BP_B - BP_A))
ld_svs$distance <- as.numeric(ld_svs$distance)

ld_bins_svs <- ld_svs %>%
  group_by(cut(distance, breaks = seq(0, max(distance), by = 1000))) %>%
  summarize(mean_R2 = mean(R2), distance = mean(distance))

ld_svs_100k <- ld_bins_svs[ld_bins_svs$distance < 100000,] #Cut the maximum distance to 100 Kb

##c. For SVs-SNPs scenario

ld_all <- read_table("/path/to/SNP_SV.ld", col_types = cols()) %>%
  mutate(distance = abs(BP_B - BP_A))
ld_all$distance <- as.numeric(ld_all$distance)

ld_bins_all <- ld_all %>%
  group_by(cut(distance, breaks = seq(0, max(distance), by = 1000))) %>% # Grouping and averaging
  summarize(mean_R2 = mean(R2), distance = mean(distance))

ld_all_100k <- ld_bins_all[ld_bins_all$distance < 100000,] #Cut the maximum distance to 100 Kb

#3. Combine the datasets in a single one indicating to which scenario belongs

##a. Indicate the scenarios for each dataset

ld_snps_100k$Scenario <- "SNPs-SNPs"
ld_svs_100k$Scenario <- "SVs-SVs"
ld_all_100k$Scenario <- "SVs-SNPs"

##b. Combine the datasets

ld_combined <- rbind(ld_svs_100k, ld_snps_100k, ld_all_100k)

#4. Plot the results

##a. Set the plot colors

colors <- brewer.pal(8, "Dark2")[c(4,3,1)]

##b. Plot the results

ggplot(ld_combined, aes(x = mean_R2, fill = Scenario, colour = Scenario)) +
  geom_density(alpha = 0.30, linewidth = 0.5) +
  labs(x = "Mean LD (R²)", y = "") +
  scale_fill_manual(values = colors) +   # fill colours
  scale_color_manual(values = colors) +  # outline colours
  theme_minimal()

#5. Clean the environment

rm(list = ls())
gc()
