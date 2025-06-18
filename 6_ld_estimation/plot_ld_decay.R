###This script is used to calculate the LD decay and plot it###

#1. Load the necessary packages

library('tidyverse')
library('ggplot2')
library('dplyr')

#2. Load the LD scenarios results and calculate the distance between the variants

ld_data <- read_table("/path/to/scenario.ld", col_types = cols()) %>%
  mutate(distance = abs(BP_B - BP_A)) #Calculate the distance between the pair of variants
ld_data$distance <- as.numeric(ld_data$distance) #Assure that the distance value is numeric

#3. Calculate mean LD for 1000 bp bins

ld_bins <- ld_data %>%
  group_by(cut(distance, breaks = seq(0, max(distance), by = 1000))) %>% #Group the variants in 1000 bp bins
  summarize(mean_R2 = mean(R2), distance = mean(distance)) #Calculate the mean LD for each bin. Create a sequence from 0 to the maximum distance in steps of 1000.

#4. Calculate the half maximum mean LD for bins and the distance where it drops to this value

max_ld <- max(ld_bins$avg_R2)
min_ld <- min(ld_bins$avg_R2)
half_ld <- (max_ld+min_ld)/2

closest_row <- ld_bins[which.min(abs(ld_bins$avg_R2 - half_ld)), ]
x_at_half_ld <- closest_row$distance

#5. Plot the LD decay curve

ggplot(ld_bins, aes(x = distance, y = mean_R2)) +
  geom_line(color = "black") +
  geom_hline(yintercept = half_ld, linetype = "dashed", color = "red") +
  annotate("text", x = max(ld_bins$distance) * 0.95, y = half_ld, label = paste("R² =", round(half_ld,2)), color = "red", hjust = 1, vjust = -0.5, size = 4.5) +
  geom_vline(xintercept = x_at_half_ld, linetype = "dotted", color = "blue") +
  annotate("text", x = x_at_half_ld, y = 0.20, label = paste(round(x_at_half_ld,0), "bp"), color = "blue", angle = 90, vjust = -0.5, size = 4.5) +
  scale_x_continuous(limits = c(0, 400000), breaks = seq(0, 400000, 100000), labels = scales::comma) +
  ylim(0, 1) +
  labs(x = "Distance (bp)", y = "Mean LD (R²)") +
  theme_minimal(base_size = 14)

#6. Clean the environment

rm(list = ls())
gc()
