###This script is used to compare the phenotype values between the accessions that carry candidate variants and the accessions that don't###

#1. Load the necessary packages

library('data.table')
library('ggplot2')
library('dplyr')
library('ggpubr')

#2. Load the SVs+SNPs dataset, the map and the phenotype files

geno <- fread("/path/to/both_japonica_transversed.txt", header = FALSE) #Only the both dataset is used as this contains all the variants

map <- fread("/path/to/both_japonica_map_imputed.txt", header = TRUE)
variants <- map$SNP
colnames(geno) <- variants

pheno <- fread("/path/to/phenotype.txt")

#3. Create data frame for a given trait

candidate <- pheno %>% select(taxa, trait) #A given trait should be indicated from pheno
colnames(candidate) <- c('Sample', "Value")

#4. Indicate the given candidate variant

candidate_variant <- "candidate_variant" #A given variant should be indicated ex. DEL:3:7852323..7852367 or SNP:5:5392530

#5. Obtain only the data for the candidate variant

candidate$genotype <- geno[, candidate_variant]

candidate$presence <- ifelse(candidate$genotype == 2, '+', '-')

#6. Plot the phenotypic values for both group of accessions

ggplot(candidate, aes(x = presence, y = Value)) +
  geom_violin(fill = 'white', color = 'black', alpha = 0.7) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.3) +
  stat_summary(fun = mean, geom = "crossbar", size = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 1)), vjust = -0.5, hjust = 2, color = "black", size = 3.5) +
  stat_compare_means(method = "t.test", label = "p.format", comparisons = list(c("+", "-")), label.y = max(candidate$Value, na.rm = TRUE)) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 16, face = "bold"),  # Increase size & bold
    axis.text.y = element_text(size = 12)                  # Optional: size for y-axis
  )

#7. Clean the environment

rm(list = ls())
gc()
