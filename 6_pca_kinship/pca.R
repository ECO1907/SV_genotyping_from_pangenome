##This script is used to run PCA with the imputed genotype dataset###

#1. Load the necessary packages

library("data.table")
library("ggplot2")
library("ggfortify")
library("RColorBrewer")

#2. Load the imputed genotype file and assign the groups

##a. Load the imputed dataset

var_df <- fread("/path/to/dataset_japonica_imputed.txt")

##b. Load the japonica accessions information

samples <- fread("/path/to/japonica_accessions_info.txt", header = FALSE)
colnames(samples) <- c("Sample", "K_group")
samples$Group <- ifelse(samples$K_group == 2, "Medium", ifelse(samples$K_group == 4, "Large", "Admixture"))

##c. Combine both dataframes

final_var_df <- cbind(samples, var_df)

#3. Run the PCA

##a. Prepare the dataset (should only contain genotype data [012])

pca_var <- final_var_df[,-(1:3)]

##b. Run the PCA

pca_res <- prcomp(pca_var, scale. = FALSE)

##c. Save the PCA results

pca_df <- as.data.frame(pca_res$x)
fwrite(pca_df, "/path/to/working_directory/dataset_japonica_pca.txt", sep = "\t")

#4. Plot the results

##a. Set the plot colors

colors <- brewer.pal(8, "Dark2")[c(3,1,4)]

##b. Plot the results

autoplot(pca_res, data = final_var_df) +
  geom_point(aes(color = Group), size = 3) +
  stat_ellipse(level = 0.95, aes(color = Group), size = 1) +
  scale_color_manual(values = colors) +
  theme_classic()

#5. Clean the environment

rm(list = ls())
gc()
