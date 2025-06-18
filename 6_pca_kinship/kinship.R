###This script is used to obtain a kinship matrix for the imputed datasets###

#1. Load the necessary packages

library("data.table")
library("AGHmatrix")
library("pheatmap")

#2. Load the imputed genotype file and assign the accession information

##a. Load the imputed dataset

var_df <- fread("/path/to/dataset_japonica_imputed.txt")

##b. Load the japonica accessions information

samples <- fread("/path/to/japonica_accessions_info.txt", header = FALSE)
colnames(samples) <- c("Sample", "K_group")

#3. Run the Kinship analysis

##a. Run the kinship analysis

kinship_matrix <- Gmatrix(SNPmatrix = as.matrix(var_df), method = "VanRaden", missingValue = NA)

##c. Save the results

fwrite(kinship_matrix, "/path/to/working_directory/dataset_japonica_kin.txt", sep = "\t", col.names = FALSE, row.names = FALSE)

#4. Plot the results

##a. Assign the samples ID to the kinship matrix

colnames(kinship_matrix) <- samples$Sample
rownames(kinship_matrix) <- samples$Sample

##b. Set the limits of plot axis

max_abs <- 1
breaks <- seq(-max_abs, max_abs, length.out = 101)

##c. Plot the results

pheatmap(kinship_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = breaks,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")

#5. Clean the environment

rm(list = ls())
gc()

