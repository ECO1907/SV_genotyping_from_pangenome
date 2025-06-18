###This script is used to impute the genotype datasets (012 files)###

#1. Load the necessary packages

library("data.table")
library("snpReady")

#2. Load the genotype dataset

##a. Load the dataset

var_df <- fread("/path/to/dataset_japonica_recoded.raw")

##b. Get the variants IDs (for processing)

var_ids <- t(fread("/path/to/dataset_japonica.bim", header = FALSE)[,2])

#3. Clean the data frame to just keep the genotypic data

colnames(var_df) <- c("V1", "V2", "V3", "V4", "V5", "V6", var_ids)
var_df_clean <- var_df[,-(1:6)]

#4. Impute the dataset

##a. Impute the dataset

var_df_imputed <- raw.data(as.matrix(var_df_clean), frame = "wide", base = FALSE,
                           sweep.sample = 1, call.rate = 0.9, maf = 0.03, imput = TRUE,
                           imput.type = "wright", outfile = "012", plot = FALSE)

##b. Save the imputed dataset

var_df_ready <- as.data.frame(var_df_imputed$M.clean)
fwrite(var_df_ready, "/path/to/working_directory/dataset_japonica_imputed.txt", sep = "\t")

#5. Clean the environment

rm(list = ls())
gc()
