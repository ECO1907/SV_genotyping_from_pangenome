###This script is used to calculate the herediability explained by dataset with RHKS from the BGLR package###

#1. Load the necessary packages

library("BGLR")
library("data.table")

#2. Load the phenotype data (this contains the information of all the traits)

phenotype <- fread("/path/to/phenotype.txt", header = TRUE)
data <- cbind(phenotype[,2:9])
list_phenotype <- list(phenotype[,2:9])
traits <- names(data)

#3. Load the genotype dataset, kinship matrix and pca

dataset <- as.matrix(t(fread("/path/to/dataset_japonica_transversed.txt", header = FALSE)))

kin_dataset <- as.matrix(fread("/path/to/dataset_japonica_kin.txt"))

pca_dataset <- fread("/path/to/dataset_japonica_pca.txt")
dataset_cov <- pca_dataset[,1:3]

#4. Compute the heritability

##a. Define the number of iterations

nIter = 100000

##b. Create a data frame to store the results

heritability_dataset <- data.frame(Trait=colnames(data), h2_dataset=NA, varE_dataset=NA)

##c. Compute the heritability

for (i in 1:8) {
  
  y <- scale(data[[i]], center=TRUE, scale=TRUE) #Extract and scale the phenotype data
  
  fm <- BGLR(y = y,
             ETA = list(
               ETA1 = list(K = kin_dataset, model = 'RKHS'), #Kinship matrix
               PCA = list(X = dataset_cov, model = "FIXED") #PCA matrix
             ),
             nIter = nIter,
             verbose = TRUE)
  
  varG <- scan(sprintf("ETA_ETA1_varU.dat")) #Store the Genetic variance
  varE <- scan(sprintf("varE.dat")) #Store the residual variance
  
  h2 <- mean(varG/(varG+varE)) #Calculate the heritability
  
  heritability_dataset[i, "h2_dataset"] <- h2 #Store the heritability in the created data frame
  heritability_dataset[i, "varE_dataset"] <- mean(varE) #Store the residual variance in the created data frame
  
}

#6.Save the results in the working directory

setwd("/path/to/working_directory")

fwrite(heritability_dataset, "heritability_dataset.txt", sep = "\t", row.names = FALSE)

#7. Clean the environment

rm(list = ls())
gc()
