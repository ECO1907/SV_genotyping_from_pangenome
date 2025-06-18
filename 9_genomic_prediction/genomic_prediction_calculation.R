###This script is used to calculate Genomic Prediction accuracy using RKHS###

#1. Load the necessary packages

library("data.table")
library("BGLR")
library("gridExtra")
library("dplyr")

#2. Load phenotype data

phenotype <- fread("/path/to/phenotype.txt", header = TRUE)
traits <- names(phenotype)[2:9]  # Select traits

#3. Load genotype dataset

dataset <- as.matrix(fread('/path/to/dataset_japonica_transversed.txt', header = FALSE))

#4. Load kinship matrix

kin_dataset <- as.matrix(fread("path/to/dataset_japonica_kin.txt"))

#5. Load the PCA

pca_dataset <- fread("/path/to/dataset_japonica_pca.txt")
dataset_cov <- pca_dataset[,1:3]

#6. Set the working directory

setwd("/path/to/working_directory")

#7. Calculate genomic prediction accuracy

nIter = 100000
results_list <- list()

for (rep in 1:10) {
  set.seed(100 + rep)
  n_samples <- nrow(phenotype)
  test_indices <- sample(1:n_samples, size = round(0.25 * n_samples))
  
  sample_ids <- phenotype[[1]]
  test_samples <- sample_ids[test_indices]
  write.table(test_samples, file = paste0("test_samples_rep", rep, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  all.phenotypes.COR <- NULL 
  all.phenotypes.MSE <- NULL 
  
  pdf(file=paste0("rkhs_dataset_rep", rep, ".pdf"), onefile=TRUE)
  par(mfrow=c(2,2))
  
  fm = list()
  for (i in 1:length(traits)) {
    y <- scale(phenotype[[traits[i]]], center=TRUE, scale=TRUE)
    yNA <- y
    yNA[test_indices] <- NA
    
    fm2 <- BGLR(y=yNA, 
                ETA=list(ETA1=list(X=dataset_cov, model="FIXED"), 
                         ETA2=list(K=kin_dataset, model='RKHS')), 
                nIter=nIter, 
                saveAt=sprintf("./rep%d_fm_RKHS_dataset_y%.0f_", rep, i), 
                verbose=FALSE)
    
    yHat <- fm2$yHat
    y_test <- y[test_indices]
    yHat_test <- yHat[test_indices]
    
    corel2 <- cor(yHat_test, y_test, use="complete.obs")
    mse <- mean((yHat_test - y_test)^2, na.rm=TRUE)
    
    all.phenotypes.COR <- rbind(all.phenotypes.COR, corel2)
    all.phenotypes.MSE <- rbind(all.phenotypes.MSE, mse)
    
    tmp <- range(c(yHat_test, y_test), na.rm=TRUE) 
    plot(yHat_test ~ y_test, col="red", xlim=tmp, ylim=tmp,
         xlab="Observed", ylab="Predicted",
         main=paste(traits[i], sprintf("Rep %d RKHS cor=%.02f", rep, corel2)))
    abline(lm(yHat_test ~ y_test))
    
    fm[[i]] <- fm2
  }
  dev.off()
  
  accuracy_results <- data.frame(Replicate=rep, Trait=traits, 
                                 Correlation=all.phenotypes.COR, 
                                 MSE=all.phenotypes.MSE)
  
  results_list[[rep]] <- accuracy_results
}

#8. Obtain a whole results table

final_results <- do.call(rbind, results_list)
fwrite(final_results, "RKHS_dataset_Accuracy_10reps.txt", sep = "\t", row.names=FALSE)

#9. Obtain a summary report

summary_stats <- final_results %>%
  group_by(Trait) %>%
  summarise(Correlation_Mean = mean(Correlation),
            Correlation_SD = sd(Correlation),
            MSE_Mean = mean(MSE),
            MSE_SD = sd(MSE))

fwrite(summary_stats, "RKHS_dataset_Accuracy_Summary.txt", sep = "\t", row.names=FALSE)

#10. Clean the environment

rm(list = ls())
gc()
