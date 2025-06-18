###This script is used to run GWAS with the BLINK method###

#1. Load the necessary packages and functions

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

library("BLINK")
library("data.table")
library("ggfortify")

#2. Load the genetic transversed dataset, the genetic map and the PCA results 

myGD = read.big.matrix("/path/to/dataset_japonica_transversed.txt", header=FALSE, sep="\t", type="char")

myGM = read.table("/path/to/dataset_japonica_map.txt", header = TRUE)

cov = fread("/path/to/dataset_japonica_pca.txt")
myCV = cov[,1:3] #Only the first 3 PCs are necessary

#2. Run the GWAS for grain width

##a. Load the phenotype file

grain_width = read.table("/path/to/grain_width.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/grain_width")

myBlink_gw = Blink(Y = grain_width, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_gw_blink <- myBlink_gw$GWAS

fwrite(gwas_gw_blink, "grain_width_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_gw <- 0.05 / nrow(myBlink_gw$GWAS)
significant_str_var_gw <- gwas_gw_blink[gwas_gw_blink$P.value < bonf_threshold_gw, ]

fwrite(significant_str_var_gw, "grain_width_signigicant_str_var_blink.txt", sep = "\t")

#3. Run the GWAS for grain length

##a. Load the phenotype file

grain_length = read.table("/path/to/grain_length.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/grain_length/")

myBlink_gl = Blink(Y = grain_length, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_gl_blink <- myBlink_gl$GWAS

fwrite(gwas_gl_blink, "grain_length_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_gl <- 0.05 / nrow(myBlink_gl$GWAS)
significant_str_var_gl <- gwas_gl_blink[gwas_gl_blink$P.value < bonf_threshold_gl, ]

fwrite(significant_str_var_gl, "grain_length_signigicant_str_var_blink.txt", sep = "\t")

#4. Run the GWAS for allelopathy index

##a. Load the phenotype file

allelopathy = read.table("/path/to/allelopathy.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/allelopathy/")

myBlink_ale = Blink(Y = ale, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_al_blink <- myBlink_al$GWAS

fwrite(gwas_al_blink, "allelopathy_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_al <- 0.05 / nrow(myBlink_al$GWAS)
significant_str_var_al <- gwas_al_blink[gwas_al_blink$P.value < bonf_threshold_al, ]

fwrite(significant_str_var_al, "allelopathy_signigicant_str_var_blink.txt", sep = "\t")

#5. Run the GWAS for heading time

##a. Load the phenotype file

days2heading = read.table("/path/to/days2heading.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/days2heading/")

myBlink_dh = Blink(Y = days2heading, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_dh_blink <- myBlink_dh$GWAS

fwrite(gwas_dh_blink, "days2heading_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_dh <- 0.05 / nrow(myBlink_dh$GWAS)
significant_str_var_dh <- gwas_dh_blink[gwas_dh_blink$P.value < bonf_threshold_dh, ]

fwrite(significant_str_var_dh, "days2heading_signigicant_str_var_blink.txt", sep = "\t")

#6. Run the GWAS for grains per panicle

##a. Load the phenotype file

grainsxpanicle = read.table("/path/to/grainsxpanicle.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/grainsxpanicle/")

myBlink_gxp = Blink(Y = grainsxpanicle, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_gxp_blink <- myBlink_gxp$GWAS

fwrite(gwas_gxp_blink, "grainsxpanicle_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_gxp <- 0.05 / nrow(myBlink_gxp$GWAS)
significant_str_var_gxp <- gwas_gxp_blink[gwas_gxp_blink$P.value < bonf_threshold_gxp, ]

fwrite(significant_str_var_gxp, "grainsxpanicle_signigicant_str_var_blink.txt", sep = "\t")

#7. Run the GWAS for height

##a. Load the phenotype file

height = read.table("/path/to/height.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/height/")

myBlink_h = Blink(Y = height, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_h_blink <- myBlink_h$GWAS

fwrite(gwas_h_blink, "height_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_h <- 0.05 / nrow(myBlink_h$GWAS)
significant_str_var_h <- gwas_h_blink[gwas_h_blink$P.value < bonf_threshold_h, ]

fwrite(significant_str_var_h, "height_signigicant_str_var_blink.txt", sep = "\t")

#8. Run the GWAS for n_panicles

##a. Load the phenotype file

n_panicles = read.table("/path/to/n_panicles.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/n_panicles/")

myBlink_np = Blink(Y = n_panicles, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_np_blink <- myBlink_np$GWAS

fwrite(gwas_np_blink, "n_panicles_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_np <- 0.05 / nrow(myBlink_np$GWAS)
significant_str_var_np <- gwas_np_blink[gwas_np_blink$P.value < bonf_threshold_np, ]

fwrite(significant_str_var_np, "n_panicles_signigicant_str_var_blink.txt", sep = "\t")

#9. Run the GWAS for panicle_length

##a. Load the phenotype file

panicle_length = read.table("/path/to/panicle_length.phe", header = TRUE)

##b. Run the GWAS

setwd("/path/to/working_directory/panicle_length/")

myBlink_pl = Blink(Y = panicle_length, GD = myGD, GM = myGM, CV = myCV, maxLoop = 50, time.cal=T) 
gwas_pl_blink <- myBlink_pl$GWAS

fwrite(gwas_pl_blink, "panicle_length_str_var_blink.txt", sep = "\t")

##c. Obtain a file with only the signficant variants

bonf_threshold_pl <- 0.05 / nrow(myBlink_pl$GWAS)
significant_spls_pl <- gwas_pl_blink[gwas_pl_blink$P.value < bonf_threshold_pl, ]

fwrite(significant_spls_pl, "panicle_length_signigicant_str_var_blink.txt", sep = "\t")

#10. Clean the environment

rm(list = ls())
gc()
