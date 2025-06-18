##This script is used to run the GWAS with the GLM method###

#1. Load the necessary packages

library("rMVP")
library("data.table")

#2. Load the imputed genetic dataset and the map file

##a. Load the genetic dataset and transverse it

recoded <- t(fread("/path/to/dataset_japonica_imputed.txt", header = TRUE))
fwrite(recoded, "/path/to/working_directory/dataset_japonica_transversed.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

##b. Create the genetic dataset map and save it

map <- fread("/path/to/str_var_japonica.bim", header = FALSE)
map <- map[,-3]
colnames(map) <- c("CHROM","SNP","POS","REF","ALT")
map <- map[, c("SNP", "CHROM", "POS", "REF", "ALT")]
fwrite(map, "/path/to/working_directory/dataset_japonica_map.txt", sep = "\t")

#3. Generate an MVP genotype file from the tranversed recoded file and the map file

##a. Generate a MVP genotype file

MVP.Data.Numeric2MVP("/path/to/working_directory/dataset_japonica_transversed.txt",
                     out="/path/to/working_directory/dataset_japonica_numeric",
                     map_file = "/path/to/working_directory/dataset_japonica_map.txt",
                     maxLine=1e4,
                     auto_transpose=F)

##b. Load the MVP genotype file as a big matrix

genotype <- attach.big.matrix("/path/to/working_directory/dataset_japonica_numeric.geno.desc")

#4. Generate a MVP kinship file from the kinship matrix

##a. Generate a MVP kinship file

MVP.Data.Kin("/path/to/dataset_japonica_kin.txt",
             out="/path/to/working_directory/dataset_japonica_kinship",
             maxLine=1e4, sep='\t')

##b. Load the MVP kinship file as a big matrix

kinship <- attach.big.matrix("/path/to/working_directory/dataset_japonica_kinship.kin.desc")

#5. Run the GWAS for the traits

##a. Run the GWAS for grain width

grain_width <- read.table("/path/to/grain_width.phe", header = T)

setwd("/path/to/working_directory/grain_width")

imMVP <- MVP(
  phe=grain_width, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##b. Run the GWAS for grain length

grain_length <- read.table("/path/to/grain_length.phe", header = T)

setwd("/path/to/working_directory/grain_length/")

imMVP <- MVP(
  phe=grain_length, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##c. Run the GWAS for allelopathy index

allelopathy <- read.table("/path/to/allelopathy.phe", header = T)

setwd("/path/to/working_directory/allelopathy/")

imMVP <- MVP(
  phe=allelopathy, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##d. Run the GWAS for heading time

days2heading <- read.table("/path/to/days2heading.phe", header = T)

setwd("/path/to/working_directory/days2heading/")

imMVP <- MVP(
  phe=days2heading, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##e. Run the GWAS for grains per panicle

grainsxpanicle <- read.table("/path/to/grainsxpanicle.phe", header = T)

setwd("/path/to/working_directory/grainsxpanicle/")

imMVP <- MVP(
  phe=grainsxpanicle, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##f. Run the GWAS for height

height <- read.table("/path/to/height.phe", header = T)

setwd("/path/to/working_directory/height/")

imMVP <- MVP(
  phe=height, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##g. Run the GWAS for number of panicles

n_panicles <- read.table("/path/to/n_panicles.phe", header = T)

setwd("/path/to/working_directory/n_panicles/")

imMVP <- MVP(
  phe=n_panicles, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

##h. Run the GWAS for panicle length

panicle_length <- read.table("/path/to/panicle_length.phe", header = T)

setwd("/path/to/working_directory/panicle_length/")

imMVP <- MVP(
  phe=panicle_length, #NA is acceptable in phenotype
  geno=genotype,
  map=map,
  K=kinship, #if you have pre-computed GRM, please keep there open, otherwise rMVP will compute it automatically
  nPC.GLM=5,
  maxLine=10000, #smaller value would reduce the memory cost
  threshold=0.05,
  method=c("GLM"),
  file.output=c("pmap", "pmap.signal", "plot", "log")
)

#6. Clean the environment

rm(list = ls())
gc()
