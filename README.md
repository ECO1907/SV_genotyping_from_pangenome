# SV_genotyping_from_pangenome
This script contains the scripts used to detect SVs in a japonica rice population using a graph-based pangenome and short reads, run genome-wide association studies (GWAS) with GLM and BLINK methods, and perform genomic prediction. The results obtained with the SVs were also compared to the ones obtained with SNPs and the combination of SVs and SNPs. All of them were employed in the project "GENOMIC STRUCTURAL VARIANTS, A PROMISING TOOL FOR RICE GENOMIC-ASSISTED BREEDING", which is divided in several steps:
## 1. Pangenome SVs detection in the japonica rice population
Short reads mapping to the pangenome and SVs call for each rice accession.
## 2. SVs merging and filtering
Merging of the variant calling format (VCF) files containing the SVs and filtering of the merged VCF.
## 3. SNPs VCF processing
Filtering of the VCF containing the SNPs.
## 4. SVs+SNPs VCF preparation
Concatenation of the both VCFs (SVs and SNPs) into a single one (SVs+SNPs).
## 5. PCA and Kinship
Principal Component Analyses (PCA) and kinship relationship analyses for the three datasets (SVs, SNPs and SVs+SNPs).
## 6. LD estimation
Linkage Disequilibrium (LD) estimation of three scenarios (SV-SV, SNP-SNP and SV-SNP).
## 7. GWAS
Genome-Wide Association Studies (GWAS) using the three datasets to find associations with phenotypic data for eight traits (grain width, grain length, grains per panicle, panicle length, number of panicles, plant height, heading time and allelopathy index) with GLM and BLINK methods.
## 8. Heritability calculation
Calculation of phenotypic variability explained by the genetic variability of the three datasets for each trait using.
## 9. Genomic Prediction accuracy assessment
Calculation of genomic prediction accuracy with the three datasets for each trait using a model that captures non-additive effects.
