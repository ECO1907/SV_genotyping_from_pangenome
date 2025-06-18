# SV_genotyping_from_pangenome
This script contains the scripts used to detect SVs in a japonica rice population using a graph-based pangenome and short reads, run genome-wide association studies (GWAS) with GLM and BLINK methods, and perform genomic prediction. The results obtained with the SVs were also compared to the ones obtained with SNPs and the combination of SVs and SNPs. All of them were employed in the project "GENOMIC STRUCTURAL VARIANTS, A PROMISING TOOL FOR RICE GENOMIC-ASSISTED BREEDING", which is divided in several steps:
1. Short reads mapping to the pangenome and SVs call for each rice accession
   
3. Merging of the variant calling format (VCF) files containing the SVs and filtering of the merged VCF
4. Filtering of the VCF containing the SNPs
5. Concatenation of the both VCFs (SVs and SNPs) into a single one
6. Principal Component Analyses (PCA) and kinship relationship analyses for the three datasets
7. Linkage Disequilibrium (LD) estimation
8. Genome-Wide Association Studies
9. Heritability calculation
10. Genomic Prediction accuracy assessment
