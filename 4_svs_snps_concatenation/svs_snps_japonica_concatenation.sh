###This script is used to concatenate the SV and SNPs (both) datasets and process it###

#1. Concantenate the SVs and SNPs VCF

module load samtools
bcftools concat -a --threads 8 -Oz -o both_japonica.vcf.gz /path/to/str_var_japonica.vcf.gz /path/to/snps_japonica.vcf.gz

#2. Filter out the rare (MAF < 0.03) and missing data variants (< 10%) and recode it to a 012 file

##a. Filter the variants

conda activate plink
plink --vcf both_japonica.vcf.gz --maf 0.03 --geno 0.1 --make-bed --out both_japonica

##b. Recode to obtain the 012 file

plink --bfile both_japonica --recode A --out both_japonica_recoded
