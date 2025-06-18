###This script is used to retain the japonica accessions from the SNPs VCF and process it###

#1. Filter out the samples that are not japonica from the SNPs VCF

module load samtools
bcftools view -S /path/to/japonica_accessions.txt --threads 8 -Oz -o snps_japonica_temp.vcf.gz /path/to/rice_snps_pruned.vcf.gz

#2. Change the IDs for the format SNP:CHROM:POS

##a. Obtain the IDs components and order them

bcftools query -f '%CHROM\t%POS\n' snps_japonica_temp.vcf.gz > id_components.txt
awk 'BEGIN{OFS="\t"} {print "SNP" ":" $1 ":" $2}' id_components.txt > new_ids.txt

##b. Create a file that contains the genomic positions with the new IDs

paste <(bcftools query -f '%CHROM\t%POS\n' snps_japonica_temp.vcf.gz) new_ids.txt > id_map.txt

##c. Set the new IDs

bcftools annotate -a id_map.txt -c CHROM,POS,ID --threads 8 -Ov -o snps_japonica.vcf snps_japonica_temp.vcf.gz

#3. Change the heterozygotes for homozygotes. This is done for practical purposes. Then compress and index the final VCF, and remove intermediate files

sed -i -E 's/(^|[^0-9])0\/1([^0-9]|$)/\11\/1\2/g' snps_japonica.vcf

bgzip -c snps_japonica.vcf > snps_japonica.vcf.gz
tabix -p vcf snps_japonica.vcf.gz

#4. Clean the working directory

rm snps_japonica.vcf id_map.txt new_ids.txt id_components.txt snps_japonica_temp.vcf.gz

#5. Filter out the rare (MAF < 0.03) and missing data variants (< 10%) and recode it to a 012 file

##a. Filter the variants

conda activate plink
plink --vcf snps_japonica.vcf.gz --maf 0.03 --geno 0.1 --make-bed --out snps_japonica

##b. Recode to obtain the 012 file

plink --bfile snps_japonica --recode A --out snps_japonica_recoded
