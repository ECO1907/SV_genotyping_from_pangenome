###This script is used to merge the VCF files obtained from SV calling, and process the merged VCF

#1. Merge all the VCF files

module load samtools
bcftools merge -m all --threads 8 *.vcf.gz -Oz -o japonica_str_var_merged.vcf.gz

#2. Remove multiallelic SVs

##a. Remove all the SVs with more than 2 alleles

bcftools view -M2 --threads 8 -Oz -o japonica_str_var_biallelic_temp.vcf.gz japonica_str_var_merged.vcf.gz

##b. Remove the variants with duplicated ID (these are considered as multiallelic)

bcftools query -f '%ID\n' japonica_str_var_biallelic_temp.vcf.gz | sort | uniq -d > multiallelic_var.txt
bcftools view -e ID=@multiallele_var.txt -Oz -o japonica_str_var_biallelic.vcf.gz japonica_biallelic_temp.vcf.gz

#4. Obtain the stats of the biallelic VCF and remove the intermediate files

bcftools stats --thread 8 japonica_str_var_biallelic.vcf > japonica_str_var_biallelic_stats.txt

rm japonica_str_var_merged.vcf.gz japonica_str_var_biallelic_temp.vcf.gz multiallelic_var.txt

#5. Annotate the merged VCF using the information of the pangenome VCF

##a. Obtain the desired information from the pangenome VCF (TEs annotation) in TSV files and compress them

bcftools query -f '%CHROM\t%POS\t%SVTYPE\t%TE_ANNOT\n' /path/to/pangenome.vcf > pangenome_info.tsv
bgzip -c pangenome_info.tsv > pangenome_info.tsv.gz
tabix -s 1 -b 2 -e 2 pangenome_info.tsv.gz

bcftools query -f '%CHROM\t%POS\t%END\n' /path/to/pangenome.vcf > end_info.tsv
bgzip -c end_info.txt > end_info.txt.gz
tabix -s 1 -b 2 -e 2 end_info.txt.gz

##b. Create a VCF header that descripts the annotation you're adding to the new VCF

touch new_headers.txt
echo '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' > new_headers.txt
echo '##INFO=<ID=TE_ANNOT,Number=1,Type=String,Description="Transposable element annotation">' >> new_headers.txt

touch end_header.txt
echo '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">' > end_header.txt

##c. Annotate the SVs in the new VCF using the created files by comparing the CHROM and POS information

bcftools annotate --threads 8 -h new_headers.txt -a pangenome_info.tsv.gz -c CHROM,POS,SVTYPE,TE_ANNOT -Oz -o japonica_str_var_annot_temp1.vcf.gz japonica_str_var_biallelic.vcf.gz
bcftools annotate --threads 8 -h end_header.txt -a end_info.txt.gz -c CHROM,POS,INFO/END -Oz -o japonica_str_var_annot_temp2.vcf.gz japonica_str_var_annot_temp1.vcf.gz

#6. Replace the IDs for the next format SVTYPE:CHROM:POS..END

##a. Obtain the IDs components and order them

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\n' japonica_str_var_annot_temp2.vcf.gz > id_components.txt
awk 'BEGIN{OFS="\t"} {print $4 ":" $1 ":" $2 ".." $3}' id_components.txt > new_ids.txt

##b. Create a file that contains the genomic positions with the new IDs

paste <(bcftools query -f '%CHROM\t%POS\n' japonica_str_var_annot_temp2.vcf.gz) new_ids.txt > id_map.txt

##c. Set the new IDs

bcftools annotate -a id_map.txt -c CHROM,POS,ID --threads 8 -Ov -o japonica_str_var_annot.vcf japonica_str_var_annot_temp2.vcf.gz

#7. Change the heterozygotes for homozygotes. This is done for practical purposes. Then compress and index the final VCF, and remove intermediate files

sed -i -E 's/(^|[^0-9])0\/1([^0-9]|$)/\11\/1\2/g; s/(^|[^0-9])1\/0([^0-9]|$)/\11\/1\2/g' japonica_str_var_annot.vcf

bgzip -c japonica_str_var_annot.vcf > str_var_japonica.vcf.gz
tabix -p vcf str_var_japonica.vcf.gz

rm japonica_str_var_biallelic.vcf.gz japonica_str_var_annot.vcf japonica_str_var_annot_temp2.vcf.gz japonica_str_var_annot_temp1.vcf.gz pangenome_info.tsv* end_info.tsv* new_headers.txt end_header.txt id_components.txt new_ids.txt id_map.txt

#8. Filter out the rare (MAF < 0.03) and missing data variants (< 10%) and recode it to a 012 file

##a. Filter the variants

conda activate plink
plink --vcf str_var_japonica.vcf.gz --maf 0.03 --geno 0.1 --make-bed --out str_var_japonica

##b. Recode to obtain the 012 file

plink --bfile str_var_japonica --recode A --out str_var_japonica_recoded
