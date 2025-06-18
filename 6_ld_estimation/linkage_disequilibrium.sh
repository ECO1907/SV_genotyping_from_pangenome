###This script is used to calculate the linkage disequilibrium between variants###

#1 Keep only the japonica samples from the pre-pruned SNPs VCF

module load samtools
bcftools view -S /path/to/japonica_accessions.txt --threads 8 -Oz -o prepruned_snps_japonica.vcf.gz /path/to/prepruned_snps_rice.vcf

#2 Obtain the binary files from pre-pruned japonica SNPs VCF

conda activate plink
plink --vcf prepruned_snps_japonica.vcf.gz --maf 0.03 --geno 0.1 --thin 0.1 --make-bed --out prepruned_snps_japonica

#3 Merge the binary files of the pre-pruned SNPs and SVs datasets

plink --bfile /path/to/str_var_japonica --bmerge prepruned_snps_japonca --make-bed prepruned_both_japonica

#4 Calculate the linkage disequilibrium

plink --bfile prepruned_both_japonica --allow-extra-chr --r2 --ld-window 100 --ld-window-kb 400 --ld-window-r2 0 --out both_japonica_ld

#5 Split the LD file in the 3 scenarios SVs-SVs, SNPs-SNPs and SVs-SNPs

##a. Split the LD file in the 3 scenarios

awk '$3 ~ /^SNP/ && $6 ~ /^SNP/' both_japonica_ld.ld > SNP_SNP.ld
awk '($3 ~ /^SNP/ && $6 !~ /^SNP/) || ($3 !~ /^SNP/ && $6 ~ /^SNP/)' both_japonica_ld.ld > SNP_SV.ld
awk '$3 !~ /^SNP/ && $6 !~ /^SNP/' both_japonica_ld.ld > SV_SV.ld

##b. Add the headers to SV_SV.ld and SNP_SV.ld

(head -n 1 SNP_SNP.ld && cat SV_SV.ld) > temp1.ld && mv temp1.ld SV_SV.ld
(head -n 1 SNP_SNP.ld && cat SNP_SV.ld) > temp2.ld && mv temp2.ld SNP_SV.ld

#6 Remove the temporal files

rm temp1.ld temp2.ld
