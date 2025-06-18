#!/bin/bash
#SBATCH --job-name="vacall"
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --partition=all
#SBATCH --time=24:00:00
#SBATCH --output=%x_%A_%a.log
#SBATCH --error=%x_%A_%a.err

###This script is used to map short reads to a pangenome and perform SV calling### 

#1. Indicate the folder where are the reads

folder=$1
folder_name=$(basename "$folder")

#2. Create a new working directory

##a. Create the directory and change to it

mkdir -p $folder_name
cd $folder_name

##b. Copy the pangenome files to the working directory

cp /path/to/ricepangenome.vg . #Pangenome file
cp /path/to/ricepangenome.xg . #Pangenome index file

#3. Trim the reads

##a. Concatanate the reads if necessary

cat $folder/*_1.fq.gz > ${folder_name}_notrimmed_1.fq.gz
cat $folder/*_2.fq.gz > ${folder_name}_notrimmed_2.fq.gz

##b. Trim the reads

module load conda
conda activate bbmap
bbduk.sh in=${folder_name}_notrimmed_1.fq.gz in2=${folder_name}_notrimmed_2.fq.gz out=${folder_name}_trimmed_1.fq.gz out2=${folder_name}_trimmed_2.fq.gz ref=/path/to/adapters.fasta ktrim=r k=23 mink=11 hdist=1 tpe tbo

#4. Find the variants from the pangenome in the reads

##a. Map the reads to the pangenome

conda activate vg
vg giraffe -t 8 -x ricepangenome.xg -f ${folder_name}_trimmed_1.fq.gz -f ${folder_name}_trimmed_2.fq.gz > ${folder_name}_pang_mapping.gam

##b. Obtain the mapping stats

vg stats -a ${folder_name}_pang_mapping.gam > ${folder_name}_mapping_stats.txt

##c. Pack the gam file

vg pack -x ricepangenome.xg -g ${folder_name}_pang_mapping.gam -o ${folder_name}_pang_mapping.pack -t 8

##d. Call the variant

vg call -t 8 -a -k ${folder_name}_pang_mapping.pack ricepangenome.xg > ${folder_name}_variants_no_name.vcf

#5. Process the VCF file

##a. Change the sample name

module load samtools
bcftools reheader -s <(echo "${folder_name}") -o ${folder_name}_variants.vcf ${folder_name}_variants_no_name.vcf

##b. Obtain the VCF stats

bcftools stats ${folder_name}_variants.vcf > ${folder_name}_variants_stats.txt

##c. Compress and index the VCF

bgzip -c ${folder_name}_variants.vcf > ${folder_name}_variants_comp.vcf.gz
tabix -p vcf ${folder_name}_variants_comp.vcf.gz

#6. Cleaning the working directory

##a. Remove every intermediate file

rm ${folder_name}_notrimmed_1.fq.gz ${folder_name}_notrimmed_2.fq.gz ${folder_name}_trimmed_1.fq.gz ${folder_name}_trimmed_2.fq.gz ${folder_name}_pang_mapping.gam ${folder_name}_pang_mapping.pack ${folder_name}_variants_no_name.vcf ricepangenome.*

##b. Move the desired files to a new directory

mv ${folder_name}_variants_comp.vcf.gz* /path/to/new_directory

##c. Return to the working directory

cd ..
