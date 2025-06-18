#!/bin/bash
#SBATCH --job-name="vc_master"
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=all
#SBATCH --time=01:00:00
#SBATCH --output=vc_master.log
#SBATCH --error=vc_master.err

###This script is used to call the variant_calling.sh script to do multiple parallel mapping and variant calling###

#1. Set the directory that contains the directories with the short reads (each one for each rice accession)

PARENT_DIR=/path/to/directory

#2. Call the variant_calling.sh for each accession directory

for folder in "$PARENT_DIR"/*; do
	folder_name=$(basename "$folder")
	sbatch --output=${folder_name}.log --error=${folder_name}.err variant_calling.sh "$folder"
done
