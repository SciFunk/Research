#!/bin/sh
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --time=4:00:00
#SBATCH --mem=90G

output_file_suffix=test

module load bcftools/1.9

ls -v chr* > files.list
bcftools concat -f files.list -o allchr_$output_file_suffix.vcf.gz -O z
