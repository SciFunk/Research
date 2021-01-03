#!/bin/sh
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=4:00:00
#SBATCH --mem=90G
#SBATCH --array=1-22

chrpos_sorted_filename=perlegen_classa_chrpos_sorted.tsv
output_file_suffix=test
path_to_1kg_files=/users/afunk2/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502

module load bcftools/1.9

bcftools view $path_to_1kg_files/ALL.chr${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
-T $chrpos_sorted_filename > chr${SLURM_ARRAY_TASK_ID}-$output_file_suffix.vcf