#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=90G
#SBATCH --array=1-22

snp_list_file=perlegen_classa_snp.list
output_chrpos_filename=perlegen_classa_chrpos.tsv
path_to_1kg_files=/users/afunk2/data/data/modern_genomes/1000genomes/ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502

module load bcftools

bcftools view --include ID==@snp_list_file \
$path_to_1kg_files/ALL.chr${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
| cut -f 1-2 >> output_chrpos_filename

### remove comment lines after
# sed '/^#/d' < perlegen_classa_chrpos.tsv > perlegen_classa_chrpos_nocomments.tsv

### sort chr/pos file for bcftools (human number sort)
# sort -k1,1n -k2,2n perlegen_classa_chrpos_nocomments.tsv > perlegen_classa_chrpos_nocomments_sorted.tsv