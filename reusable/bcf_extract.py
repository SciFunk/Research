##### using a combination of python and bcftools to pull specific sites (lines) from the 1KG data

### 1. create a tab-delimited file of chr and pos for sites you want, make sure they're human-readable sorted for bcftools
### if you have rsID instead of chr/pos, you can extract the chr/pos directly from the 1KG files using (input: your_snps.list, output: chrpos.tsv)

#!/bin/sh
#SBATCH --array=1-22

module load bcftools/1.9
bcftools view --include ID==@your_snps.list /path/to/1000genomes/release/20130502/ALL.chr${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | cut -f 1-2 >> chrpos.tsv

### remove comment lines after
sed '/^#/d' < chrpos.tsv > chrpos_nocomments.tsv

### sort chr/pos file for bcftools (human number sort)
sort -k1,1n -k2,2n chrpos_nocomments.tsv > chrpos_sorted.tsv

###
