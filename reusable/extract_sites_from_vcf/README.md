## Scripts pull a subset of specific sites from 1000 Genomes vcf.gz files. 

1. get_chrpos_from_rsID.sh: create a tab-delimited file of chr and pos for sites you want, make sure they're human-readable sorted for bcftools
- if you have rsID instead of chr/pos, you can extract the chr/pos directly from the 1KG files using (input: your_snps.list, output: chrpos.tsv)
