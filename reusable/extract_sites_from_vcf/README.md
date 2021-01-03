## Scripts to pull a subset of specific sites from 1000 Genomes vcf.gz files. 

1. create a tab-delimited file of chr and pos for sites you want, make sure they're human-readable sorted for bcftools
   - if you have rsID instead of chr/pos, you can extract the chr/pos directly from the 1KG files using get_chrpos_from_rsID.sh (input: your_snps.list, output: chrpos.tsv)

2. extract sites from 1000 Genomes vcf.gz files, bcf_extract.sh extracts from all chromosome files at once and creates one file per chromosome
   - if your sites are on just one chromosome, remove `#SBATCH --array=1-22` line from the .sh script and change the filename to the appropriate chromosome
  
3. combine chromosome files from previous step with combine_chromfiles.sh
