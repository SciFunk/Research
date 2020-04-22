### count all consequences across all chromosomes
"""""""""""""""""""""""""""""""""""""""""""""""""
### .sh file used to submit script on hpc:
#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --mem=100G
#SBATCH --array=1-22

python count_annotation.py /path/to/files/ALL.chr${SLURM_ARRAY_TASK_ID}.phase3_shapeit2_mvncall_integrated_v5_func_anno.20130502.sites.vcf.gz
"""""""""""""""""""""""""""""""""""""""""""""""""

import pandas as pd
import sys

### this file is just the original 1KG header, it's needed so each chunk of vcf has a header
with open("1000.genome.header.txt", "r") as header:
    for line in header:
        header_line = line.strip().split()
        
anno_iter = pd.read_table(sys.argv[1], compression='gzip', skiprows=246, sep='\t', header=None, names=header_line, 
                   engine='python', chunksize=50000)
for anno in anno_iter:
    anno['vep'] = anno['INFO'].str.split("|", expand = True)[4]
    anno['vep'].value_counts().reset_index().to_csv('counts_raw.csv', header=False, index=False, mode='a')


### this should be run in a different script, once all 22 chromosomes have completed in the above job
### clean above file and get counts for each category
# counts = pd.read_csv('counts_raw.csv', header=None, names=['conseq', 'counts'])
# final = pd.melt(counts, id_vars = ['counts'], value_name = 'consequence', value_vars = ['conseq'])
# final = final.drop('variable', axis=1).groupby('consequence').sum()
# final.sort_values(['counts'], ascending=False).reset_index().to_csv('counts_final.csv', index=False)
