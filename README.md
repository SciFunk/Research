Data:
1. all SNPs <1% Afr >5% other pop: sample info, rsID, chr, pos, allele - IDs/allchr_allSamples_rsID_pos_pop_CLEAN.csv
2. all SNPs <1% Afr >5% other pop: unique rsID (per pop) only: sample info, rsID, chr, pos, allele, pop, superpop - IDs/unique_allchr_allSamples_rsID_pos_pop_CLEAN.csv
3. all SNPs <1% Afr >5% other pop: rsID, chr, pos, allelle, pop, superpop - /IDs/allchrPops_rsID.csv
4. all SNPs <1% Afr >5% other pop: unique rsIDs only: rsID, chr, pos, allelle, pop, superpop - /IDs/unique_allchrPops_rsID


More info below:



Necessary files for the next two scripts are in the scripts folder (1000.genome.header.txt, pop_locations.txt)

* To create a file with all samples, rsID, chr, pos, allele (<1% in Afr pops, >5% in other pop): IDs_perSamplewithPos.py (submission file: IDs.qsub, runtime ~few hours)
  * only counts bi-allelic sites, checks both 0's and 1's to see if either are <1% in Afr (then carries through to count 0's or 1's for >5% in other pops)
  * creates a separate file for each chr, combine with combineFiles.py (adds a column with filename title)
  * combined file messy with extra characters, use cleanSamples.py to clean it up 
      * file has columns: chr, pos, rsID, sample, allele
      * *stored as Working/IDs/allchrSampleswithPos_CLEAN.csv*
      * *file with pop/superpop info and only unique rsIDs (see below): Working/IDs/unique_allchr_allSamples_rsID_pos_pop_CLEAN.csv*
        * 559,065 rows in unique file, 12.6 million rows in non-unique (*allchr_allSamples_rsID_pos_pop_CLEAN.csv*)
      
* To create a file with only rsIDs, pop, allele (<1% in Afr pops, >5% in other pop): IDs_perPop.py (IDs.qsub works for this one too)
  * creates a separate file for each chr, once combined use `df.drop_duplicates(subset = 'rs', keep = False)` to remove non-unique SNPs
  * file has columns: rsID, pop, allele
  * *stored as Working/IDs/allchrPops_rsID.csv*
  
* To find significant tissue expression using GTEx data:
