Necessary files for the next two scripts are in the scripts folder (1000.genome.header.txt, pop_locations.txt)

* To create a file with all samples, rsID, chr, pos, allele (<1% in Afr pops, >5% in other pop): IDs_perSamplewithPos.py (submission file: IDs.qsub, runtime ~few hours)
  * creates a separate file for each chr, combine with combineFiles.py (adds a column with filename title)
  * combined file messy with extra characters, use cleanSamples.py to clean it up 
      * file has columns: chr, pos, rsID, sample, allele
      * *stored as Working/IDs/allchrSampleswithPos_CLEAN.csv*
      
* To create a file with only rsIDs, pop, allele (<1% in Afr pops, >5% in other pop): IDs_perPop.py (IDs.qsub works for this one too)
